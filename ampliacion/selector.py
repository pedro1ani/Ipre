"""Bloque 5: elección por metabolito dentro del entrenamiento, prueba posterior por grupos."""
import json

import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import nan_euclidean_distances
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from .diseno_5_6 import DESTINO, PANELES, TIPOS_PRUEBA, diseno, huella, panel_elegible
from .incertidumbre import perfiles
from .medir_bloques import medir
from .modelos import CONFIG, predecir
from .preparar import SALIDA
from .validar import ARCHIVO_MASCARAS

METODOS = list(CONFIG)
MIN_CELDAS, MIN_GRUPOS, MIN_BLOQUES = 10, 5, 2
EMPATE_MM = 1e-12


def decidir(p, objetivos):
    """No acepta predicciones externas. Los patrones repetidos no multiplican el soporte."""
    indice = ["interno", "grupo", "metabolito", "celda", "evento"]
    pares = p.assign(error=(p.pred_mM-p.verdad_mM).abs()).pivot(
        index=indice, columns="metodo", values="error").reindex(columns=METODOS).dropna()
    celdas = pares.groupby(level=indice[:-1]).mean()
    filas = []
    for m in objetivos:
        t = celdas[celdas.index.get_level_values("metabolito") == m]
        grupos = t.index.get_level_values("grupo").nunique()
        por_bloque = t.groupby(level="interno").mean()
        mae = t.mean()
        ganador, victorias = None, 0
        if len(t) < MIN_CELDAS or grupos < MIN_GRUPOS or len(por_bloque) < MIN_BLOQUES:
            razon = "soporte_interno_insuficiente"
        else:
            orden = mae.sort_values()
            if orden.iloc[1] - orden.iloc[0] <= EMPATE_MM:
                razon = "empate"
            else:
                mejor = orden.index[0]
                otros = por_bloque.drop(columns=mejor)
                victorias = int((otros.min(axis=1) - por_bloque[mejor] > EMPATE_MM).sum())
                if victorias >= MIN_BLOQUES:
                    ganador, razon = mejor, "ventaja_repetida"
                else:
                    razon = "ventaja_no_repetida"
        fila = dict(metabolito=m, elegido=ganador, razon=razon, celdas=len(t), grupos=grupos,
                    internos=len(por_bloque), victorias=victorias)
        for metodo in METODOS:
            fila[f"MAE_{metodo}_mM"] = mae[metodo]
            for i in range(3):
                fila[f"MAE_{metodo}_interno_{i}"] = por_bloque[metodo].get(i, np.nan)
        filas.append(fila)
    return pd.DataFrame(filas)


def diagnosticar(x, ish, q, cols, objetivos, grupos, k):
    """Descripciones del soporte, nunca variables para decidir el ganador."""
    datos = []
    pares = x[cols].notna().astype(int)
    n_pares = pares.T @ pares
    correlaciones = x[cols].corr(min_periods=5)
    for m in objetivos:
        s = x[m].dropna()
        row = dict(metabolito=m, n_train=len(x), medidos=len(s), grupos_medidos=grupos.loc[s.index].nunique(),
                   fraccion_faltante=1-len(s)/len(x), IQR_train_mM=s.quantile(.75)-s.quantile(.25))
        if m in cols:
            r = correlaciones[m].drop(m).dropna()
            if not r.empty:
                for nombre, valor in r.items():
                    datos.append(dict(**row, predictor=nombre, pares=int(n_pares.at[m, nombre]),
                                      r=valor, mas_fuerte=nombre == r.abs().idxmax()))
                continue
        datos.append(dict(**row, predictor=None, pares=0, r=np.nan, mas_fuerte=False))
    z = StandardScaler().fit(x[cols].to_numpy())
    consulta, i, j = perfiles(ish, q, cols)
    train = x[cols].to_numpy()
    dist = nan_euclidean_distances(z.transform(consulta), z.transform(train))
    vecinos = []
    for n, (f, m) in enumerate(zip(i, j)):
        if m < 0:
            continue
        validos = np.flatnonzero(np.isfinite(train[:, m]) & np.isfinite(dist[f]))
        cercanos = validos[np.argsort(dist[f, validos], kind="stable")[:k]]
        if len(cercanos):
            comunes = (np.isfinite(train[cercanos]) & np.isfinite(consulta[f])).sum(axis=1)
            vecinos.append(dict(tipo=q.tipo.iat[n], evento=q.evento.iat[n], metabolito=q.metabolito.iat[n],
                                donantes_disponibles=len(validos), vecinos_examinados=len(cercanos),
                                distancia_media=dist[f, cercanos].mean(), compartidos_medios=comunes.mean()))
    return pd.DataFrame(datos), pd.DataFrame(vecinos)


def aplicar(fijos, decisiones):
    """Entrega la elección si responde; en su defecto una mediana explícita, o ninguna predicción."""
    elegido = decisiones.set_index("metabolito").elegido
    med = fijos[fijos.metodo.eq("Mediana")].copy().set_index("evento")
    med["eleccion"] = med.metabolito.map(elegido)
    med["politica"] = "mediana_general_sin_eleccion"
    campos = [c for c in ["pred_mM", "motivo", "n_predictores", "n_donantes"] if c in med]
    for metodo in METODOS:
        datos = fijos[fijos.metodo.eq(metodo)].set_index("evento")
        sel = med.eleccion.eq(metodo)
        si = sel & datos.pred_mM.notna()
        med.loc[si, campos] = datos.loc[si, campos]
        med.loc[si, "politica"] = "eleccion_interna"
        med.loc[sel & ~si, "politica"] = "mediana_general_por_falta_de_soporte_externo"
    med.loc[med.pred_mM.isna(), "politica"] = "abstencion"
    return med.reset_index().assign(metodo="Selector")


def ejecutar(destino=DESTINO):
    ish, mmc, grupos, ref, opciones, tramos = diseno(destino)
    detalle = destino / "detalle"
    resultados, decisiones, perfiles_train, vecinos, interno_total = [], [], [], [], []
    for b, ((train, qext), internas) in enumerate(tramos):
        print(f"Selector: bloque {b+1}/5", flush=True)
        iqr = train.quantile(.75) - train.quantile(.25)
        for panel, estudio in PANELES.items():
            for esc in ("A", "B"):
                etiquetas = dict(panel=panel, escenario=esc, bloque=b)
                opciones_fold = {m: int(opciones.loc[estudio, b, esc, m]) for m in METODOS}
                def predecir_particion(sub, q, interno):
                    cols = panel_elegible(sub, grupos, ref, panel)
                    x = sub if esc == "A" else pd.concat([sub, mmc])
                    if set(grupos.loc[x.index]) & set(grupos.loc[q.muestra]):
                        raise ValueError("Un grupo de evaluación aparece en entrenamiento.")
                    k = CONFIG["KNN"][opciones_fold["KNN"]]["n_neighbors"]
                    diagnostico, vec = diagnosticar(x, ish, q, cols, ref.index, grupos, k)
                    perfiles_train.append(diagnostico.assign(**etiquetas, interno=interno))
                    vecinos.append(vec.assign(**etiquetas, interno=interno))
                    partes = []
                    for metodo in METODOS:
                        p, _ = predecir(x, ish, q, cols, metodo, opciones_fold[metodo])
                        partes.append(p.assign(metodo=metodo, grupo=p.muestra.map(grupos), interno=interno))
                    return pd.concat(partes, ignore_index=True)
                internos = pd.concat([predecir_particion(sub, q, i)
                                      for i, (sub, q) in enumerate(internas)], ignore_index=True)
                interno_total.append(internos.assign(**etiquetas))
                # Congelar TODAS las decisiones antes de predecir o leer errores externos.
                reglas = {tipo: decidir(internos[internos.tipo.eq(tipo)], ref.index) for tipo in TIPOS_PRUEBA}
                externos = predecir_particion(train, qext, -1)
                for tipo, regla in reglas.items():
                    fijos = externos[externos.tipo.eq(tipo)].copy()
                    decisiones.append(regla.assign(**etiquetas, tipo=tipo))
                    if fijos.empty:
                        continue
                    operativo = aplicar(fijos, regla)
                    conjunto = pd.concat([fijos.assign(politica="metodo_fijo"), operativo], ignore_index=True)
                    t = conjunto.assign(**etiquetas, experimento=f"{panel}_{tipo}", alternativa="misma_regla",
                                        rep=0, IQR_train=conjunto.metabolito.map(iqr),
                                        referencia_mM=conjunto.metabolito.map(ref))
                    resultados.append(t)
    pred = pd.concat(resultados, ignore_index=True)
    pred.to_csv(detalle / "selector_predicciones.csv.gz", index=False)
    pd.concat(interno_total).to_csv(detalle / "selector_predicciones_internas.csv.gz", index=False)
    reglas = pd.concat(decisiones, ignore_index=True)
    reglas.to_csv(detalle / "decisiones_internas.csv.gz", index=False)
    perfiles = pd.concat(perfiles_train, ignore_index=True)
    perfiles.to_csv(detalle / "correlaciones_entrenamiento.csv.gz", index=False)
    pd.concat(vecinos).to_csv(detalle / "soporte_vecinos.csv.gz", index=False)
    comun = perfiles[perfiles.interno.ge(0) & perfiles.predictor.notna()].copy()
    comun["signo"] = np.sign(comun.r)
    estabilidad = comun.groupby(["panel","escenario","bloque","metabolito","predictor"]).agg(
        internos=("r","size"), pares_min=("pares","min"), r_min=("r","min"), r_max=("r","max"),
        signos=("signo","nunique"), veces_mas_fuerte=("mas_fuerte","sum")).reset_index()
    estabilidad.to_csv(detalle / "estabilidad_correlaciones.csv.gz", index=False)
    glob, met, por_grupo = medir(pred)
    glob.to_csv(destino / "selector_global.csv", index=False)
    met.to_csv(detalle / "selector_por_metabolito.csv.gz", index=False)
    por_grupo.to_csv(detalle / "selector_por_grupo.csv.gz", index=False)
    claves = ["panel","tipo","escenario"]
    uso = pred[pred.metodo.eq("Selector")].groupby(claves + ["politica"], dropna=False).agg(
        eventos=("evento","size"), celdas=("celda","nunique"), grupos=("grupo","nunique")).reset_index()
    uso.to_csv(destino / "uso_selector.csv", index=False)
    reglas.assign(elegido=reglas.elegido.fillna("ninguno")).groupby(
        claves + ["elegido","razon"]).size().rename("decisiones").reset_index().to_csv(
            destino / "resumen_decisiones.csv", index=False)
    (destino / "selector.json").write_text(json.dumps(dict(
        celdas_minimas=MIN_CELDAS, grupos_minimos=MIN_GRUPOS, bloques_minimos=MIN_BLOQUES,
        empate_mM=EMPATE_MM, sha256_individual=huella(ARCHIVO_MASCARAS),
        sha256_conjunta=huella(detalle / "mascaras_conjuntas.json"),
        sha256_opciones=huella(SALIDA / "ajustes.csv")), indent=2) + "\n")


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        ejecutar()
