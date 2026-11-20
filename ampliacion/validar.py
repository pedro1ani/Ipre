"""Evaluación anidada en los bloques archivados. Ejecutar: python -m ampliacion.validar."""
import hashlib
import json
import platform
from importlib.metadata import version

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from .modelos import CONFIG, SEMILLA, diagnostico_extremos, elegibilidad, predecir
from .preparar import RAIZ, SALIDA, preparar

ARCHIVO_MASCARAS = RAIZ / "evaluacion" / "mascaras.json"


def bloques(registros, cohorte, grupos, objetivos, n_bloques):
    """Comprueba verdades y partición completa antes de separar entrenamiento y pruebas."""
    q = pd.DataFrame(registros)
    conocidos = {(s, m) for m in objetivos for s in cohorte.index[cohorte[m].notna()]}
    assert set(zip(q.muestra, q.metabolito)) == conocidos and q.celda.is_unique
    assert (q.celda == q.muestra + "|" + q.metabolito).all()
    np.testing.assert_allclose(q.verdad_mM, [cohorte.at[s, m] for s, m in zip(q.muestra, q.metabolito)],
                               atol=1e-14, rtol=1e-12)
    assert set(q.bloque) == set(range(n_bloques)) and set(q.muestra) == set(cohorte.index)
    assert q.groupby(q.muestra.map(grupos)).bloque.nunique().eq(1).all()
    for bloque, test in q.groupby("bloque", sort=True):
        excluidos = grupos.loc[test.muestra].unique()
        train = cohorte.loc[~grupos.loc[cohorte.index].isin(excluidos)]
        assert set(grupos.loc[train.index]).isdisjoint(excluidos)
        yield bloque, train, test


def elegir(train, registros, mmc, grupos, ref, ampliado, etiqueta, busqueda):
    """Una elección con MAE en mM y otra histórica con MAE/R; las verdades son internas."""
    panel = train.columns if ampliado else ref.index
    partes = []
    for interno, sub, q in bloques(registros, train, grupos, ref.index, 3):
        e = elegibilidad(sub[panel], grupos)
        columnas = e.index[e.elegible]
        for escenario, x in [("A", sub), ("B", pd.concat([sub, mmc]))]:
            for metodo in ("KNN", "Iterative"):
                for k in range(len(CONFIG[metodo])):
                    p, info = predecir(x, train, q, columnas, metodo, k)
                    partes.append(p[["celda", "metabolito", "verdad_mM", "pred_mM"]].assign(
                        metodo=metodo, escenario=escenario, k=k, interno=interno,
                        aviso=info["aviso_convergencia"]))
    tabla = pd.concat(partes)
    elecciones = {}
    for (escenario, metodo), g in tabla.groupby(["escenario", "metodo"]):
        comunes = g.pivot(index="celda", columns="k", values="pred_mM").notna().all(axis=1)
        g = g[g.celda.isin(comunes.index[comunes])].copy()
        if g.empty:
            raise ValueError(f"Sin celdas comparables en selección: {etiqueta}, {escenario}, {metodo}")
        g["mM"] = (g.pred_mM - g.verdad_mM).abs()
        g["R"] = g.mM / g.metabolito.map(ref)
        puntajes = g.groupby(["k", "metabolito"])[["mM", "R"]].mean().groupby("k").mean()
        for k, valores in puntajes.iterrows():
            busqueda.append(dict(bloque=etiqueta, panel="ampliado" if ampliado else "comun",
                                 escenario=escenario, metodo=metodo, k=k, macro_mM=valores.mM,
                                 macro_R=valores.R, n_comunes=int(comunes.sum()),
                                 n_intentos=len(comunes),
                                 avisos_en_bloques=int(g.loc[g.k.eq(k)].groupby("interno").aviso.any().sum())))
        for criterio in ["mM", "R"]:
            elecciones[escenario, metodo, criterio] = int(puntajes[criterio].idxmin())
    for escenario in ["A", "B"]:
        for criterio in ["mM", "R"]:
            elecciones[escenario, "Mediana", criterio] = 0
    return elecciones


def aporte_individual(ish, mmc, grupos, ref, archivo, pred, ajustes):
    """Añade una sola columna a la vez, sin reabrir la selección de opciones ni de objetivos."""
    base = pred[pred.diseno.eq("seleccion_mM")].set_index(["bloque", "escenario", "metodo", "celda"]).sort_index()
    opciones = ajustes[ajustes.diseno.eq("seleccion_mM")].set_index(["bloque", "escenario", "metodo"]).k
    partes = []
    for b, train, q in bloques(archivo["externas"], ish, grupos, ref.index, 5):
        comun = elegibilidad(train[ref.index], grupos)
        for nueva in ish.columns.difference(ref.index, sort=False):
            soporte = elegibilidad(train[list(ref.index) + [nueva]], grupos)
            columnas = list(comun.index[comun.elegible])
            if soporte.at[nueva, "elegible"]:
                columnas.append(nueva)
            for esc, x in [("A", train), ("B", pd.concat([train, mmc]))]:
                for metodo in ["KNN", "Iterative"]:
                    llave = (b, esc, metodo)
                    p, info = predecir(x, ish, q, columnas, metodo, opciones.loc[llave])
                    p["anterior"] = p.celda.map(base.loc[llave, "pred_mM"])
                    ok = p.pred_mM.notna() & p.anterior.notna()
                    p["MAE_comun"] = (p.anterior - p.verdad_mM).abs().where(ok)
                    p["MAE_con_predictor"] = (p.pred_mM - p.verdad_mM).abs().where(ok)
                    partes.append(p.assign(predictor=nueva, escenario=esc, metodo=metodo,
                                           grupo=p.muestra.map(grupos).where(ok),
                                           activo=nueva in columnas, aviso=info["aviso_convergencia"]))
    claves = ["predictor", "escenario", "metodo"]
    t = pd.concat(partes, ignore_index=True)
    por_metabolito = t.groupby(claves + ["metabolito"]).agg(
        intentos=("celda", "size"), n_pares=("MAE_comun", "count"), n_grupos=("grupo", "nunique"),
        MAE_comun=("MAE_comun", "mean"), MAE_con_predictor=("MAE_con_predictor", "mean"))
    por_metabolito["delta_mM"] = por_metabolito.MAE_con_predictor - por_metabolito.MAE_comun
    resumen = por_metabolito.groupby(level=claves)[["MAE_comun", "MAE_con_predictor", "delta_mM"]].mean()
    resumen["n_pares"] = por_metabolito.groupby(level=claves).n_pares.sum()
    por_bloque = t.groupby(claves + ["bloque"])[["activo", "aviso"]].any()
    resumen["bloques_activos"] = por_bloque.groupby(level=claves).activo.sum()
    resumen["bloques_sin_convergencia"] = por_bloque.groupby(level=claves).aviso.sum()
    por_metabolito.to_csv(SALIDA / "predictor_por_metabolito.csv")
    resumen.to_csv(SALIDA / "predictor_individual.csv")


def ejecutar():
    ish, mmc, wt, meta, ref = preparar()
    grupos = meta.grupo
    assert (meta.loc[mmc.index, "growth_rate"] - 0.2).abs().idxmin() == "MMC2_GALACTOSE"
    assert meta.loc["MMC2_GLUCOSE", "carbono"] == "glucose"
    archivo = json.loads(ARCHIVO_MASCARAS.read_text())
    historicos = pd.read_csv(RAIZ / "evaluacion/parametros_seleccionados.csv")
    historicos = historicos.set_index(["bloque_externo", "escenario", "metodo"]).configuracion
    validacion, ajustes, soporte, extremos, escalas, busqueda = [], [], [], [], [], []
    for b, train, q in bloques(archivo["externas"], ish, grupos, ref.index, 5):
        interior = archivo["internas"][str(b)]
        comun = elegir(train, interior, mmc, grupos, ref, False, b, busqueda)
        amplio = elegir(train, interior, mmc, grupos, ref, True, b, busqueda)
        # El requisito en Ishii garantiza soporte propio. Añadir filas solo puede sumar
        # observaciones, grupos y pares; se comprueba que B nunca reduce ese conjunto.
        paneles = {}
        for nombre, panel in [("comun", ref.index), ("ampliado", ish.columns)]:
            for etapa, sub in [("externa", train)] + [(f"interna_{i}", s) for i, s, _ in
                                                       bloques(interior, train, grupos, ref.index, 3)]:
                for minimo in (8, 10, 12):
                    a = elegibilidad(sub[panel], grupos, minimo)
                    e_b = elegibilidad(pd.concat([sub, mmc])[panel], grupos, minimo)
                    assert (~a.elegible | e_b.elegible).all()
                    soporte.append(a.assign(bloque=b, etapa=etapa, panel=nombre, minimo=minimo,
                                            elegible_B_sin_interseccion=e_b.elegible,
                                            n_B=e_b.n).rename_axis("metabolito").reset_index())
                    if etapa == "externa":
                        if minimo == 10:
                            paneles[nombre] = a.index[a.elegible]
                        elif nombre == "ampliado":
                            paneles[f"minimo_{minimo}"] = a.index[a.elegible]
        ex = diagnostico_extremos(train[paneles["ampliado"]])
        extremos.append(ex.assign(bloque=b))
        # Se prueba RobustScaler solo si en este entrenamiento una SD cae al menos 25 %
        # por un valor fuera de las cercas. Nunca se descarta ese punto del ajuste.
        usar_robusto = ((ex.caida_SD >= 0.25) & (ex.fuera_cercas_IQR > 0)).any()
        iqr = train.quantile(.75) - train.quantile(.25)
        for esc in ["A", "B"]:
            x = train if esc == "A" else pd.concat([train, mmc])
            comunes = paneles["comun"]
            uno = StandardScaler().fit_transform(x[comunes])
            dos = StandardScaler().fit_transform(x[comunes].div(ref[comunes]))
            escalas.append(dict(bloque=b, escenario=esc, max_delta_Z=np.nanmax(np.abs(uno - dos))))
            for metodo in CONFIG:
                previo = 0 if metodo == "Mediana" else int(historicos.loc[str(b), esc, metodo])
                fija_A = amplio["A", metodo, "mM"]
                # Cada control cambia una sola decisión frente al comparador documentado.
                disenos = [
                    ("ratio_fijo", "comun", previo, x, ref, False),
                    ("mM_fijo", "comun", previo, x, None, False),
                    ("seleccion_R", "comun", comun[esc, metodo, "R"], x, None, False),
                    ("seleccion_mM", "comun", comun[esc, metodo, "mM"], x, None, False),
                    ("panel_fijo", "ampliado", comun[esc, metodo, "mM"], x, None, False),
                    ("principal", "ampliado", amplio[esc, metodo, "mM"], x, None, False),
                    ("mismas_opciones_A", "ampliado", fija_A, x, None, False),
                    ("minimo_8", "minimo_8", amplio[esc, metodo, "mM"], x, None, False),
                    ("minimo_12", "minimo_12", amplio[esc, metodo, "mM"], x, None, False),
                ]
                # La sensibilidad de WT añade solo datos, preserva columnas y parámetros.
                wt_train = wt if "WT" in set(grupos.loc[train.index]) else wt.iloc[:0]
                disenos.append(("WT_extra", "ampliado", amplio[esc, metodo, "mM"],
                                pd.concat([x, wt_train]), None, False))
                if usar_robusto:
                    disenos.append(("robusto", "ampliado", amplio[esc, metodo, "mM"], x, None, True))
                if esc == "B":
                    for nombre, aux in [("solo_glucosa", mmc.loc[["MMC2_GLUCOSE"]]),
                                        ("solo_galactosa", mmc.loc[["MMC2_GALACTOSE"]])]:
                        disenos.append((nombre, "ampliado", fija_A, pd.concat([train, aux]), None, False))
                cache = {}
                for nombre, panel, k, entrenamiento, referencia, robusto in disenos:
                    llave = (panel, k, tuple(entrenamiento.index), referencia is not None, robusto)
                    if llave not in cache:
                        assert set(grupos.loc[entrenamiento.index]).isdisjoint(grupos.loc[q.muestra])
                        cache[llave] = predecir(entrenamiento, ish, q, paneles[panel], metodo,
                                                 k, referencia, robusto)
                    p, info = cache[llave]
                    validacion.append(p.assign(diseno=nombre, escenario=esc, metodo=metodo,
                                               grupo=p.muestra.map(grupos), IQR_train=p.metabolito.map(iqr),
                                               referencia_mM=p.metabolito.map(ref)))
                    ajustes.append(dict(diseno=nombre, bloque=b, escenario=esc, metodo=metodo,
                                        k=k, n_train=len(entrenamiento), n_variables=len(paneles[panel]),
                                        **info))
        print(f"Bloque externo {b + 1}/5 terminado", flush=True)
    pred = pd.concat(validacion, ignore_index=True)
    tabla_ajustes = pd.DataFrame(ajustes)
    pred.to_csv(SALIDA / "predicciones.csv.gz", index=False,
                compression={"method": "gzip", "mtime": 0})
    tabla_ajustes.to_csv(SALIDA / "ajustes.csv", index=False)
    pd.concat(soporte).to_csv(SALIDA / "elegibilidad_bloques.csv", index=False)
    pd.concat(extremos).to_csv(SALIDA / "extremos.csv", index=False)
    pd.DataFrame(escalas).to_csv(SALIDA / "equivalencia_Z.csv", index=False)
    from .resumir import resumir
    resumir(pred)
    aporte_individual(ish, mmc, grupos, ref, archivo, pred, tabla_ajustes)
    elegidos = elegir(ish, archivo["internas"]["final"], mmc, grupos, ref, True, "final", busqueda)
    pd.DataFrame(busqueda).to_csv(SALIDA / "busqueda.csv", index=False)
    pd.DataFrame([dict(escenario=esc, metodo=metodo, k=k)
                  for (esc, metodo, criterio), k in elegidos.items() if criterio == "mM"]).to_csv(
        SALIDA / "parametros_finales.csv", index=False)
    entradas = [ARCHIVO_MASCARAS, RAIZ / "sources/ishii.xlsx", RAIZ / "sources/mmc2.xlsx",
                *[RAIZ / "auditoria" / f for f in ("panel_documental_mM.csv", "metadatos.csv",
                  "referencias_glucosa.csv", "homologacion_maestra.csv", "mediciones_auditadas.csv")],
                RAIZ / "evaluacion/parametros_seleccionados.csv", RAIZ / "evaluacion/validacion_celdas.csv",
                SALIDA / "predicciones.csv.gz", SALIDA / "parametros_finales.csv", SALIDA / "matriz_mM.csv"]
    registro = dict(python=platform.python_version(), paquetes={p: version(p) for p in [
        "pandas", "numpy", "scikit-learn", "scipy", "openpyxl", "threadpoolctl"]},
        semilla=SEMILLA, hilos_BLAS=1, configuraciones=CONFIG,
        sha256={str(p.relative_to(RAIZ)): hashlib.sha256(p.read_bytes()).hexdigest() for p in entradas})
    (SALIDA / "ejecucion.json").write_text(json.dumps(registro, indent=2), encoding="utf-8")


def validar():
    with threadpool_limits(limits=1):
        ejecutar()


def completar():
    """Aplica en una sesión nueva o en el notebook, solo con una validación compatible."""
    from .resumir import aplicar

    informe = SALIDA / "ejecucion.json"
    if not informe.exists():
        raise ValueError("Ejecute primero la sección de validación.")
    for nombre, esperado in json.loads(informe.read_text())["sha256"].items():
        archivo = RAIZ / nombre
        if not archivo.exists() or hashlib.sha256(archivo.read_bytes()).hexdigest() != esperado:
            raise ValueError(f"Cambió {nombre}; vuelva a ejecutar la validación.")
    ish, mmc, _, meta, ref = preparar()
    pred = pd.read_csv(SALIDA / "predicciones.csv.gz")
    parametros = pd.read_csv(SALIDA / "parametros_finales.csv")
    elegidos = {(f.escenario, f.metodo, "mM"): f.k for f in parametros.itertuples()}
    with threadpool_limits(limits=1):
        aplicar(ish, mmc, meta.grupo, ref, pred, elegidos)


if __name__ == "__main__":
    validar()
    completar()
