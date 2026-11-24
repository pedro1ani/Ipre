"""Decisiones para una casilla. Solo los errores de los cultivos de entrenamiento enseñan la regla."""
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from .datos import PRIV, leer_matriz
from .diseno import bloques, consultas, elegibilidad, escoger_companeros
from .modelos import METODOS, estimar, matriz_consultas
from .validar import leer_split, particiones

DET = PRIV / "revision"
RES = DET / "resultados"
FIJOS = list(METODOS)
PRINCIPALES = [*FIJOS, "Original", "Respaldo_KNN", "Respaldo_interno", "Jerarquica"]
# Se ensayan primero las reglas más simples. La elección definitiva usa solo pruebas internas.
CONFIG = {
    "Global": ("", 10, 1), "Similares": ("s", 10, 1), "Metabolito": ("m", 10, 1),
    "Similar_met": ("sm", 10, 1), "Completa": ("smp", 10, 1),
    "Completa_n8": ("smp", 8, 1), "Completa_n12": ("smp", 12, 1),
    "Completa_sin_SE": ("smp", 10, 0),
}
CLAVES = ["cultivo", "objetivo", "perfil", "grupo", "fold", "tipo"]


def guardar(tabla, nombre):
    tabla.to_csv(DET / nombre, index=False, compression={"method":"gzip","mtime":0} if nombre.endswith(".gz") else None)


def comunes(p):
    return p.loc[p.opcion.eq("A_amplio")].pivot(index="evento", columns="metodo",
                                                values="pred").reindex(columns=FIJOS)


def errores_internos(q, p):
    z = q[q.tipo.ne("completa")].drop_duplicates(["cultivo", "objetivo", "perfil"])
    z = z.join(p, on="evento")
    z[FIJOS] = z[FIJOS].sub(z.verdad, axis=0).abs()
    return z[CLAVES + FIJOS].dropna(subset=FIJOS)


def subinterno(ish, meta, ids, test, etiqueta):
    """Genera evidencia para ensayar configuraciones sin conocer los cultivos del bloque interno."""
    tr = ish.loc[ids]
    e = elegibilidad(tr, meta.grupo)
    eleg = e.index[e.elegible]
    q = consultas(tr, ish.loc[test])
    q = q[q.tipo.ne("completa")].drop_duplicates(["cultivo", "objetivo", "perfil"]).copy()
    q["grupo"], q["fold"] = q.cultivo.map(meta.grupo), etiqueta
    datos = []
    for y, g in q.groupby("objetivo"):
        if y not in eleg:
            continue
        cols = [y] + escoger_companeros(tr, eleg, y)
        if len(cols) < 3:
            continue
        fit = estimar(tr[cols].to_numpy(), matriz_consultas(ish, g, cols), len(cols))
        verdad = ish.loc[g.cultivo, y].to_numpy()
        z = g[CLAVES].copy()
        for m in FIJOS:
            z[m] = np.abs(fit[m][0] - verdad)
        datos.append(z.dropna(subset=FIJOS))
    return pd.concat(datos, ignore_index=True) if datos else pd.DataFrame(columns=CLAVES + FIJOS)


def clases(tr):
    # Separar por cobertura y variabilidad relativas evita agrupar por magnitud en mM.
    falta = tr.isna().mean()
    cv = tr.std().div(tr.mean().where(tr.mean().gt(1e-12)))
    return {y: (int(falta[y] > falta.median()), int(cv[y] > cv.median()))
            for y in tr if pd.notna(cv[y])}


def error_estandar(diferencia, grupo):
    """Error estándar de una media por casilla con conglomerados de tamaño desigual."""
    d = np.asarray(diferencia, dtype=float)
    codigo, g = pd.factorize(grupo)
    if len(g) < 2 or not len(d):
        return np.inf
    suma = np.bincount(codigo, weights=d - d.mean())
    return float(np.sqrt(len(g) / (len(g) - 1) * (suma @ suma)) / len(d))


def evidencia(z):
    if z.empty:
        return None
    errores = z[FIJOS].to_numpy()
    g = z.grupo.to_numpy()
    media = errores.mean(axis=0)
    mejor = int(np.argmin(media))
    por_bloque = z.groupby("fold")[FIJOS].mean().to_numpy()
    se = [error_estandar(errores[:, i] - errores[:, mejor], g) for i in range(len(FIJOS))]
    return (mejor, len(z), z.grupo.nunique(), len(por_bloque), media - media[mejor],
            (por_bloque[:, [mejor]] < por_bloque - 1e-12).sum(axis=0), se)


def preparar_evidencia(e, tr):
    individual = e[e.tipo.eq("individual")]
    cat = clases(tr)
    agrupados = {y: g for y, g in individual.groupby("objetivo")}
    simil = {}
    for y, c in cat.items():
        otros = [k for k in agrupados if k != y and cat.get(k) == c]
        if len(otros) >= 3:
            simil[y] = evidencia(pd.concat([agrupados[k] for k in otros], ignore_index=True))
    return {"g": evidencia(individual), "s": simil,
            "m": {y: evidencia(g) for y, g in agrupados.items()},
            "p": {k: evidencia(g) for k, g in e.groupby(["objetivo", "perfil"])}}


def elegir(ev, y, perfil, config):
    niveles, minimo, margen = CONFIG[config]
    met, nivel = 0, "mediana_sin_evidencia"
    global_ = ev["g"]
    if global_ and global_[1] >= minimo and global_[2] >= 5 and global_[3] >= 2:
        met, nivel = global_[0], "global"
    for clave, nombre in (("s", "similares"), ("m", "metabolito"), ("p", "metabolito_perfil")):
        if clave not in niveles:
            continue
        a = ev[clave].get((y, perfil) if clave == "p" else y)
        if not a or a[1] < minimo or a[2] < 5 or a[3] < 2:
            continue
        candidato, _, _, _, ganancia, victorias, se = a
        if candidato != met and victorias[met] >= 2 and ganancia[met] > margen * se[met] + 1e-12:
            met, nivel = candidato, nombre
    return FIJOS[met], nivel


def elegir_config(errores):
    # Cada configuración se prueba sobre idénticas celdas internas; regla de un error estándar pareado.
    a = errores.dropna(subset=list(CONFIG))
    if a.empty or a.grupo.nunique() < 5:
        return "Global"
    mejor = a[list(CONFIG)].mean().idxmin()
    for c in CONFIG:
        d = a[c] - a[mejor]
        if d.mean() <= error_estandar(d, a.grupo) + 1e-12:
            return c
    return mejor


def aprender_fold(args):
    nombre, ids = args
    ish, meta = leer_matriz("ishii"), leer_matriz("metadatos")
    intern, evaluacion = [], []
    for i, (tr, te) in enumerate(bloques(ids, meta.grupo, 3)):
        q, p = leer_split(f"{nombre}_i{i}", ish, meta)
        assert set(q.cultivo) <= set(te)
        errores = errores_internos(q, comunes(p))
        intern.append(errores)
        sub = [subinterno(ish, meta, a, b, j) for j, (a, b) in enumerate(bloques(tr, meta.grupo, 3))]
        con_datos = [z for z in sub if not z.empty]
        ev = preparar_evidencia(pd.concat(con_datos) if con_datos else sub[0], ish.loc[tr])
        z = errores[errores.tipo.eq("individual")].copy()
        for c in CONFIG:
            met = [elegir(ev, y, pf, c)[0] for y, pf in zip(z.objetivo, z.perfil)]
            z[c] = z[FIJOS].to_numpy()[np.arange(len(z)), [FIJOS.index(m) for m in met]]
        evaluacion.append(z)
    riesgo = pd.concat(evaluacion)
    config = elegir_config(riesgo)
    resumen = pd.DataFrame({"configuracion": list(CONFIG), "MAE_interno_mM": riesgo[list(CONFIG)].mean(),
                            "casillas": len(riesgo), "grupos": riesgo.grupo.nunique(),
                            "elegida": config}, index=list(CONFIG))
    errores = pd.concat(intern)
    ev = preparar_evidencia(errores, ish.loc[ids])
    # Tamaño de muestra antes del requisito de ganarle a todos.
    conteos = errores.groupby(["objetivo", "perfil"]).agg(n=("cultivo","nunique"),
            grupos=("grupo","nunique"), bloques=("fold","nunique")).reset_index()
    return nombre, ev, config, resumen.assign(particion=nombre), conteos.assign(particion=nombre)


def soporte(ish, q, fold):
    ajustes = pd.read_csv(PRIV / "particiones" / fold / "ajustes.csv")
    salida = []
    for aj in ajustes[ajustes.opcion.eq("A_amplio")].itertuples():
        z = q[q.objetivo.eq(aj.objetivo)]
        cols = aj.columnas.split("|")
        tr = ish.loc[aj.ids.split("|"), cols].to_numpy()
        vis = np.isfinite(matriz_consultas(ish, z, cols)[:, 1:])
        pares = vis.astype(int) @ np.isfinite(tr[:, 1:]).T.astype(int)
        n = ((pares >= 2) & np.isfinite(tr[:, 0])).sum(axis=1)
        salida.append(pd.DataFrame({"evento": z.evento, "predictores": vis.sum(axis=1),
                                    "vecinos_disponibles": n}))
    return pd.concat(salida).set_index("evento")


def aplicar_estrategias(q, p, original, ev, config):
    """No consulta verdades ni estadísticas de prueba; todos los valores proceden de ajustes previos."""
    fijos = comunes(p).reindex(q.evento)
    originales = original.set_index("evento").reindex(q.evento)
    reglas = {(y, pf, c): elegir(ev, y, pf, c) for y, pf in
              q[["objetivo", "perfil"]].drop_duplicates().itertuples(index=False, name=None) for c in CONFIG}
    filas = []
    for f, valores, ori in zip(q.itertuples(), fijos.to_dict("records"), originales.itertuples()):
        def fila(nombre, m, valor, nivel, estado):
            if not np.isfinite(valor):
                m, estado = "", "pendiente"
            filas.append((f.evento, nombre, m, valor, nivel, estado))
        for m in FIJOS:
            fila(m, m, valores[m], "fijo", "fijo")
        clave = ori.clave if isinstance(ori.clave, str) else ""
        elegido = clave.split("/")[-1]
        fila("Original", elegido, ori.pred, ori.estado, ori.estado)
        for nom, c in [("Respaldo_KNN", None), ("Respaldo_interno", "Global"),
                       ("Jerarquica", config), *((k, k) for k in CONFIG)]:
            hibrido = nom.startswith("Respaldo")
            if hibrido and ori.estado == "validada" and np.isfinite(ori.pred):
                fila(nom, elegido, ori.pred, "metabolito_perfil", "validada")
                continue
            m, nivel = ("KNN3", "knn") if c is None else reglas[f.objetivo, f.perfil, c]
            if not np.isfinite(valores[m]):
                m, nivel = "Mediana", "mediana_sin_soporte_multivariado"
            estado = "validada" if not hibrido and nivel in ("metabolito", "metabolito_perfil") else "respaldo"
            fila(nom, m, valores[m], nivel, estado)
    return pd.DataFrame(filas, columns=["evento", "estrategia", "metodo", "pred", "nivel", "estado"])


def ejecutar(procesos=4):
    RES.mkdir(exist_ok=True, parents=True)
    ish, meta = leer_matriz("ishii"), leer_matriz("metadatos")
    _, externos = particiones(ish, meta)
    entrenamiento = {n: tr for n, (tr, _) in externos.items()} | {"F": ish.index}
    tablas, riesgos, soportes = [], [], []
    with ProcessPoolExecutor(max_workers=procesos) as pool:
        for n, ev, conf, riesgo, cantidad in pool.map(aprender_fold, entrenamiento.items()):
            riesgos.append(riesgo)
            soportes.append(cantidad)
            fold = "real" if n == "F" else n
            q, p = leer_split(fold, ish, meta, con_verdad=False)
            archivo = "auditoria_huecos.csv" if n == "F" else "externas.csv.gz"
            orig = pd.read_csv(PRIV / archivo, float_precision="round_trip")
            if n != "F":
                orig = orig[orig.fold.eq(n)]
            estimados = aplicar_estrategias(q, p, orig, ev, conf)
            q = q.join(soporte(ish, q, fold), on="evento")
            tr = ish.loc[entrenamiento[n]]
            q["IQR"] = q.objetivo.map(tr.quantile(.75) - tr.quantile(.25))
            q["SD_train"] = q.objetivo.map(tr.std())
            q["faltante_train_pct"] = q.objetivo.map(tr.isna().mean() * 100)
            z = q.merge(estimados, on="evento", validate="one_to_many")
            if n == "F":
                z = z.merge(orig[["evento","ausencia_original","source_cell","source_sheet"]], on="evento")
                no = z.ausencia_original.ne("not_detected_unspecified")
                z.loc[no, ["pred","metodo","estado","nivel"]] = [np.nan,"","pendiente","ausencia_no_autorizada"]
                guardar(z, "huecos.csv.gz")
            else:
                # La verdad externa se incorpora únicamente después de fijar todas las estimaciones.
                z = z.join(orig.set_index("evento").verdad, on="evento")
                assert z.verdad.notna().all()
                tablas.append(z)
            print(n, conf, flush=True)
    guardar(pd.concat(tablas), "validacion.csv.gz")
    pd.concat(riesgos).to_csv(RES / "seleccion_interna.csv", index=False, float_format="%.9g")
    guardar(pd.concat(soportes), "soporte_interno.csv.gz")


if __name__ == "__main__":
    ejecutar()

