"""Sensibilidades sobre la misma validación de Ishii; no escribe matrices imputadas."""
import hashlib
import json
from importlib.metadata import version

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from .fisiologia import REVISION, auditar
from .modelos import CONFIG, SEMILLA, elegibilidad, predecir
from .preparar import RAIZ, SALIDA, preparar
from .validar import ARCHIVO_MASCARAS, bloques


def sensibilidad(ish, mmc, wt, meta, ref, fis, mascaras, opciones):
    """La unidad de comparación es una familia: mismos cultivos, columnas y opciones."""
    ish, mmc, wt = ish[ref.index], mmc[ref.index], wt[ref.index]
    grupos = meta.grupo
    predicciones, ajustes = [], []
    for b, train, q in bloques(mascaras["externas"], ish, grupos, ref.index, 5):
        # El IQR de evaluación se congela ANTES de las alternativas del bloque:
        # siempre las 24 filas de entrenamiento Ishii; jamás el grupo retirado ni MMC2.
        iqr = train.quantile(.75) - train.quantile(.25)
        for familia in ("escala", "captacion", "flujo", "dilucion_WT"):
            sub = train
            aux = mmc
            consulta = q
            if familia == "captacion":
                sub = train.loc[fis.loc[train.index, "q_C"].notna()]
                consulta = q[q.muestra.map(fis.q_C).notna()]
            elif familia == "flujo":
                disponible = fis[["v_rel", "v_abs"]].notna().all(axis=1)
                sub = train.loc[disponible.loc[train.index]]
                consulta = q[q.muestra.map(disponible)]
            elif familia == "dilucion_WT" and "WT" in set(grupos.loc[train.index]):
                sub = pd.concat([train, wt])
            e = elegibilidad(sub, grupos)
            cols = list(e.index[e.elegible])
            for esc in (("A", "B") if familia in ("escala", "captacion") else ("A",)):
                x = sub if esc == "A" else pd.concat([sub, aux])
                assert set(grupos.loc[x.index]).isdisjoint(grupos.loc[q.muestra])
                estrategias = ["mM_Z"]
                if familia == "escala":
                    estrategias += ["Gerosa_Z", "robusta"]
                elif familia == "captacion":
                    estrategias += ["cov_q_C"]
                    if esc == "A":
                        estrategias += ["indice_q"]  # Nunca se reescalan filas de Gerosa por captación.
                    else:
                        estrategias += ["solo_origen"]  # Control: el estudio es conocido sin medir tasas.
                elif familia == "flujo":
                    estrategias += ["cov_v_rel", "cov_v_abs"]
                else:
                    estrategias += ["cov_D"]
                for metodo in CONFIG:
                    k = int(opciones.loc[b, esc, metodo])
                    for estrategia in estrategias:
                        datos, preguntas = x, ish
                        cov = None
                        mediana_q = np.nan
                        if estrategia.startswith("cov_"):
                            cov = fis[estrategia.removeprefix("cov_")]
                        if estrategia == "solo_origen":
                            cov = fis.cohorte.eq("MMC2").astype(float)
                        if estrategia == "indice_q":
                            # Se ensaya proporcionalidad con q, no se asume una ley cinética.
                            mediana_q = fis.loc[sub.index, "q_sustrato"].median()
                            factor = fis.q_sustrato / mediana_q
                            datos = x.div(factor.loc[x.index], axis=0)
                            preguntas = ish.div(factor.loc[ish.index], axis=0)
                        p, info = predecir(datos, preguntas, consulta, cols, metodo, k,
                                           referencia=ref if estrategia == "Gerosa_Z" else None,
                                           robusto=estrategia == "robusta", covariable=cov)
                        if estrategia == "indice_q":
                            p["pred_mM"] *= p.muestra.map(factor)
                        # Se conservan los 457 intentos, incluidos los de muestras sin tasa.
                        pred = q.merge(p[["celda", "pred_mM", "motivo"]], on="celda", how="left",
                                       validate="one_to_one")
                        pred["motivo"] = pred.motivo.fillna("sin_fisiologia_requerida")
                        predicciones.append(pred.assign(
                            familia=familia, estrategia=estrategia, escenario=esc, metodo=metodo,
                            grupo=pred.muestra.map(grupos), IQR_train=pred.metabolito.map(iqr),
                            referencia_mM=pred.metabolito.map(ref)))
                        ajustes.append(dict(
                            bloque=b, familia=familia, estrategia=estrategia, escenario=esc,
                            metodo=metodo, k=k, n_train_Ishii=len(sub), n_train_MMC2=len(x)-len(sub),
                            cultivos_prueba=q.muestra.nunique(), cultivos_con_dato=consulta.muestra.nunique(),
                            ids_entrenamiento="|".join(x.index), metabolitos="|".join(cols),
                            covariable_activa=(cov is not None and metodo != "Mediana"
                                               and cov.loc[x.index].nunique() > 1),
                            min_cov_train=cov.loc[x.index].min() if cov is not None else np.nan,
                            max_cov_train=cov.loc[x.index].max() if cov is not None else np.nan,
                            referencia_q_Ishii_train=mediana_q, **info))
    return pd.concat(predicciones, ignore_index=True), pd.DataFrame(ajustes)


def revisar():
    from .revision_metricas import figura_pesos, resumir_revision

    ish, mmc, wt, meta, ref = preparar()
    fis = auditar(ish.index, wt.index, mmc.index, meta)
    manifiesto = json.loads((SALIDA / "ejecucion.json").read_text(encoding="utf-8"))
    for nombre, huella in manifiesto["sha256"].items():
        if hashlib.sha256((RAIZ / nombre).read_bytes()).hexdigest() != huella:
            raise ValueError(f"Ha cambiado {nombre}; ejecute primero python -m ampliacion.validar.")
    mascaras = json.loads(ARCHIVO_MASCARAS.read_text())
    prev = pd.read_csv(SALIDA / "predicciones.csv.gz", float_precision="round_trip")
    opciones = pd.read_csv(SALIDA / "ajustes.csv").query("diseno == 'seleccion_mM'").set_index(
        ["bloque", "escenario", "metodo"]).k
    # Las opciones deben provenir de los tres bloques internos en el panel común, no de la prueba.
    busqueda = pd.read_csv(SALIDA / "busqueda.csv").query("panel == 'comun' and bloque != 'final'").copy()
    busqueda["bloque"] = busqueda.bloque.astype(int)
    for llave, datos in busqueda.groupby(["bloque", "escenario", "metodo"]):
        assert opciones.loc[llave] == datos.loc[datos.macro_mM.idxmin(), "k"]
    assert len(opciones) == 30 and opciones.index.is_unique
    with threadpool_limits(limits=1):
        pred, ajustes = sensibilidad(ish, mmc, wt, meta, ref, fis, mascaras, opciones)
    llaves = ["celda", "metodo", "escenario"]
    base = pred.query("familia == 'escala' and estrategia == 'mM_Z'").set_index(llaves).sort_index()
    previo = prev[prev.diseno.eq("seleccion_mM")].set_index(llaves).sort_index()
    assert base.index.equals(previo.index)
    np.testing.assert_allclose(base.pred_mM, previo.pred_mM, atol=1e-10, rtol=1e-10, equal_nan=True)
    pred.to_csv(REVISION / "predicciones.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    ajustes.to_csv(REVISION / "ajustes.csv", index=False)
    resumir_revision(pred, prev, ish, ref)
    figura_pesos()
    archivos = [ARCHIVO_MASCARAS, RAIZ / "sources/ishii.xlsx", RAIZ / "sources/mmc2.xlsx",
                RAIZ / "sources/ecoli_auditada_v3.db", SALIDA / "ejecucion.json", SALIDA / "ajustes.csv",
                SALIDA / "busqueda.csv", SALIDA / "predicciones.csv.gz"]
    ficha = {"semilla": SEMILLA, "hilos": 1, "paquetes": {p: version(p) for p in
             ("numpy", "pandas", "scikit-learn", "scipy", "openpyxl", "threadpoolctl")},
             "objetivos": list(ref.index), "hiperparametros": "seleccion_mM, panel común, 3 bloques internos",
             "sha256": {str(p.relative_to(RAIZ)): hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in archivos}}
    (REVISION / "ejecucion.json").write_text(json.dumps(ficha, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    revisar()
