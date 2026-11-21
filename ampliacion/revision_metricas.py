"""Errores sobre celdas pareadas y dispersión del MAE entre grupos de perturbación."""
import numpy as np
import pandas as pd

from .fisiologia import REVISION
from .resumir import metricas

CLAVES = ["alcance", "familia", "escenario", "metodo", "estrategia"]
INDICADORES = ["MAE_mM", "RMSE_mM", "MAE_IQR", "sesgo_mM", "MAE_R"]


def emparejar(t):
    """Todos los métodos ven las mismas verdades; se distinguen respuesta y comparación."""
    llave = ["familia", "escenario", "metodo", "estrategia"]
    assert not t.duplicated(llave + ["celda"]).any()
    assert t.groupby("celda")[
        ["verdad_mM", "grupo", "bloque", "IQR_train", "referencia_mM"]
    ].nunique(dropna=False).eq(1).all().all()
    total = t.celda.nunique()
    assert t.groupby(llave).size().eq(total).all()
    visible = t.pred_mM.notna()
    familia = visible.groupby([t.familia, t.celda]).transform("all")
    todas = visible.groupby(t.celda).transform("all")
    return pd.concat([t.assign(alcance=nombre, respondio=visible, pred_mM=t.pred_mM.where(mask))
                      for nombre, mask in [("comunes_familia", familia), ("comunes_todas", todas)]],
                     ignore_index=True)


def estabilidad(t):
    v = t.assign(error=(t.pred_mM-t.verdad_mM).abs()).groupby(CLAVES + ["metabolito", "grupo"])
    por_grupo = v.error.agg(n_pares="count", MAE_mM="mean").reset_index()
    # Global de cada grupo: primero MAE por metabolito y luego el promedio de esos MAE.
    macro = por_grupo.groupby(CLAVES + ["grupo"]).agg(n_pares=("n_pares", "sum"),
                                                      MAE_mM=("MAE_mM", "mean")).reset_index()
    filas = pd.concat([por_grupo, macro.assign(metabolito="GLOBAL_MACRO")], ignore_index=True)
    base = filas[filas.estrategia.eq("mM_Z")].drop(columns=["n_pares", "estrategia"])
    llaves = [c for c in CLAVES if c != "estrategia"] + ["grupo", "metabolito"]
    filas = filas.merge(base.rename(columns={"MAE_mM": "MAE_control"}), on=llaves, validate="many_to_one")
    filas["delta_mM"] = filas.MAE_mM - filas.MAE_control
    filas.to_csv(REVISION / "por_grupo.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    g = filas.groupby(CLAVES + ["metabolito"])
    s = g.MAE_mM.agg(grupos_con_pares="count", mediana_MAE_grupos="median")
    s["Q25_MAE_grupos"] = g.MAE_mM.quantile(.25).where(s.grupos_con_pares >= 3)
    s["Q75_MAE_grupos"] = g.MAE_mM.quantile(.75).where(s.grupos_con_pares >= 3)
    cuentas = filas.assign(mejora=filas.delta_mM.lt(-1e-10), empeora=filas.delta_mM.gt(1e-10))
    s = s.join(cuentas.groupby(CLAVES + ["metabolito"])[["mejora", "empeora"]].sum())
    s["empatan"] = s.grupos_con_pares - s.mejora - s.empeora
    return s.reset_index()


def resumir_revision(pred, anterior, ish, ref):
    primario = anterior[anterior.diseno.eq("principal")].assign(familia="panel_ampliado", estrategia="mM_Z")
    t = emparejar(pd.concat([pred, primario[pred.columns]], ignore_index=True))
    met = metricas(t, CLAVES + ["metabolito"])
    micro = metricas(t, CLAVES)
    disp = t.groupby(CLAVES).respondio.sum().rename("respuestas_disponibles")
    macro = met.groupby(CLAVES)[INDICADORES].mean().add_suffix("_macro")
    glob = micro.rename(columns={c: c + "_micro" for c in INDICADORES}).merge(macro, on=CLAVES)
    glob = glob.merge(disp, on=CLAVES).merge(
        met.groupby(CLAVES).MAE_mM.count().rename("metabolitos_evaluados"), on=CLAVES)
    met = met.merge(t.groupby(CLAVES + ["metabolito"]).respondio.sum().rename("respuestas_disponibles"),
                    on=CLAVES + ["metabolito"])
    for tabla in (met, glob):
        tabla.rename(columns={"respondidas": "n_pares", "cobertura": "fraccion_comun"}, inplace=True)
        tabla["cobertura_modelo"] = tabla.respuestas_disponibles / tabla.intentos
    est = estabilidad(t)
    glob = glob.merge(est[est.metabolito.eq("GLOBAL_MACRO")].drop(columns="metabolito"), on=CLAVES)
    met = met.merge(est[est.metabolito.ne("GLOBAL_MACRO")], on=CLAVES + ["metabolito"])
    glob.to_csv(REVISION / "global.csv", index=False)
    met.to_csv(REVISION / "por_metabolito.csv", index=False)
    # Estas tres tablas permiten revisar los resultados sin distribuir registros de cultivos.
    columnas = CLAVES + ["n_pares", "grupos_con_pares", "metabolitos_evaluados", "MAE_mM_macro",
                         "MAE_mM_micro", "RMSE_mM_micro", "MAE_IQR_macro", "sesgo_mM_micro",
                         "MAE_R_macro", "Q25_MAE_grupos", "Q75_MAE_grupos", "mejora", "empeora"]
    glob[columnas].to_csv(REVISION / "resumen_comparaciones.csv", index=False, float_format="%.6g")
    principal = met.query("alcance == 'comunes_familia' and familia == 'panel_ampliado'")
    columnas = ["escenario", "metodo", "metabolito", "intentos", "n_pares", "grupos_con_pares",
                *INDICADORES, "mediana_MAE_grupos", "Q25_MAE_grupos", "Q75_MAE_grupos"]
    principal[columnas].to_csv(REVISION / "error_panel_ampliado.csv", index=False, float_format="%.6g")
    captacion = met.query("alcance == 'comunes_familia' and familia == 'captacion'")
    ids = ["escenario", "metodo", "metabolito"]
    control = captacion[captacion.estrategia.eq("mM_Z")].set_index(ids).MAE_mM
    con_tasa = captacion[captacion.estrategia.eq("cov_q_C")].set_index(ids)
    aporte = con_tasa[["n_pares", "grupos_con_pares", "MAE_mM", "mejora", "empeora"]].rename(
        columns={"MAE_mM": "MAE_con_captacion_mM"})
    aporte["MAE_control_mM"] = control
    aporte["delta_mM"] = aporte.MAE_con_captacion_mM - aporte.MAE_control_mM
    aporte.to_csv(REVISION / "aporte_captacion.csv", float_format="%.6g")
    # Contribuciones al macro: la referencia es fija en todos los métodos; nunca divide por la verdad.
    pesos = met.query("alcance == 'comunes_familia' and estrategia == 'mM_Z' "
                      "and familia in ['panel_ampliado','escala']").copy()
    pesos["referencia_mM"] = pesos.metabolito.map(ref)
    pesos["mediana_Ishii_mM"] = pesos.metabolito.map(ish.median())
    pesos["mediana_Ishii_sobre_Gerosa"] = pesos.mediana_Ishii_mM / pesos.referencia_mM
    pesos["peso_por_mM"] = 1 / pesos.referencia_mM
    for indicador in ("MAE_mM", "MAE_R"):
        total = pesos.groupby(CLAVES)[indicador].transform("sum")
        pesos["aporte_pct_" + indicador] = 100 * pesos[indicador] / total.where(total > 0)
    pesos[["familia", "escenario", "metodo", "metabolito", "n_pares", "referencia_mM",
           "mediana_Ishii_mM", "mediana_Ishii_sobre_Gerosa", "peso_por_mM", "MAE_mM",
           "MAE_R", "aporte_pct_MAE_mM", "aporte_pct_MAE_R"]].to_csv(REVISION / "pesos_Gerosa.csv",
                                                                    index=False)
    copia = pred[pred.familia.eq("escala")].pivot(
        index=["escenario", "metodo", "bloque", "celda"], columns="estrategia", values="pred_mM")
    assert copia.mM_Z.isna().equals(copia.Gerosa_Z.isna())
    np.testing.assert_allclose(copia.mM_Z, copia.Gerosa_Z, atol=1e-10, rtol=1e-10, equal_nan=True)
    (copia.mM_Z - copia.Gerosa_Z).abs().groupby(level=["escenario", "metodo"]).agg(
        n_pares="count", max_delta_escala_mM="max").to_csv(REVISION / "equivalencia.csv")


def figura_pesos():
    import matplotlib.pyplot as plt

    g = pd.read_csv(REVISION / "global.csv").query(
        "alcance == 'comunes_familia' and familia == 'panel_ampliado' and escenario == 'B' "
        "and estrategia == 'mM_Z'").set_index("metodo").loc[["Mediana", "KNN", "Iterative"]]
    fig, ejes = plt.subplots(1, 2, figsize=(8.2, 3.1), layout="constrained")
    for eje, col, etiqueta in zip(ejes, ["MAE_mM_macro", "MAE_R_macro"],
                                  ["Error físico (mM)", "Error / referencia de Gerosa"]):
        barras = eje.barh(g.index, g[col], color=["#777777", "#176881", "#ad5025"])
        eje.invert_yaxis()
        eje.set_xlim(0, 1.24*g[col].max())
        eje.set_title(etiqueta)
        eje.bar_label(barras, labels=[f"{v:.3f}".replace(".", ",") for v in g[col]], padding=3)
        eje.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Mismas {int(g.n_pares.iloc[0])} celdas · panel ampliado con MMC2 · menor es mejor")
    fig.savefig(REVISION / "efecto_metrica.png", dpi=150)
    plt.close(fig)
