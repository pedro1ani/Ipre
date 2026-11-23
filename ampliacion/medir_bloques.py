"""Error y cobertura emparejados sin contar una misma concentración como varias muestras."""
import numpy as np
import pandas as pd


def resumen_celdas(t, por):
    claves = ["experimento", "alternativa", "escenario", "metodo", "alcance"] + por
    g = t.groupby(claves, sort=False)
    r = g.agg(intentos=("intentos", "sum"), eventos_predichos=("eventos", "sum"),
              celdas_intentadas=("celda", "nunique"), celdas_predichas=("abs", "count"),
              grupos_con_intentos=("grupo", "nunique"), grupos_evaluados=("grupo_valido", "nunique"),
              n_IQR=("rel", "count"),
              MAE_mM=("abs", "mean"), RMSE_mM=("cuad", "mean"), sesgo_mM=("sesgo", "mean"),
              MAE_IQR=("rel", "mean"), MAE_R=("gerosa", "mean")).reset_index()
    r["RMSE_mM"] = np.sqrt(r.RMSE_mM)
    r["cobertura_eventos"] = r.eventos_predichos / r.intentos
    return r


def medir(predicciones):
    partes = []
    casos = ["bloque", "rep", "evento"]
    modelos = ["alternativa", "escenario", "metodo"]
    for nombre, tabla in predicciones.groupby("experimento", sort=False):
        # Un intento permanece en el denominador aunque el modelo se abstenga.
        wide = tabla.pivot(index=casos, columns=modelos, values="pred_mM")
        assert len(tabla) == wide.size
        comunes = wide.notna().all(axis=1)
        t = tabla.copy()
        t["comun"] = comunes.reindex(pd.MultiIndex.from_frame(t[casos])).to_numpy()
        for alcance in ("disponibles", "comunes"):
            v = t.copy()
            if alcance == "comunes":
                v.loc[~v.comun, "pred_mM"] = np.nan
            err = v.pred_mM - v.verdad_mM
            v = v.assign(abs=err.abs(), cuad=err.pow(2), sesgo=err,
                         rel=err.abs() / v.IQR_train.where(v.IQR_train > 0),
                         gerosa=err.abs() / v.referencia_mM)
            # Se promedian errores absolutos/cuadráticos por celda, no predicciones.
            c = v.groupby(modelos + ["metabolito", "grupo", "celda"], sort=False).agg(
                intentos=("evento", "size"), eventos=("pred_mM", "count"),
                **{k: (k, "mean") for k in ["abs", "cuad", "sesgo", "rel", "gerosa"]}).reset_index()
            partes.append(c.assign(experimento=nombre, alcance=alcance,
                                   grupo_valido=c.grupo.where(c.eventos.gt(0))))
    celdas = pd.concat(partes, ignore_index=True)
    met = resumen_celdas(celdas, ["metabolito"])
    grupos = resumen_celdas(celdas, ["grupo", "metabolito"])
    glob = resumen_celdas(celdas, [])
    clave = ["experimento", "alternativa", "escenario", "metodo", "alcance"]
    indicadores = ["MAE_mM", "RMSE_mM", "sesgo_mM", "MAE_IQR", "MAE_R"]
    macro = met.groupby(clave)[indicadores].mean().add_suffix("_macro")
    macro["metabolitos_evaluados"] = met.groupby(clave).MAE_mM.count()
    glob = glob.rename(columns={k: k+"_micro" for k in indicadores}).merge(
        macro.reset_index(), on=clave, validate="one_to_one")
    return glob, met, grupos
