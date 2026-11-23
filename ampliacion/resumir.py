"""Métricas sobre mediciones ocultas y, después de evaluarlas, aplicación a huecos reales."""
import numpy as np
import pandas as pd

from .modelos import CONFIG, elegibilidad, predecir
from .preparar import RAIZ, SALIDA


def metricas(tabla, claves):
    t = tabla.copy()
    t["sesgo_mM"] = t.pred_mM - t.verdad_mM
    t["MAE_mM"] = t.sesgo_mM.abs()
    t["MAE_IQR"] = t.MAE_mM / t.IQR_train.where(t.IQR_train > 0)
    t["MAE_R"] = t.MAE_mM / t.referencia_mM
    t["cuadrado"] = t.sesgo_mM.pow(2)
    t["grupo_evaluado"] = t.grupo.where(t.pred_mM.notna())
    g = t.groupby(claves, sort=False)
    r = g.agg(intentos=("celda", "size"), respondidas=("pred_mM", "count"),
              n_grupos=("grupo", "nunique"), grupos_evaluados=("grupo_evaluado", "nunique"),
              n_IQR=("MAE_IQR", "count"), MAE_mM=("MAE_mM", "mean"),
              RMSE_mM=("cuadrado", "mean"), sesgo_mM=("sesgo_mM", "mean"),
              MAE_IQR=("MAE_IQR", "mean"), MAE_R=("MAE_R", "mean")).reset_index()
    r["RMSE_mM"] = np.sqrt(r.RMSE_mM)
    r["cobertura"] = r.respondidas / r.intentos
    return r


def aporte_pareado(t):
    """Cada diferencia usa exactamente las mismas celdas A/B dentro de un grupo y analito."""
    indices = ["diseno", "metodo", "bloque", "grupo", "metabolito", "celda", "verdad_mM"]
    pares = t.pivot(index=indices, columns="escenario", values="pred_mM").reset_index()
    pares["pareada"] = pares.A.notna() & pares.B.notna()
    for esc in ["A", "B"]:
        pares[f"sesgo_{esc}"] = (pares[esc] - pares.verdad_mM).where(pares.pareada)
        pares[f"MAE_{esc}"] = pares[f"sesgo_{esc}"].abs()
    pares["delta_MAE_B_menos_A"] = pares.MAE_B - pares.MAE_A
    pares["delta_sesgo_B_menos_A"] = pares.sesgo_B - pares.sesgo_A
    indices = indices[:5]
    resultado = pares.groupby(indices).agg(
        intentos=("celda", "size"), n_pares=("pareada", "sum"), MAE_A=("MAE_A", "mean"),
        MAE_B=("MAE_B", "mean"), delta_MAE_B_menos_A=("delta_MAE_B_menos_A", "mean"),
        delta_sesgo_B_menos_A=("delta_sesgo_B_menos_A", "mean")).reset_index()
    return resultado



def contrastes(pred):
    """Dos alternativas por vez: siempre se empareja por celda, incluso si cambia la cobertura."""
    resultados = []
    for metodo, datos in pred.groupby("metodo"):
        for esc in ["A", "B"]:
            pares = [
                ("escala", "ratio_fijo", esc, "mM_fijo", esc),
                ("criterio", "seleccion_R", esc, "seleccion_mM", esc),
                ("panel", "seleccion_mM", esc, "panel_fijo", esc),
                ("reajuste", "panel_fijo", esc, "principal", esc),
                ("4_WT", "principal", esc, "WT_extra", esc),
                ("minimo_8", "principal", esc, "minimo_8", esc),
                ("minimo_12", "principal", esc, "minimo_12", esc),
                ("escalador_robusto", "principal", esc, "robusto", esc),
            ]
            if esc == "B":
                pares += [("8_MMC2_sin_reajuste", "mismas_opciones_A", "A", "mismas_opciones_A", "B"),
                          ("8_MMC2_con_reajuste", "principal", "A", "principal", "B"),
                          ("glucosa_sola", "mismas_opciones_A", "A", "solo_glucosa", "B"),
                          ("galactosa_sola", "mismas_opciones_A", "A", "solo_galactosa", "B")]
            for nombre, dis1, esc1, dis2, esc2 in pares:
                base = datos[(datos.diseno == dis1) & (datos.escenario == esc1)].set_index("celda")
                otra = datos[(datos.diseno == dis2) & (datos.escenario == esc2)].set_index("celda")
                if otra.empty:
                    continue
                base = base.loc[otra.index].copy()
                np.testing.assert_array_equal(base.verdad_mM, otra.verdad_mM)
                base["pareada"] = base.pred_mM.notna() & otra.pred_mM.notna()
                base["p2"] = otra.pred_mM
                for marca, col in [("1", "pred_mM"), ("2", "p2")]:
                    error = (base[col] - base.verdad_mM).where(base.pareada)
                    base["sesgo_" + marca] = error
                    base["MAE_" + marca] = error.abs()
                base["delta_MAE"] = base.MAE_2 - base.MAE_1
                base["delta_sesgo"] = base.sesgo_2 - base.sesgo_1
                for estrato, sub in [("todos", base), ("WT", base[base.grupo.eq("WT")]),
                                     ("KO", base[base.grupo.ne("WT")])]:
                    r = sub.groupby("metabolito").agg(
                        intentos=("pareada", "size"), n_pares=("pareada", "sum"),
                        n_1=("pred_mM", "count"), n_2=("p2", "count"),
                        MAE_1=("MAE_1", "mean"), MAE_2=("MAE_2", "mean"),
                        delta_MAE=("delta_MAE", "mean"), delta_sesgo=("delta_sesgo", "mean"))
                    r = r.reset_index().assign(contraste=nombre, metodo=metodo, escenario=esc,
                                               estrato=estrato, primero=dis1+"_"+esc1, segundo=dis2+"_"+esc2)
                    resultados.append(r)
    pd.concat(resultados, ignore_index=True).to_csv(SALIDA / "contrastes.csv", index=False)


def resumir(pred):
    claves = ["diseno", "escenario", "metodo"]
    assert not pred.duplicated(claves + ["celda"]).any()
    # En un bloque sin indicios de influencia, el control robusto no se ejecuta.
    completo = pred[~pred.diseno.eq("robusto")]
    celdas = set(completo.celda)
    assert completo.groupby(claves).celda.apply(lambda v: set(v) == celdas).all()
    principal = pred[pred.diseno.eq("principal")].copy()
    comun = principal.pivot(index="celda", columns=["escenario", "metodo"], values="pred_mM")
    assert comun.shape[1] == 6
    todas = principal.copy()
    todas.loc[~todas.celda.isin(comun.index[comun.notna().all(axis=1)]), "pred_mM"] = np.nan
    t = pd.concat([pred.assign(alcance="disponibles"), todas.assign(alcance="comunes_seis")])
    claves = ["alcance"] + claves
    met = metricas(t, claves + ["metabolito"])
    met.to_csv(SALIDA / "por_metabolito.csv", index=False)
    micro = metricas(t, claves)
    valores = ["MAE_mM", "RMSE_mM", "sesgo_mM", "MAE_IQR", "MAE_R"]
    macro = met.groupby(claves)[valores].mean().add_suffix("_macro").reset_index()
    n_met = met.groupby(claves).MAE_mM.count().rename("metabolitos_evaluados")
    micro.rename(columns={v: v + "_micro" for v in valores}).merge(macro, on=claves).merge(
        n_met, on=claves).to_csv(SALIDA / "global.csv", index=False)
    metricas(principal, ["escenario", "metodo", "bloque", "metabolito"]).to_csv(
        SALIDA / "por_bloque.csv", index=False)
    aporte_pareado(pred[pred.diseno.isin(["principal", "mismas_opciones_A"])]).to_csv(
        SALIDA / "aporte_por_grupo.csv", index=False)
    contrastes(pred)
    principal.groupby(["escenario", "metodo", "metabolito", "motivo"]).size().rename(
        "n").to_csv(SALIDA / "motivos.csv")

    indices = ["escenario", "metodo", "celda"]
    antiguo = pd.read_csv(RAIZ / "evaluacion/validacion_celdas.csv").set_index(indices)
    fijo = pred[pred.diseno.isin(["mM_fijo", "ratio_fijo"])]
    tabla = fijo.pivot(index=indices, columns="diseno", values="pred_mM")
    assert tabla.index.sort_values().equals(antiguo.index.sort_values())
    tabla["historico"] = antiguo.pred_mM
    tabla["diferencia_escala_mM"] = (tabla.mM_fijo - tabla.ratio_fijo).abs()
    tabla["diferencia_archivo_mM"] = (tabla.ratio_fijo - tabla.historico).abs()
    tabla["cambio_cobertura_escala"] = tabla.mM_fijo.isna() != tabla.ratio_fijo.isna()
    tabla["cambio_cobertura_archivo"] = tabla.ratio_fijo.isna() != tabla.historico.isna()
    tabla.reset_index().groupby(indices[:2]).agg(
        n_pares=("diferencia_escala_mM", "count"), max_delta_escala_mM=("diferencia_escala_mM", "max"),
        max_delta_archivo_mM=("diferencia_archivo_mM", "max"),
        cambios_cobertura_escala=("cambio_cobertura_escala", "sum"),
        cambios_cobertura_archivo=("cambio_cobertura_archivo", "sum"),
    ).to_csv(SALIDA / "equivalencia_predicciones.csv")


def aplicar(ish, mmc, grupos, ref, validacion, elegidos):
    q = [(s, m) for s in ish.index for m in ref.index if pd.isna(ish.at[s, m])]
    huecos = pd.DataFrame(q, columns=["muestra", "metabolito"])
    huecos["celda"] = huecos.muestra + "|" + huecos.metabolito
    ext = validacion[validacion.diseno.eq("principal") & validacion.pred_mM.notna()]
    soporte = ext.groupby(["escenario", "metodo", "metabolito"]).agg(
        medidas=("celda", "size"), bloques=("bloque", "nunique"))
    e = elegibilidad(ish, grupos)
    resultados = []
    for esc, train in [("A", ish), ("B", pd.concat([ish, mmc]))]:
        for metodo in CONFIG:
            k = elegidos[esc, metodo, "mM"]
            p, info = predecir(train, ish, huecos, e.index[e.elegible], metodo, k)
            datos = soporte.loc[esc, metodo]
            p["n_validacion"] = p.metabolito.map(datos.medidas).fillna(0).astype(int)
            p["bloques_validacion"] = p.metabolito.map(datos.bloques).fillna(0).astype(int)
            sin_prueba = (p.n_validacion < 10) | (p.bloques_validacion < 3)
            p.loc[sin_prueba & p.pred_mM.notna(), "motivo"] = "validacion_externa_insuficiente"
            p.loc[sin_prueba, "pred_mM"] = np.nan
            p["aviso_revision"] = ""
            if info["aviso_convergencia"]:
                p.loc[p.pred_mM.notna(), "aviso_revision"] = "ajuste_final_sin_convergencia"
            resultados.append(p.assign(escenario=esc, metodo=metodo, k=k, **info))
    p = pd.concat(resultados, ignore_index=True)
    p.to_csv(SALIDA / "imputaciones.csv", index=False)
    comun = p.pivot(index="celda", columns=["escenario", "metodo"], values="pred_mM").notna().all(axis=1)
    estad = []
    for (esc, metodo), g in p.groupby(["escenario", "metodo"]):
        original = ish.loc[:, ref.index]
        for alcance in ["observados", "disponibles", "comunes_seis"]:
            matriz = original.copy()
            usar = g[g.pred_mM.notna()]
            if alcance == "comunes_seis":
                usar = usar[usar.celda.map(comun)]
            if alcance != "observados":
                for f in usar.itertuples():
                    assert pd.isna(original.at[f.muestra, f.metabolito])
                    matriz.at[f.muestra, f.metabolito] = f.pred_mM
            pd.testing.assert_frame_equal(matriz.where(original.notna()), original)
            if alcance == "disponibles":
                matriz.to_csv(SALIDA / f"matriz_{esc}_{metodo}_mM.csv")
            for nombre, x in matriz.items():
                n = original[nombre].count()
                media, sd = x.mean(), x.std()
                cv_valido = x.count() >= 2 and pd.notna(media) and media > 0
                estad.append(dict(escenario=esc, metodo=metodo, alcance=alcance, metabolito=nombre,
                                  observados=n, imputados=x.count()-n, pendientes=x.isna().sum(),
                                  media=media, mediana=x.median(), SD=sd,
                                  IQR=x.quantile(.75)-x.quantile(.25),
                                  CV_pct=100*sd/media if cv_valido else np.nan,
                                  CV_interpretable=cv_valido))
    pd.DataFrame(estad).to_csv(SALIDA / "descriptivos.csv", index=False)
