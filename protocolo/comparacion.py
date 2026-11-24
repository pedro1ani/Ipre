"""Evaluación externa pareada. El límite ideal solo existe cuando se conoce la medición."""
from itertools import combinations

import numpy as np
import pandas as pd

from .datos import PRIV, SEMILLA
from .estrategias import DET, RES, FIJOS, PRINCIPALES, guardar


def publicar(tabla, nombre):
    tabla.to_csv(RES / (nombre + ".csv"), index=False, float_format="%.9g")


def metricas(z):
    t = z[z.pred.notna()]
    e = t.pred - t.verdad
    a = e.abs()
    iqr = t.IQR.where(t.IQR > 1e-12)
    return dict(intentos=len(z), predichas=len(t), abstenciones=len(z)-len(t),
                cobertura=len(t)/len(z) if len(z) else np.nan,
                cultivos=t.cultivo.nunique(), grupos=t.grupo.nunique(), objetivos=t.objetivo.nunique(),
                MAE_micro_mM=a.mean(), MAE_macro_mM=a.groupby(t.objetivo).mean().mean(),
                RMSE_mM=np.sqrt((e*e).mean()), sesgo_mM=e.mean(), medAE_mM=a.median(),
                P90_mM=a.quantile(.90), P95_mM=a.quantile(.95), P99_mM=a.quantile(.99),
                MAE_IQR=(a/iqr).mean(), n_IQR=iqr.count(), SD_error_mM=e.std(),
                SD_MAE_bloques_mM=a.groupby(t.fold).mean().std(),
                SD_MAE_grupos_mM=a.groupby(t.grupo).mean().std())


def bootstrap(d, grupos):
    """Re-muestrea perturbaciones completas; la unidad de la media sigue siendo una consulta."""
    a = pd.DataFrame({"d": d, "g": grupos}).dropna()
    if a.g.nunique() < 2:
        return np.nan, np.nan
    sumas = a.groupby("g",sort=True).d.agg(["sum","size"])
    rng = np.random.default_rng(SEMILLA)
    idx = rng.integers(len(sumas), size=(2000, len(sumas)))
    replica = sumas["sum"].to_numpy()[idx].sum(axis=1) / sumas["size"].to_numpy()[idx].sum(axis=1)
    return tuple(np.quantile(replica, [.025, .975]))


def amplio(t):
    proy = t.pivot(index=["fold","evento"], columns="estrategia", values="pred")
    info = t[t.estrategia.eq("Original")].set_index(["fold","evento"])
    return info.drop(columns=["pred","estrategia"]).join(proy, validate="one_to_one")


def ideal(w):
    a = w[FIJOS].sub(w.verdad, axis=0).abs()
    r = w.copy()
    r["pred"] = np.nan
    tiene = a.notna().any(axis=1)
    mejor = a.loc[tiene].idxmin(axis=1)
    r.loc[tiene, "pred"] = w.loc[tiene, FIJOS].to_numpy()[np.arange(len(mejor)),
                                                         [FIJOS.index(m) for m in mejor]]
    r["estrategia"] = "Ideal_no_aplicable"
    return r.reset_index()


def diagnosticar(w, etiqueta):
    w = w.copy()
    w["diferencia"] = (w.Original-w.verdad).abs() - (w.KNN3-w.verdad).abs()
    origen = pd.read_csv(PRIV / "externas.csv.gz", usecols=["evento","fold","clave"])
    w = w.reset_index().merge(origen, on=["evento","fold"], validate="one_to_one")
    w["resultado"] = np.select([w.diferencia.gt(1e-12), w.diferencia.lt(-1e-12)],
                                ["peor_que_KNN3","mejor_que_KNN3"], default="empate")
    w["metodo"] = w.clave.str.partition("/")[2]
    filas = []
    for campo in ["estado", "metodo", "clave", "objetivo", "cultivo", "grupo", "perfil", "resultado", "estado_resultado"]:
        if campo == "estado_resultado":
            w[campo] = w.estado + "/" + w.resultado
        for valor, g in w.groupby(campo):
            delta = g.diferencia.sum() / len(w)
            filas.append(dict(prueba=etiqueta, dimension=campo, valor=valor, n=len(g),
                              MAE_regla_mM=(g.Original-g.verdad).abs().mean(),
                              MAE_KNN3_mM=(g.KNN3-g.verdad).abs().mean(),
                              diferencia_condicional_mM=g.diferencia.mean(), aporte_al_MAE_mM=delta,
                              fraccion_del_exceso=delta/w.diferencia.mean() if abs(w.diferencia.mean())>1e-12 else np.nan))
    for m in FIJOS:
        w["error_"+m] = (w[m]-w.verdad).abs()
    w["error_original"] = (w.Original-w.verdad).abs()
    return pd.DataFrame(filas), w


def ejecutar():
    todo = pd.read_csv(DET/"validacion.csv.gz", float_precision="round_trip", low_memory=False)
    todo["prueba"] = np.where(todo.fold.str.startswith("E"), "perturbacion", "serie")
    resultados, pares, diag, detalles, cond, por_y, frecuencia, motivos = [], [], [], [], [], [], [], []
    for (cv,tipo), z in todo.groupby(["prueba","tipo"]):
        w = amplio(z)
        ids = w.index[w[PRINCIPALES].notna().all(axis=1)]
        evaluar = pd.concat([z, ideal(w)[z.columns]], ignore_index=True)
        for met, t in evaluar.groupby("estrategia"):
            for alcance, g in [("disponibles", t), ("comunes", t[t.set_index(["fold","evento"]).index.isin(ids)])]:
                fila = metricas(g)
                resultados.append(dict(prueba=cv,tipo=tipo,alcance=alcance,estrategia=met,**fila))
                if alcance == "comunes" and met in PRINCIPALES:
                    for y, gy in g.groupby("objetivo"):
                        a = (gy.pred-gy.verdad)
                        por_y.append(dict(prueba=cv,tipo=tipo,estrategia=met,objetivo=y,n=len(a),
                                grupos=gy.grupo.nunique(), MAE_mM=a.abs().mean(),
                                RMSE_mM=np.sqrt(a.pow(2).mean()),sesgo_mM=a.mean()))
            if met in PRINCIPALES:
                f = t.groupby(["estado","nivel","metodo"],dropna=False).size().rename("n").reset_index()
                frecuencia.append(f.assign(estrategia=met,prueba=cv,tipo=tipo))
        for a,b in combinations(PRINCIPALES, 2):
            comparables = w.loc[w[[a,b]].notna().all(axis=1)]
            d = (comparables[a]-comparables.verdad).abs()-(comparables[b]-comparables.verdad).abs()
            lo,hi = bootstrap(d.to_numpy(),comparables.grupo.to_numpy())
            pares.append(dict(prueba=cv,tipo=tipo,a=a,b=b,n=len(d),grupos=comparables.grupo.nunique(),
                              delta_MAE_mM=d.mean(),IC95_inf=lo,IC95_sup=hi))
        if len(ids):
            dg, detalle = diagnosticar(w.loc[ids], cv+"/"+tipo)
            diag.append(dg)
            detalles.append(detalle.assign(prueba=cv))
            for (estado,met), g in detalle.groupby(["estado","metodo"]):
                for alt in ["original", *FIJOS]:
                    cond.append(dict(prueba=cv,tipo=tipo,estado=estado,seleccion=met,alternativa=alt,
                             n=len(g), MAE_mM=g["error_"+alt].mean()))
        # Registrar también la falta de predicción, incluso si no hubo ajuste para el objetivo.
        for fold, t in z[z.estrategia.isin(FIJOS)].groupby("fold"):
            p = pd.read_csv(PRIV/"particiones"/fold/"pred.csv.gz")
            p = p[p.opcion.eq("A_amplio")][["evento","metodo","motivo"]]
            r = t[t.pred.isna()].merge(p,left_on=["evento","estrategia"],right_on=["evento","metodo"],
                                       how="left",suffixes=("","_ajuste"))
            motivos.append(r.assign(motivo=r.motivo.fillna("objetivo_sin_ajuste"),
                                    prueba=cv)[["prueba","tipo","estrategia","motivo"]])
    publicar(pd.DataFrame(resultados), "comparacion_estrategias")
    publicar(pd.DataFrame(pares), "diferencias_pareadas")
    publicar(pd.concat(diag), "descomposicion_regla")
    publicar(pd.DataFrame(cond), "alternativas_mismas_casillas")
    publicar(pd.concat(frecuencia), "frecuencia_eleccion")
    publicar(pd.concat(motivos).value_counts().rename("n").reset_index(), "abstenciones_fijas")
    pd.DataFrame(por_y).to_csv(RES/"error_por_metabolito.csv.gz",index=False,
                              compression={"method":"gzip","mtime":0},float_format="%.9g")
    guardar(pd.concat(detalles), "alternativas_por_casilla.csv.gz")
    # Una fila por objetivo, incluidos los que no tienen verdades comparables.
    from .datos import leer_matriz
    t = todo[todo.prueba.eq("perturbacion") & todo.tipo.eq("individual") & todo.estrategia.eq("KNN3")]
    a = t.assign(error=(t.pred-t.verdad).abs()).groupby("objetivo")
    evidencia = a.agg(intentos=("evento","size"), predichas=("pred","count"), MAE_KNN3_mM=("error","mean"))
    disponible = t[t.pred.notna()].groupby("objetivo")
    evidencia["grupos"] = disponible.grupo.nunique()
    evidencia["bloques"] = disponible.fold.nunique()
    evidencia = evidencia.reindex(leer_matriz("ishii").columns)
    evidencia["apoyo_evaluable_KNN3"] = (evidencia.predichas.ge(10) & evidencia.grupos.ge(5) & evidencia.bloques.ge(2))
    # No equivale a una garantía para el perfil concreto de un hueco real.
    diagnostico_detallado(todo)
    publicar(evidencia.rename_axis("objetivo").reset_index(), "evidencia_por_metabolito")


def diagnostico_detallado(todo):
    """Relaciona el soporte del perfil, el método realmente utilizado y sus alternativas."""
    from .datos import leer_matriz
    ish = leer_matriz("ishii")
    sop = pd.read_csv(DET/"soporte_interno.csv.gz").rename(columns={"particion":"fold"})
    origen = pd.read_csv(PRIV/"externas.csv.gz")
    ori = origen.merge(sop,on=["fold","objetivo","perfil"],how="left",validate="many_to_one")
    ori["soporte_cuatro"] = ori.n.ge(10) & ori.grupos.ge(5) & ori.bloques.ge(2)
    ori["prueba"] = np.where(ori.fold.str.startswith("E"),"perturbacion","serie")
    tabla = ori.groupby(["prueba","tipo","estado","soporte_cuatro"]).size().rename("intentos").reset_index()
    publicar(tabla,"soporte_regla_original")
    parejas=[]
    for fold,g in origen.groupby("fold"):
        p=pd.read_csv(PRIV/"particiones"/fold/"pred.csv.gz",float_precision="round_trip")
        z=g[["evento","objetivo","cultivo","grupo","perfil","pred","verdad","clave","estado","tipo"]].dropna(subset=["pred"])
        z=z.rename(columns={"pred":"elegida"}).merge(p.dropna(subset=["pred"]),on="evento")
        z["error_elegida"]=(z.elegida-z.verdad).abs()
        z["error_alternativa"]=(z.pred-z.verdad).abs()
        z["delta_mM"]=z.error_elegida-z.error_alternativa
        parejas.append(z.assign(fold=fold))
    guardar(pd.concat(parejas,ignore_index=True),"todas_las_alternativas_pareadas.csv.gz")
    z=todo[todo.estrategia.isin(PRINCIPALES)].copy()
    clave=origen.set_index(["fold","evento"]).clave
    z["opcion"]="A_amplio"
    usa_regla=z.estrategia.eq("Original") | (z.estrategia.str.startswith("Respaldo") & z.estado.eq("validada"))
    seleccion=clave.reindex(pd.MultiIndex.from_frame(z.loc[usa_regla,["fold","evento"]]))
    z.loc[usa_regla,"opcion"]=seleccion.str.partition("/")[0].to_numpy()
    z=z[z.prueba.eq("perturbacion") & z.tipo.ne("completa") & z.pred.notna()].drop_duplicates(
        ["estrategia","cultivo","objetivo","perfil","opcion","metodo"])
    llave=["estrategia","objetivo","perfil","opcion","metodo"]
    s=z.assign(error=(z.pred-z.verdad).abs()).groupby(llave).agg(
        n_evaluadas=("cultivo","nunique"),grupos_comparables=("grupo","nunique"),
        bloques_comparables=("fold","nunique"),MAE_comparable_mM=("error","mean")).reset_index()
    h=pd.read_csv(DET/"huecos.csv.gz",float_precision="round_trip")
    h=h[h.estrategia.isin(PRINCIPALES)].copy()
    a=pd.read_csv(PRIV/"auditoria_huecos.csv").set_index("evento")
    h["opcion"]="A_amplio"
    usa=h.estrategia.eq("Original") | (h.estrategia.str.startswith("Respaldo") & h.estado.eq("validada"))
    h.loc[usa,"opcion"]=h.loc[usa,"evento"].map(a.opcion)
    h=h.merge(s,on=llave,how="left",validate="many_to_one")
    h[["n_evaluadas","grupos_comparables","bloques_comparables"]]=h[["n_evaluadas","grupos_comparables","bloques_comparables"]].fillna(0).astype(int)
    apoyo=h.n_evaluadas.ge(10) & h.grupos_comparables.ge(5) & h.bloques_comparables.ge(2)
    h["evidencia_en_prueba"]=np.select([h.pred.isna(),apoyo,h.n_evaluadas.gt(0)],
                            ["pendiente","comparable","escasa"],default="sin_casos")
    pp = pd.read_csv(PRIV/"particiones/real/pred.csv.gz")
    motivos = pp[["evento","opcion","metodo","motivo"]]
    h = h.merge(motivos, on=["evento","opcion","metodo"], how="left", validate="many_to_one")
    h.loc[h.estrategia.eq("Original"),"motivo"] = h.loc[h.estrategia.eq("Original"),"evento"].map(a.motivo)
    for met in FIJOS:
        faltan = h.estrategia.eq(met) & h.pred.isna()
        por_evento = pp[pp.opcion.eq("A_amplio") & pp.metodo.eq(met)].set_index("evento").motivo
        h.loc[faltan,"motivo"] = h.loc[faltan,"evento"].map(por_evento).fillna("objetivo_sin_ajuste")
    h["motivo"] = h.motivo.fillna("sin_alternativa_aplicable")
    guardar(h,"evidencia_por_hueco.csv.gz")
    columnas = ["cultivo","objetivo","perfil","estrategia","opcion","metodo","nivel","estado","pred",
                "evidencia_en_prueba","n_evaluadas","grupos_comparables","bloques_comparables","MAE_comparable_mM",
                "motivo","predictores","vecinos_disponibles","ausencia_original","source_sheet","source_cell"]
    h[columnas].to_csv(RES/"estimaciones_por_casilla.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    publicar(h.groupby(["objetivo","estrategia","estado","evidencia_en_prueba"]).size().rename("casillas").reset_index(),
             "evidencia_huecos_por_objetivo")
    for nombre,g in h.groupby("estrategia"):
        propuesta=g.pivot(index="cultivo",columns="objetivo",values="pred").reindex_like(ish)
        assert not (propuesta.notna() & ish.notna()).any().any()
        destino=DET/"matrices";destino.mkdir(exist_ok=True)
        ish.fillna(propuesta).to_csv(destino/(nombre+".csv"))
    individual=todo[todo.prueba.eq("perturbacion") & todo.tipo.eq("individual") & todo.estrategia.isin(PRINCIPALES)]
    w=amplio(individual).dropna(subset=PRINCIPALES)
    filas=[]
    for m in PRINCIPALES:
        err=w[m]-w.verdad
        for q,g in w.groupby(pd.qcut(w.verdad,5,labels=False,duplicates="drop"),observed=True):
            e=err.loc[g.index]
            filas.append(dict(estrategia=m,quintil=int(q)+1,n=len(e),mediana_verdad_mM=g.verdad.median(),
                MAE_mM=e.abs().mean(),sesgo_mM=e.mean(),SD_error_mM=e.std()))
    publicar(pd.DataFrame(filas),"error_por_concentracion")



if __name__ == "__main__":
    ejecutar()

