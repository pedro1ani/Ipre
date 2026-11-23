"""Aplicación posterior a la validación: lee y verifica la regla congelada."""
import json

import numpy as np
import pandas as pd

from .datos import PRIV, PUB, guardar_json, huella, leer_matriz
from .ensayos import cargar, evaluar_split
from .modelos import matriz_consultas
from sklearn.metrics.pairwise import nan_euclidean_distances
from sklearn.preprocessing import StandardScaler
from .regla import FIJAS
from .validar import huellas_entradas, leer_split, seleccionar, versiones


def comprobar():
    cierre=json.loads((PRIV/"congelado.json").read_text())
    if huellas_entradas()!=cierre["entradas"]:
        raise ValueError("Cambió la matriz, su procedencia o los auxiliares: repetir preparación y validación.")
    for nombre,clave in [("reglas.json","reglas"),("diseno.json","diseno"),("externas.csv.gz","externas")]:
        if huella(PRIV/nombre)!=cierre[clave]:
            raise ValueError(f"{nombre} no coincide con la validación congelada.")
    diseno=json.loads((PRIV/"diseno.json").read_text())
    if versiones()!=diseno["versiones"] or any(huella(PUB.parent/n)!=v for n,v in diseno["codigo"].items()):
        raise ValueError("Cambió el código o entorno de validación: volver a validar.")
    if any(huella(PRIV/"particiones"/n/"consultas.csv")!=v for n,v in cierre["mascaras"].items()):
        raise ValueError("Las máscaras no corresponden a esta validación.")
    return json.loads((PRIV/"reglas.json").read_text())["F"]


def comparar_descriptivos(original, matrices):
    vacios=original.isna()
    comun=vacios.copy()
    for x in matrices.values():
        comun &= x.notna()
    filas=[]
    for alcance in ("disponible","huecos_comunes"):
        for metodo,m in {"observado":original,**matrices}.items():
            x=m if alcance=="disponible" else original.where(~comun,m)
            for y,v in x.items():
                medidos=original[y].count()
                media,sd=v.mean(),v.std(ddof=1)
                cv=100*sd/media if v.count()>=2 and media>0 else np.nan
                filas.append((alcance,metodo,y,medidos,int(v.count()-medidos),int(v.isna().sum()),
                              media,v.median(),sd,v.quantile(.75)-v.quantile(.25),cv,
                              "" if np.isfinite(cv) else "menos_de_2_o_media_no_positiva"))
    return pd.DataFrame(filas,columns=["alcance","metodo","objetivo","observados","imputados","pendientes",
                                      "media_mM","mediana_mM","SD_mM","IQR_mM","CV_pct","motivo_CV"])


def identificar_donantes(z):
    """Reconstruye el conjunto de vecinos efectivos y verifica su media contra scikit-learn."""
    d=cargar()
    datos=pd.concat([d["ishii"],d["mmc2"].reindex(columns=d["ishii"].columns)])
    z["ids_donantes_KNN"]=""
    for (op,y),q in z[z.metodo.str.startswith("KNN",na=False)].groupby(["opcion","objetivo"]):
        f=q.iloc[0]
        cols=f.columnas.split("|")
        train=datos.loc[f.cultivos_entrenamiento.split("|"),cols].copy()
        aux=d["auxiliares"][f.auxiliar] if f.agregar else None
        v=matriz_consultas(d["ishii"],q,cols,aux)
        if aux is not None:
            train["aux"]=aux.reindex(train.index)
        escala=StandardScaler().fit(train.to_numpy())
        normal=escala.transform(train.to_numpy())
        for j,r in enumerate(q.itertuples()):
            apoyo=(np.isfinite(v[j,1:len(cols)]) & np.isfinite(normal[:,1:len(cols)])).sum(axis=1)
            idx=np.where((apoyo>=2) & np.isfinite(normal[:,0]))[0]
            dist=nan_euclidean_distances(escala.transform(v[[j]]),normal[idx])[0]
            k=min(int(r.metodo[-1]),len(idx))
            ids=train.index[idx[np.argpartition(dist,k-1)[:k]]]
            np.testing.assert_allclose(train.loc[ids,y].mean(),r.pred,rtol=1e-10,atol=1e-12)
            assert len(ids)==r.donantes
            z.at[r.Index,"ids_donantes_KNN"]="|".join(ids)
    return z


def aplicar():
    reglas=comprobar()
    ish,meta=leer_matriz("ishii"),leer_matriz("metadatos")
    evaluar_split("real",ish.index,ish.index,real=True)
    q,p=leer_split("real",ish,meta,con_verdad=False)
    z=seleccionar(q,p,reglas)
    ausencia=leer_matriz("ausencias")
    z["ausencia_original"]=[ausencia.at[s,y] for s,y in zip(z.cultivo,z.objetivo)]
    permitido=z.ausencia_original.eq("not_detected_unspecified")
    z.loc[~permitido,["clave","estado","motivo","pred"]]=["","pendiente","ausencia_no_autorizada",np.nan]
    z["lectura"]="exploratoria: no detectable o bajo límite sin distinción en Ishii"
    z["perfil_completo_cultivo"]=["|".join(ish.columns[ish.loc[s].isna()]) for s in z.cultivo]
    z["unidad"]="mM"
    z["opcion"] = z.clave.str.partition("/")[0]
    z["metodo"] = z.clave.str.partition("/")[2]
    a=pd.read_csv(PRIV/"particiones/real/ajustes.csv").rename(columns={"ids":"cultivos_entrenamiento"})
    z=z.merge(a,on=["opcion","objetivo"],how="left",validate="many_to_one")
    trazas=pd.read_csv(PRIV/"trazabilidad.csv.gz",low_memory=False)
    z=z.merge(trazas[["sample_key","variable","source_sheet","source_cell","analyte_id","protocol"]].rename(
        columns={"sample_key":"cultivo","variable":"objetivo"}),on=["cultivo","objetivo"],validate="one_to_one")
    ex=pd.read_csv(PRIV/"externas.csv.gz")
    ex=ex[ex.fold.str.startswith("E") & ex.tipo.ne("completa") & ex.pred.notna()].copy()
    ex=ex.drop_duplicates(["cultivo","objetivo","perfil","clave"])
    ex["error"]=(ex.pred-ex.verdad).abs()
    comp=ex.groupby(["objetivo","perfil","clave"]).agg(MAE_externo_comparable_mM=("error","mean"),
              n_celdas_comparables=("error","size"),n_grupos_comparables=("grupo","nunique")).reset_index()
    z=z.merge(comp,on=["objetivo","perfil","clave"],how="left",validate="many_to_one")
    z["n_celdas_comparables"]=z.n_celdas_comparables.fillna(0).astype(int)
    z["n_grupos_comparables"]=z.n_grupos_comparables.fillna(0).astype(int)
    z["evidencia_externa"]=np.select([z.n_celdas_comparables.eq(0),
        z.n_celdas_comparables.lt(10)|z.n_grupos_comparables.lt(5)],["sin_casos","escasa"],default="comparable")
    z=identificar_donantes(z)
    z.loc[z.metodo.ne("Iterativo"),["rondas","alerta_convergencia"]]=None
    z.to_csv(PRIV/"auditoria_huecos.csv",index=False)
    matrices={}
    for nom in ["regla",*FIJAS]:
        if nom=="regla":
            valores=z[["cultivo","objetivo","pred"]]
        else:
            opcion,metodo=nom.split("/")
            valores=q.loc[permitido.to_numpy(),["evento","cultivo","objetivo"]].merge(
                p[p.opcion.eq(opcion)&p.metodo.eq(metodo)][["evento","pred"]],on="evento",how="left")
        estimados=valores.pivot(index="cultivo",columns="objetivo",values="pred").reindex_like(ish)
        m=ish.fillna(estimados)
        assert np.array_equal(m.to_numpy()[ish.notna()],ish.to_numpy()[ish.notna()])
        matrices[nom]=m
        m.to_csv(PRIV/("matriz_"+nom.replace("/","_")+".csv"))
    completa=leer_matriz("original_34_todas").loc[ish.index].copy()
    completa.to_csv(PRIV/"original_30_completa.csv")
    completa.loc[:,ish.columns]=matrices["regla"]
    completa.to_csv(PRIV/"resultado_30_completo.csv")
    solo=ish.copy()
    for f in z[z.estado.eq("validada")].itertuples():
        solo.at[f.cultivo,f.objetivo]=f.pred
    solo.to_csv(PRIV/"matriz_solo_con_eleccion_validada.csv")
    comparar_descriptivos(ish,matrices).to_csv(PUB/"descriptivos.csv",index=False,float_format="%.8g")
    z.groupby(["objetivo","estado","opcion","metodo"],dropna=False).size().rename("huecos").to_csv(
        PUB/"aplicacion_por_objetivo.csv")
    guardar_json(PRIV/"aplicacion.json",dict(reglas=huella(PRIV/"reglas.json"),
                 candidatos=len(z),originales_intactos=int(ish.count().sum()),
                 estimaciones=int(z.pred.count()),pendientes=int(z.pred.isna().sum())))
    print((PRIV/"aplicacion.json").read_text(),flush=True)


if __name__=="__main__":
    aplicar()
