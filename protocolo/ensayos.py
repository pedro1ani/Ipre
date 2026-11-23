"""Ajustes por partición; cada alternativa conserva sus controles sobre las mismas filas."""

import numpy as np
import pandas as pd

from ampliacion.modelos import elegibilidad
from ampliacion.proteinas import PARES, elegir_sin_red
from .datos import PRIV, PUB, leer_matriz
from .diseno import consultas, escoger_companeros
from .modelos import estimar, matriz_consultas
from sklearn.impute import SimpleImputer

ARCHIVOS = ("ishii","mmc2","wt","metadatos","auxiliares")


def cargar():
    return {k:leer_matriz(k) for k in ARCHIVOS}


def variantes(d, train, y, eleg, comunes_eleg, ciegas, con_wt=False):
    ish, mmc, cov, wt = (d[k] for k in ("ishii","mmc2","auxiliares","wt"))
    a = ish.loc[train]
    fuente_genes = [(g,t) for g,t,*_ in PARES if t == y]
    for panel, campos in [("amplio",eleg),("comun",comunes_eleg)]:
        if y not in campos:
            continue
        cols = [y] + escoger_companeros(a, campos, y)
        if len(cols)<3:
            continue
        for fuente in (("A",) if panel=="amplio" else ("A","B")):
            base = a[cols]
            if fuente=="B":
                base = pd.concat([base, mmc[cols]])
            prefijo = f"{fuente}_{panel}"
            yield prefijo, base, cols, None, None
            if con_wt and prefijo=="A_amplio":
                nuevos = wt[~d["metadatos"].loc[wt.index,"grupo"].isin(
                    d["metadatos"].loc[ish.index.difference(train),"grupo"]).to_numpy()]
                if len(nuevos):
                    yield "WT_sensibilidad", pd.concat([base,nuevos[cols]]),cols,None,None
            if prefijo=="A_comun":
                continue
            columnas = [("qC", "q_C", cov.q_C.notna())]
            if panel=="amplio":
                for gene, objetivo in fuente_genes:
                    elegida = ciegas[gene,objetivo]
                    if elegida:
                        red, sin = "prot:"+gene,"prot:"+elegida
                        par = cov[red].notna() & cov[sin].notna()
                        columnas += [(f"red_{gene}",red,par),(f"ciega_{gene}",sin,par)]
            for nombre, clave, par in columnas:
                x = base.loc[par.reindex(base.index,fill_value=False)]
                if (x[y].count()<10 or d["metadatos"].loc[x.index[x[y].notna()],"grupo"].nunique()<5
                        or any(x[c].nunique()<2 or x[[y,c]].dropna().shape[0]<5 for c in cols[1:])
                        or cov.loc[x.index,clave].nunique()<2):
                    continue
                control = "ctrl_"+prefijo+"_"+(nombre.replace("ciega_","red_"))
                if not nombre.startswith("ciega_"):
                    yield control,x,cols,clave,False
                yield prefijo+"_"+nombre,x,cols,clave,True


def evaluar_split(nombre, ids, prueba, sensibilidad=False, real=False):
    d = cargar()
    ish, meta = d["ishii"], d["metadatos"]
    a = ish.loc[ids]
    comun = ish.columns.intersection(d["mmc2"].columns,sort=False)
    full = elegibilidad(a,meta.grupo)
    corto = elegibilidad(a[comun],meta.grupo)
    eleg = full.index[full.elegible]
    comp = corto.index[corto.elegible]
    q = consultas(a,ish.loc[prueba],real=real)
    salida = PRIV/"particiones"/nombre
    salida.mkdir(parents=True,exist_ok=True)
    q.to_csv(salida/"consultas.csv",index=False)
    grupo = meta.grupo
    prot = d["auxiliares"].loc[ids].filter(like="prot:").rename(columns=lambda c:c[5:])
    ciegas = {(g,t):elegir_sin_red(prot,g,t)[0] for g,t,*_ in PARES}
    resultados, ajustes = [], []
    for y, qq in q.groupby("objetivo",sort=False):
        suficiente = full.at[y,"n"] >= 10 and full.at[y,"grupos"] >= 5 and a[y].nunique() >= 2
        if suficiente and (y not in eleg or len(escoger_companeros(a,eleg,y)) < 2):
            # El objetivo admite referencia univariada aunque no tenga covariables defendibles.
            mediana=SimpleImputer(strategy="median").fit(a[[y]]).statistics_[0]
            resultados.append(pd.DataFrame(dict(evento=qq.evento,opcion="A_amplio",metodo="Mediana",
                  pred=mediana,motivo="general_entrenamiento",visibles=0,donantes=0)))
            ajustes.append((y,"A_amplio","|".join(ids),y,"",False,int(a[y].count()),
                            int(full.at[y,"grupos"]),0,False))
        for nom, x, cols, clave, agregar in variantes(d,ids,y,eleg,comp,ciegas,sensibilidad):
            xi = x.copy()
            auxiliar = d["auxiliares"][clave] if clave else None
            if agregar:
                xi["aux"] = auxiliar.reindex(x.index)
            v = matriz_consultas(ish,qq,cols,auxiliar if agregar else None)
            clave_valida = auxiliar.reindex(qq.cultivo).notna().to_numpy() if clave else np.ones(len(qq),bool)
            if clave and nom.startswith(("ctrl_A_amplio_red_", "A_amplio_red_", "A_amplio_ciega_")):
                gen = nom.rsplit("_",1)[-1]
                par = d["auxiliares"][["prot:"+gen,"prot:"+ciegas[gen,y]]].notna().all(axis=1)
                clave_valida &= par.reindex(qq.cultivo).to_numpy()
            if not clave_valida.any():
                continue
            fit = estimar(xi.to_numpy(),v,len(cols))
            ajustes.append((y,nom,"|".join(x.index),"|".join(cols),clave or "",bool(agregar),
                            int(x[y].count()),int(grupo.loc[x.index[x[y].notna()]].nunique()),
                            fit["Iterativo"][4],fit["Iterativo"][5]))
            for met,(pred,mot,vis,don,_,_) in fit.items():
                pred[~clave_valida],mot[~clave_valida] = np.nan,"sin_auxiliar"
                resultados.append(pd.DataFrame(dict(evento=qq.evento,opcion=nom,metodo=met,
                                                      pred=pred,motivo=mot,visibles=vis,donantes=don)))
    columnas = ["objetivo","opcion","ids","columnas","auxiliar","agregar","observados","grupos",
                "rondas","alerta_convergencia"]
    pd.DataFrame(ajustes,columns=columnas).to_csv(salida/"ajustes.csv",index=False)
    pred = pd.concat(resultados,ignore_index=True) if resultados else pd.DataFrame(
        columns=["evento","opcion","metodo","pred","motivo","visibles","donantes"])
    pred.to_csv(salida/"pred.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    return nombre,len(q)


def diagnostico(splits):
    h = pd.read_csv(PUB/"inventario.csv")
    x = leer_matriz("original_34_todas")[h.loc[h.identidad.isin(["nominal_Ishii","compartida"]),"variable"]]
    meta = leer_matriz("metadatos")
    filas = []
    for nom,(train,_) in splits.items():
        for n in (8,10,12):
            e = elegibilidad(x.loc[train],meta.grupo,n)
            filas.append((nom,n,int(e.elegible.sum())))
    pd.DataFrame(filas,columns=["particion","minimo","n_elegibles"]).to_csv(PUB/"umbrales.csv",index=False)
