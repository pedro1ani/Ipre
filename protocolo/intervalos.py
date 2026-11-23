"""Percentiles predictivos por celda para la opción iterativa realmente utilizada."""
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from .datos import PRIV, SEMILLA, guardar_json, huella
from .ensayos import cargar
from .modelos import iterador, matriz_consultas
from .validar import NPROC


def simular(tarea):
    fold,y,op,qs,aj=tarea
    with threadpool_limits(limits=1):
        d=cargar()
        cols=aj["columnas"].split("|")
        tr=pd.concat([d["ishii"],d["mmc2"].reindex(columns=d["ishii"].columns)])
        tr=tr.loc[aj["ids"].split("|"),cols].copy()
        aux=None
        if aj["agregar"]:
            aux=d["auxiliares"][aj["auxiliar"]]
            tr["aux"]=aux.reindex(tr.index)
        v=matriz_consultas(d["ishii"],qs,cols,aux)
        esc=StandardScaler().fit(tr.to_numpy())
        z,query=esc.transform(tr.to_numpy()),esc.transform(v)
        # Estos son sorteos independientes de la imputación, no el IC de la media.
        n=100 if "_red_" in op or "_ciega_" in op else 200
        muestra=np.empty((n,len(qs)))
        sin_faltantes=np.isfinite(z).all()
        for i in range(n):
            semilla=SEMILLA+1000+i
            if not sin_faltantes or i==0:
                modelo,_=iterador(z,esc,posterior=True,semilla=semilla)
            else:
                # Con entrenamiento completo el ajuste no consume sorteos; equivale al reajuste.
                modelo.random_state_=np.random.RandomState(semilla)
            muestra[i]=modelo.transform(query)[:,0]*esc.scale_[0]+esc.mean_[0]
        assert np.isfinite(muestra).all() and muestra.min()>=-1e-10
        lo,hi=np.quantile(np.maximum(0,muestra),[.025,.975],axis=0)
        return pd.DataFrame(dict(fold=fold,evento=qs.evento,objetivo=y,lo_mM=lo,hi_mM=hi,sorteos=n))


def intervalos(procesos=NPROC):
    externas=pd.read_csv(PRIV/"externas.csv.gz")
    real=pd.read_csv(PRIV/"auditoria_huecos.csv").assign(fold="real")
    base=pd.concat([externas,real],ignore_index=True)
    selec=base[base.clave.str.endswith("/Iterativo",na=False)]
    etiquetas=["fold","evento","objetivo","lo_mM","hi_mM","sorteos"]
    tareas=[]
    for (fold,y,clave),q in selec.groupby(["fold","objetivo","clave"]):
        op=clave.split("/")[0]
        a=pd.read_csv(PRIV/"particiones"/fold/"ajustes.csv")
        aj=a[a.objetivo.eq(y)&a.opcion.eq(op)].iloc[0].to_dict()
        tareas.append((fold,y,op,q[["evento","cultivo","ocultas"]].copy(),aj))
    print(f"{len(tareas)} ajustes con sorteos; {len(selec)} predicciones iterativas.",flush=True)
    con=[]
    with ProcessPoolExecutor(max_workers=procesos) as pool:
        futuros=[pool.submit(simular,t) for t in tareas]
        for n,f in enumerate(as_completed(futuros),1):
            con.append(f.result())
            if n%30==0:
                print(n,len(tareas),flush=True)
    i=pd.concat(con,ignore_index=True) if con else pd.DataFrame(columns=etiquetas)
    i.sort_values(["fold","evento"]).to_csv(PRIV/"intervalos.csv.gz",index=False,
                                            compression={"method":"gzip","mtime":0})
    r=real.drop(columns=["fold","lo_mM","hi_mM","sorteos"],errors="ignore").merge(
        i[i.fold.eq("real")].drop(columns="fold"),on=["objetivo","evento"],how="left",validate="one_to_one")
    assert r.lo_mM.notna().equals(r.metodo.eq("Iterativo"))
    r.to_csv(PRIV/"auditoria_huecos.csv",index=False)
    guardar_json(PRIV/"sorteos.json",dict(congelado=huella(PRIV/"congelado.json"),
                 intervalos=huella(PRIV/"intervalos.csv.gz"),semilla_inicial=SEMILLA+1000,
                 sorteos_base=200,sorteos_proteinas=100,implementacion=huella(Path(__file__))))



if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--procesos",type=int,default=NPROC)
    intervalos(p.parse_args().procesos)
