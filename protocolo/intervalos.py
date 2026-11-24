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
    fold,y,op,qs,aj,*destino=tarea
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
        muestra = np.maximum(0, muestra)
        carpeta = PRIV / (destino[0] if destino else "posteriores") / fold
        carpeta.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(carpeta / f"{y}__{op}.npz",
                            eventos=qs.evento.to_numpy(dtype=str), valores=muestra)
        resultado = pd.DataFrame(dict(fold=fold, evento=qs.evento, objetivo=y, sorteos=n,
                                       sd_posterior_mM=muestra.std(axis=0, ddof=1)))
        for nivel in (50, 80, 90, 95):
            cola = (100 - nivel) / 200
            lo, hi = np.quantile(muestra, [cola, 1 - cola], axis=0)
            sufijo = "mM" if nivel == 95 else str(nivel)
            resultado["lo_" + sufijo] = lo
            resultado["hi_" + sufijo] = hi
        return resultado


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
    r=real.drop(columns=[c for c in i if c not in ("evento","objetivo")],errors="ignore").merge(
        i[i.fold.eq("real")].drop(columns="fold"),on=["objetivo","evento"],how="left",validate="one_to_one")
    assert r.lo_mM.notna().equals(r.metodo.eq("Iterativo"))
    r.to_csv(PRIV/"auditoria_huecos.csv",index=False)
    guardar_json(PRIV/"sorteos.json",dict(congelado=huella(PRIV/"congelado.json"),
                 intervalos=huella(PRIV/"intervalos.csv.gz"),semilla_inicial=SEMILLA+1000,
                 sorteos_base=200,sorteos_proteinas=100,implementacion=huella(Path(__file__))))



def intervalos_iterativo_fijo(procesos=NPROC):
    """Control sin selección de método: todos los casos individuales por perturbación y huecos aplicables."""
    from .datos import leer_matriz
    from .validar import leer_split

    ish, meta = leer_matriz("ishii"), leer_matriz("metadatos")
    tareas, consultas = [], []
    for fold in [*(f"E{i}" for i in range(5)), "real"]:
        q, p = leer_split(fold, ish, meta, con_verdad=fold!="real")
        q = q[q.tipo.isin(["individual","real"])]
        if fold == "real":
            a = pd.read_csv(PRIV/"auditoria_huecos.csv")
            q = q[q.evento.isin(a.loc[a.ausencia_original.eq("not_detected_unspecified"),"evento"])]
        g = q.merge(p[p.opcion.eq("A_amplio") & p.metodo.eq("Iterativo") & p.pred.notna()],
                    on="evento",validate="one_to_one")
        ajustes = pd.read_csv(PRIV/"particiones"/fold/"ajustes.csv").query("opcion=='A_amplio'").set_index("objetivo")
        consultas.append(g)
        for y,z in g.groupby("objetivo"):
            tareas.append((fold,y,"A_amplio",z[["evento","cultivo","ocultas"]],ajustes.loc[y].to_dict(),"posteriores_fijos"))
    with ProcessPoolExecutor(max_workers=procesos) as pool:
        res = pd.concat(pool.map(simular,tareas),ignore_index=True)
    base = pd.concat(consultas,ignore_index=True).merge(res,on=["fold","evento","objetivo"],validate="one_to_one")
    ruta=PRIV/"revision"
    ruta.mkdir(parents=True,exist_ok=True)
    base.sort_values(["fold","evento"]).to_csv(ruta/"posteriores_iterativo_fijo.csv.gz",index=False,
                                             compression={"method":"gzip","mtime":0})
    # La cobertura se puede medir solo en las casillas con verdad conocida.
    b=base[base.fold.ne("real")]
    filas=[]
    for y,g in [("todos",b),*b.groupby("objetivo")]:
        for nivel in (50,80,90,95):
            suf="mM" if nivel==95 else str(nivel)
            ancho=g["hi_"+suf]-g["lo_"+suf]
            error=(g.pred-g.verdad).abs()
            rho=ancho.corr(error,method="spearman") if len(g)>=8 and ancho.nunique()>1 and error.nunique()>1 else np.nan
            filas.append(dict(objetivo=y,nominal=nivel/100,n=len(g),grupos=g.grupo.nunique(),
                              observada=g.verdad.between(g["lo_"+suf],g["hi_"+suf]).mean(),
                              mediana_ancho_mM=ancho.median(),rho_ancho_error=rho))
    from .estrategias import RES
    RES.mkdir(exist_ok=True,parents=True)
    pd.DataFrame(filas).to_csv(RES/"calibracion_iterativo_fijo.csv",index=False,float_format="%.9g")



if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--procesos",type=int,default=NPROC)
    intervalos(p.parse_args().procesos)
