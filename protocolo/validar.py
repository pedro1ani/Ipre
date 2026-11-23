"""Ejecutar: python -m protocolo.validar. No completa ninguna ausencia real."""
import argparse
import hashlib
import json
import os
import platform
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import scipy
import sklearn
from threadpoolctl import threadpool_limits

from .datos import PRIV, PUB, SEMILLA, guardar_json, huella, leer_matriz
from .diseno import bloques
from .ensayos import diagnostico, evaluar_split
from .regla import aprender, decidir

NPROC = min(7, len(os.sched_getaffinity(0))) if hasattr(os,"sched_getaffinity") else 2


def particiones(ish, meta):
    grupos = meta.grupo
    externas = {f"E{i}": (a,b) for i,(a,b) in enumerate(bloques(ish.index,grupos,5))}
    for serie, prueba in meta.loc[ish.index].groupby("serie_familia"):
        a = ish.index[~grupos.loc[ish.index].isin(prueba.grupo)]
        externas[f"S{int(serie)}"] = a,prueba.index
    splits = dict(externas)
    for nombre, a in [(n,tr) for n,(tr,_) in externas.items()]+[("F",ish.index)]:
        for i,(tr,te) in enumerate(bloques(a,grupos,3)):
            splits[f"{nombre}_i{i}"]=tr,te
    for a,b in splits.values():
        assert not set(grupos.loc[a]) & set(grupos.loc[b])
    return splits,externas


def huellas_entradas():
    return {n:huella(PRIV/n) for n in ("ishii.csv","mmc2.csv","wt.csv","metadatos.csv","ausencias.csv",
                                        "auxiliares.csv","fuentes.json","original_34_todas.csv","trazabilidad.csv.gz")}


def versiones():
    return dict(python=platform.python_version(),pandas=pd.__version__,numpy=np.__version__,
                scipy=scipy.__version__,sklearn=sklearn.__version__)


def leer_split(nombre, ish, meta, con_verdad=True, para_regla=False):
    carpeta=PRIV/"particiones"/nombre
    q=pd.read_csv(carpeta/"consultas.csv")
    p=pd.read_csv(carpeta/"pred.csv.gz",float_precision="round_trip")
    q["grupo"]=q.cultivo.map(meta.grupo)
    q["fold"]=nombre
    if para_regla:
        # Evitar contar dos veces la misma verdad para un perfil; completa es una sensibilidad aparte.
        q=q[q.tipo.ne("completa")].drop_duplicates(["cultivo","objetivo","perfil"])
    if con_verdad:
        q["verdad"]=[ish.at[s,y] for s,y in zip(q.cultivo,q.objetivo)]
        assert q.verdad.notna().all()
    return q,p


def seleccionar(q,p,reglas):
    indices = {evento:g for evento,g in p.groupby("evento",sort=False)}
    filas=[]
    for f in q.itertuples():
        grupo=indices.get(f.evento)
        alts={} if grupo is None else {v.opcion+"/"+v.metodo:v for v in grupo.itertuples()}
        clave,estado,r=decidir(f,alts,reglas)
        regla=reglas.get(f.objetivo+"|"+f.perfil,{})
        primaria=alts.get(regla.get("opcion"))
        if r:
            motivo=r.motivo
        elif primaria:
            motivo=primaria.motivo+";sin_alternativa_validada"
        else:
            motivo="objetivo_o_modelo_sin_soporte_en_entrenamiento"
        filas.append((f.evento,clave,estado,r.pred if r else np.nan,
                      motivo,
                      int(r.visibles) if r else 0,int(r.donantes) if r else 0))
    return q.merge(pd.DataFrame(filas, columns=["evento","clave","estado","pred","motivo","visibles","donantes"]),
                   on="evento",validate="one_to_one")


def ajustar_regla(n, ish, meta):
    data=[]
    for i in range(3):
        q,p=leer_split(f"{n}_i{i}",ish,meta,para_regla=True)
        data.append(q.merge(p[p.pred.notna()],on="evento"))
    return aprender(pd.concat(data,ignore_index=True))


def tarea(args):
    with threadpool_limits(limits=1):
        return evaluar_split(*args)


def validar(procesos=NPROC, reutilizar=False):
    ish,meta=leer_matriz("ishii"),leer_matriz("metadatos")
    splits,externas=particiones(ish,meta)
    fuentes=huellas_entradas()
    modulos={p.name:huella(p) for p in sorted((PUB.parent).glob("*.py")) if p.name in
             {"modelos.py","diseno.py","ensayos.py","datos.py","regla.py","validar.py"}}
    manifest=dict(entradas=fuentes,codigo=modulos,semilla=SEMILLA,versiones=versiones(),
                  particiones={n:dict(train=list(tr),test=list(te)) for n,(tr,te) in splits.items()})
    archivo=PRIV/"diseno.json"
    anterior=json.loads(archivo.read_text()) if archivo.exists() else {}
    def huella_de_ajustes(m):
        # Se puede volver a aprender la decisión sin repetir predicciones idénticas.
        return {k:v for k,v in m.items() if k!="codigo"} | {"modelos":{
            n:v for n,v in m.get("codigo",{}).items() if n not in ("regla.py","validar.py")}}
    if reutilizar and huella_de_ajustes(anterior)!=huella_de_ajustes(manifest):
        raise ValueError("Las entradas o el diseño cambiaron. Ejecuta sin --reanudar.")
    guardar_json(archivo,manifest)
    diagnostico(splits)
    todo=[]
    for n,(tr,te) in splits.items():
        d=PRIV/"particiones"/n
        if reutilizar and (d/"ok.json").exists():
            ok=json.loads((d/"ok.json").read_text())
            if huella(d/"consultas.csv")!=ok["consultas"] or huella(d/"pred.csv.gz")!=ok["pred"]:
                raise ValueError(f"Resultados guardados de {n} no válidos.")
            continue
        todo.append((n,tr,te,n.startswith("E") and "_" not in n))
    # Un proceso por partición. La iteración de cada modelo se limita a un hilo.
    with ProcessPoolExecutor(max_workers=procesos) as pool:
        tareas={pool.submit(tarea,t):t[0] for t in todo}
        for f in as_completed(tareas):
            n,c=f.result()
            carpeta=PRIV/"particiones"/n
            guardar_json(carpeta/"ok.json",{"consultas":huella(carpeta/"consultas.csv"),"eventos":c,
                                            "pred":huella(carpeta/"pred.csv.gz")})
            print(n,c,flush=True)
    reglas={}
    externos=[]
    for n in [*externas,"F"]:
        reglas[n]=ajustar_regla(n,ish,meta)
        if n!="F":
            q,p=leer_split(n,ish,meta)
            z=seleccionar(q,p,reglas[n])
            tr=splits[n][0]
            escala=ish.loc[tr].quantile(.75)-ish.loc[tr].quantile(.25)
            z["IQR_Ishii_train"]=z.objetivo.map(escala)
            externos.append(z)
    seleccion=pd.concat(externos,ignore_index=True)
    seleccion.to_csv(PRIV/"externas.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    mascaras={n:huella(PRIV/"particiones"/n/"consultas.csv") for n in splits}
    guardar_json(PRIV/"reglas.json",reglas)
    guardar_json(PRIV/"congelado.json",dict(entradas=fuentes,reglas=huella(PRIV/"reglas.json"),
                 mascaras=mascaras, diseno=huella(archivo), externas=huella(PRIV/"externas.csv.gz")))
    print("Reglas congeladas. Aún no se completó ningún hueco.",flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--procesos",type=int,default=NPROC)
    p.add_argument("--reanudar",action="store_true",help="solo para entradas y diseño idénticos")
    args=p.parse_args()
    validar(args.procesos,args.reanudar)
