"""Intervalos de las comparaciones proteicas: reutiliza el mismo ajuste si sus entradas coinciden."""
import hashlib
import inspect
from concurrent.futures import ProcessPoolExecutor
from importlib.metadata import version

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

from .incertidumbre import imputador, perfiles, simular_lote
from .modelos import CONFIG

CLAVE = ["tipo", "evento"]


def entradas(x, ish, q, cols, k, cov):
    medido = x[cols].to_numpy(dtype=float)
    oculto, filas, objetivos = perfiles(ish, q, cols)
    # La proteína no satisface el requisito de dos metabolitos acompañantes.
    assert (objetivos >= 0).all() and (np.isfinite(oculto[filas]).sum(axis=1) >= 2).all()
    if cov is not None and cov.loc[x.index].nunique() > 1:
        muestras = q[["muestra","ocultar"]].drop_duplicates().muestra
        medido = np.column_stack([medido,cov.loc[x.index]])
        oculto = np.column_stack([oculto,cov.loc[muestras]])
        assert np.isfinite(medido[:,-1]).all() and np.isfinite(oculto[:,-1]).all()
    scaler = StandardScaler().fit(medido)
    return (scaler.transform(medido), scaler.transform(oculto), -scaler.mean_/scaler.scale_,
            k, filas, objetivos, scaler.mean_, scaler.scale_)


def identificador(datos):
    h = hashlib.sha256()
    h.update((inspect.getsource(imputador)+inspect.getsource(simular_lote)+str(CONFIG["Iterative"])
              +"|".join(version(p) for p in ["numpy","scipy","scikit-learn"])).encode())
    for d in datos:
        a = np.ascontiguousarray(d)
        h.update(str((a.shape,a.dtype)).encode())
        h.update(a.tobytes())
    return h.hexdigest()


def simular(casos, ish, destino, procesos):
    cache = destino/"detalle/sorteos"
    cache.mkdir(parents=True,exist_ok=True)
    grupos = {}
    for par, brazo, b, x, cols, k, cov, q in casos:
        key = (b,tuple(x.index),tuple(cols),k,None if cov is None else cov.name)
        grupos.setdefault(key,[]).append((par,brazo,x,cov,q))
    tareas, archivos = [], []
    for (b, _, cols, k, _), lista in grupos.items():
        _, _, x, cov, _ = lista[0]
        # Reunir consultas permite ajustar una vez por semilla, nunca usar la prueba para entrenar.
        q = pd.concat([r[-1] for r in lista]).drop_duplicates(CLAVE).sort_values(CLAVE).reset_index(drop=True)
        datos = entradas(x,ish,q,list(cols),k,cov)
        firma = identificador(datos)
        rutas = [cache/f"{firma}_{inicio:04}.npy" for inicio in range(1100,1200,10)]
        tareas.extend((ruta,sem,sem+10,datos) for ruta,sem in zip(rutas,range(1100,1200,10)))
        archivos.append((q,lista,rutas,firma))
    tareas.sort(key=lambda t:t[1])
    print(f"{len(grupos)} entrenamientos distintos × 100 reajustes independientes; {procesos} procesos.",flush=True)
    with ProcessPoolExecutor(max_workers=procesos) as pool:
        for i,_ in enumerate(pool.map(simular_lote,tareas),1):
            if i%100==0 or i==len(tareas):
                print(f"{i}/{len(tareas)} lotes posteriores",flush=True)
    partes = []
    for q,lista,rutas,firma in archivos:
        dibujos = np.concatenate([np.load(p,allow_pickle=False) for p in rutas])
        assert dibujos.shape == (100,len(q))
        for n in [30,100]:
            a = dibujos[:n]
            for nivel in [50,80,95]:
                alpha = (1-nivel/100)/2
                li,ls = np.quantile(a,[alpha,1-alpha],axis=0)
                z = q[CLAVE].assign(inferior_mM=li,superior_mM=ls,media_posterior_mM=a.mean(axis=0),
                    SD_posterior_mM=a.std(axis=0,ddof=1),columna_sorteo=np.arange(len(q)),lote=firma)
                for par,brazo,_,_,preguntas in lista:
                    v = preguntas.merge(z,on=CLAVE,validate="one_to_one")
                    partes.append(v.assign(nivel=nivel,sorteos=n))
    todo = pd.concat(partes,ignore_index=True)
    # La tabla de predicciones conserva el denominador, incluidos los intentos sin proteína o sin soporte.
    todo.to_csv(destino/"detalle/intervalos.csv.gz",index=False)
    return todo, len(grupos)*100
