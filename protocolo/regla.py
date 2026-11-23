"""Decisiones internas por objetivo y perfil; comparación siempre pareada por cultivo."""
import numpy as np
import pandas as pd

BASE = "A_amplio/Mediana"
FIJAS = ["A_amplio/Mediana","A_amplio/KNN3","A_amplio/Iterativo"]


def preparar_pares(tabla):
    # Una misma medición puede aparecer en varias máscaras: solo tiene un voto en ese perfil.
    indices = ["cultivo","grupo","fold"]
    x = tabla.copy()
    x["error"] = (x.pred-x.verdad).abs()
    x["clave"] = x.opcion+"/"+x.metodo
    return x.pivot_table(index=indices,columns="clave",values="error",aggfunc="mean")


def comparar(errores, candidato, controles, minimo=10):
    cols = list(dict.fromkeys([candidato,*controles]))
    if any(c not in errores for c in cols) or len(cols)<2:
        return None
    par = errores[cols].dropna()
    n, gr, bl = len(par),par.index.get_level_values("grupo").nunique(),par.index.get_level_values("fold").nunique()
    if n<minimo or gr<5 or bl<2:
        return None
    por_bloque = par.groupby(level="fold").mean()
    otros = [c for c in controles if c!=candidato]
    if not all(par[candidato].mean()<par[c].mean()-1e-12 and
               (por_bloque[candidato]<por_bloque[c]-1e-12).sum()>=2 for c in otros):
        return None
    return n,gr,bl,float(par[candidato].mean())


def aprender(datos):
    if datos.empty:
        return {}
    reglas = {}
    for (y,perfil),g in datos.groupby(["objetivo","perfil"]):
        e = preparar_pares(g)
        if BASE not in e:
            continue
        candidatos = [x for x in ("A_amplio/KNN3","A_amplio/KNN5","A_amplio/Iterativo") if x in e]
        # Exigir exactamente las mismas celdas para ordenar las alternativas básicas.
        ganadores = []
        candidatas = [c for c in candidatos if len(e[[BASE,c]].dropna())>=10]
        ee = e[[BASE,*candidatas]].dropna()
        for c in [BASE,*candidatas]:
            par = comparar(ee,c,[v for v in [BASE,*candidatas] if v!=c])
            if par:
                ganadores.append((par[-1],c,par))
        if ganadores:
            _,elegida,datos_base=min(ganadores)
            tipo="validada"
        else:
            elegida,datos_base,tipo=BASE,(0,0,0,None),"respaldo_general"
        # Una ampliación solo entra si mejora la regla Ishii y sus controles pertinentes.
        extra=[]
        evidencia={}
        for c in e:
            op,met=c.split("/")
            if op in ("A_amplio","A_comun","WT_sensibilidad") or op.startswith("ctrl_") or met=="Mediana":
                continue
            controles=[elegida,BASE]
            if op.startswith("B_comun"):
                controles.append("A_comun/"+met)
            if op!="B_comun":
                controles.append("ctrl_"+op.replace("_ciega_","_red_")+"/"+met)
                controles.append(("B_comun" if op.startswith("B_") else "A_amplio")+"/"+met)
            par=comparar(e,c,list(dict.fromkeys(controles)))
            if par:
                extra.append(c)
                evidencia[c]=par
        # Si los ganadores de distintos ensayos tienen diferente soporte, deben competir otra vez.
        if len(extra)>1:
            finalistas=[c for c in extra if comparar(e,c,[v for v in extra if v!=c])]
        else:
            finalistas=extra
        elegido=finalistas[0] if len(finalistas)==1 else elegida
        par=evidencia[elegido] if elegido!=elegida else datos_base
        reglas[y+"|"+perfil]=dict(opcion=elegido,estado="validada" if elegido!=elegida else tipo,
                                   base=elegida,estado_base=tipo,n_interno=par[0],grupos=par[1],
                                   bloques=par[2],MAE_interno_mM=par[3])
    return reglas


def decidir(evento, alternativas, reglas):
    r=reglas.get(evento.objetivo+"|"+evento.perfil,{})
    opciones=[(r.get("opcion"),r.get("estado"))]
    depende_aux = (r.get("opcion") or "").split("/")[0].count("_") > 1
    if r.get("estado_base")=="validada" or not depende_aux:
        opciones.append((r.get("base"),r.get("estado_base")))
    if not depende_aux:
        opciones.append((BASE,"respaldo_general"))
    for clave,estado in opciones:
        if clave in alternativas and np.isfinite(alternativas[clave].pred):
            f=alternativas[clave]
            return clave,estado,f
    return "", "pendiente", None
