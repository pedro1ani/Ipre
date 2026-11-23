"""Métricas sobre verdades de Ishii. Disponibilidad y cobertura de intervalos tienen bases distintas."""
import json

import numpy as np
import pandas as pd

from .datos import PRIV, PUB, guardar_json, huella, leer_matriz
from .regla import FIJAS
from .validar import leer_split


def resumen(t):
    v=t.dropna(subset=["pred"]).copy()
    err=v.pred-v.verdad
    absoluto=err.abs()
    grup=absoluto.groupby(v.grupo).mean()
    rel=absoluto.divide(v.IQR_Ishii_train.where(v.IQR_Ishii_train.gt(0)))
    return dict(intentos=len(t),predicciones=len(v),disponibilidad=len(v)/len(t) if len(t) else np.nan,
                celdas_distintas=len(v[["cultivo","objetivo"]].drop_duplicates()),
                grupos=v.grupo.nunique(),MAE_mM=absoluto.mean(),RMSE_mM=np.sqrt((err**2).mean()),
                sesgo_mM=err.mean(),MAE_sobre_IQR=rel.mean(),n_IQR=rel.count(),
                SD_MAE_grupos_mM=grup.std(),MAE_peor_grupo_mM=grup.max())


def contrastes(q,p,ish,meta,nom):
    """Todas las diferencias usan las mismas consultas y verdades, no diferencias de promedios marginales."""
    if not nom.startswith("E"):
        return []
    p=p[p.pred.notna()]
    filas=[]
    for (op,met),v in p.groupby(["opcion","metodo"]):
        if op.startswith("ctrl_") or op in ("A_amplio","A_comun"):
            continue
        if op=="B_comun":
            control="A_comun"
        elif op=="WT_sensibilidad":
            control="A_amplio"
        else:
            control="ctrl_"+op.replace("_ciega_","_red_")
        a=p[p.opcion.eq(control)&p.metodo.eq(met)][["evento","pred"]].rename(columns={"pred":"control"})
        z=q.merge(v[["evento","pred"]],on="evento").merge(a,on="evento")
        z["opcion"],z["base"],z["metodo"]=op,control,met
        filas.append(z)
    return filas


def escribir(t,nombre,indice=False):
    t.to_csv(PUB/nombre,index=indice,float_format="%.8g")


def metricas_validacion(usar_intervalos=False):
    ish,meta=leer_matriz("ishii"),leer_matriz("metadatos")
    externo=pd.read_csv(PRIV/"externas.csv.gz")
    diseno=json.loads((PRIV/"diseno.json").read_text())
    filas,pares,counts=[],[],[]
    for fold,z in externo.groupby("fold"):
        q,p=leer_split(fold,ish,meta)
        tipo_cv="perturbacion" if fold.startswith("E") else "serie"
        pares.extend(contrastes(q,p,ish,meta,fold))
        for tipo in ("individual","conjunta_anclas","completa"):
            selec=z[z.tipo.eq(tipo)].copy()
            pp=p[p.evento.isin(selec.evento)]
            counts.append((tipo_cv,fold,tipo,len(selec),selec.objetivo.nunique(),selec.grupo.nunique()))
            if selec.empty:
                continue
            fijo=pp.assign(clave_fija=pp.opcion+"/"+pp.metodo)
            fijo=fijo[fijo.clave_fija.isin(FIJAS)].pivot(index="evento",columns="clave_fija",values="pred")
            matriz=selec.set_index("evento")[["cultivo","grupo","objetivo","verdad","IQR_Ishii_train","pred"]].rename(
                columns={"pred":"Regla"}).join(fijo.reindex(columns=FIJAS))
            for alcance in ("disponible","mismas_celdas"):
                m=matriz if alcance=="disponible" else matriz.dropna(subset=["Regla",*FIJAS])
                for metodo in ("Regla",*FIJAS):
                    t=m[["cultivo","grupo","objetivo","verdad","IQR_Ishii_train",metodo]].rename(columns={metodo:"pred"})
                    filas.append(t.assign(fold=fold,prueba=tipo_cv,tipo=tipo,alcance=alcance,metodo=metodo))
    escribir(pd.DataFrame(counts,columns=["prueba","fold","tipo","eventos","objetivos","grupos"]),"mascaras_resumen.csv")
    todos=pd.concat(filas,ignore_index=True)
    columnas=["prueba","tipo","alcance","metodo"]
    por=[]
    for k,t in todos.groupby(columnas+["objetivo"]):
        por.append(dict(zip(columnas+["objetivo"],k))|resumen(t))
    por=pd.DataFrame(por)
    por.to_csv(PRIV/"metricas_por_objetivo.csv.gz",index=False)
    principal=por[por.prueba.eq("perturbacion")&por.tipo.eq("individual")&por.alcance.eq("mismas_celdas")]
    escribir(principal.drop(columns=["prueba","tipo","alcance"]),"por_objetivo.csv")
    globales=[]
    for k,t in todos.groupby(columnas):
        met=por.set_index(columnas).loc[k]
        globales.append(dict(zip(columnas,k))|resumen(t)|dict(MAE_macro_mM=met.MAE_mM.mean(),
            MAE_macro_sobre_IQR=met.MAE_sobre_IQR.mean(),metabolitos=met.MAE_mM.count()))
    escribir(pd.DataFrame(globales),"metricas.csv")
    comp=pd.concat(pares,ignore_index=True)
    comp["error"]=(comp.pred-comp.verdad).abs()
    comp["error_control"]=(comp.control-comp.verdad).abs()
    comp["sesgo"]=comp.pred-comp.verdad
    comp["sesgo_control"]=comp.control-comp.verdad
    llave=["tipo","opcion","base","metodo","objetivo"]
    detalle=comp.groupby(llave).agg(celdas=("evento","size"),grupos=("grupo","nunique"),
        MAE_mM=("error","mean"),MAE_control_mM=("error_control","mean"),sesgo_mM=("sesgo","mean"),
        sesgo_control_mM=("sesgo_control","mean"))
    por_grupo=comp.groupby(llave+["grupo"])[["error","error_control"]].mean()
    bien=por_grupo.error.lt(por_grupo.error_control-1e-12).groupby(level=llave).sum()
    detalle["grupos_con_menor_MAE"]=bien
    detalle["delta_MAE_mM"]=detalle.MAE_mM-detalle.MAE_control_mM
    detalle.reset_index().to_csv(PRIV/"contrastes_pareados.csv.gz",index=False)
    escribir(detalle.reset_index().query("tipo=='individual' and opcion=='B_comun'"),"MMC2_pareado.csv")
    escribir(comp.query("tipo=='individual'").groupby(["opcion","base","metodo"]).agg(
        celdas=("evento","size"),objetivos=("objetivo","nunique"),grupos=("grupo","nunique"),
        MAE_mM=("error","mean"),MAE_control_mM=("error_control","mean"),sesgo_mM=("sesgo","mean")).reset_index(),
        "auxiliares_pareados.csv")
    z=externo.assign(prueba=np.where(externo.fold.str.startswith("E"),"perturbacion","serie"))
    escribir(z.groupby(["prueba","tipo","estado","clave"],dropna=False).size().rename("eventos").reset_index(),"decisiones.csv")
    ajustes=pd.concat([pd.read_csv(PRIV/"particiones"/n/"ajustes.csv").assign(fold=n)
                      for n in externo.fold.unique()],ignore_index=True)
    ajustes["prueba"]=np.where(ajustes.fold.str.startswith("E"),"perturbacion","serie")
    escribir(ajustes.groupby(["prueba","opcion"]).agg(ajustes=("objetivo","size"),
             con_alerta=("alerta_convergencia","sum"),rondas_medianas=("rondas","median")).reset_index(),
             "convergencia.csv")
    if usar_intervalos:
        manifiesto=json.loads((PRIV/"sorteos.json").read_text())
        if (manifiesto["congelado"]!=huella(PRIV/"congelado.json") or
                manifiesto["intervalos"]!=huella(PRIV/"intervalos.csv.gz")):
            raise ValueError("Los intervalos no corresponden a la validación actual.")
        i=pd.read_csv(PRIV/"intervalos.csv.gz")
        e=z.merge(i,on=["fold","evento","objetivo"],how="left",validate="one_to_one")
        e["adentro"]=(e.verdad>=e.lo_mM)&(e.verdad<=e.hi_mM)
        e["ancho_mM"]=e.hi_mM-e.lo_mM
        def cobertura(v):
            ok=v.lo_mM.notna(); n=int(ok.sum())
            return dict(candidatas=len(v),con_prediccion=v.pred.count(),con_intervalo=n,
                        dentro=int(v.adentro.sum()),grupos_con_intervalo=v.loc[ok,"grupo"].nunique(),
                        disponibilidad_intervalo=n/len(v),cobertura_entre_intervalos=v.loc[ok,"adentro"].mean(),
                        ancho_mediano_mM=v.ancho_mM.median())
        for claves,archivo in [(["prueba","tipo"],"intervalos.csv"),(["prueba","tipo","objetivo"],None)]:
            vals=pd.DataFrame([dict(zip(claves,k))|cobertura(v) for k,v in e.groupby(claves)])
            if archivo:
                escribir(vals,archivo)
            else:
                vals.to_csv(PRIV/"cobertura_por_objetivo.csv.gz",index=False)


def informe():
    metricas_validacion(usar_intervalos=True)
    ish,meta=leer_matriz("ishii"),leer_matriz("metadatos")
    diseno=json.loads((PRIV/"diseno.json").read_text())
    a=pd.read_csv(PRIV/"auditoria_huecos.csv")
    inv=pd.read_csv(PUB/"inventario.csv")
    resumen_ejecucion=dict(filas_originales=len(inv),con_diez=int(inv.n_30.ge(10).sum()),objetivos=ish.shape[1],
        cultivos=len(ish),grupos=meta.loc[ish.index,"grupo"].nunique(),medidos=int(ish.count().sum()),
        huecos=len(a),por_estado={str(k):int(v) for k,v in a.estado.value_counts().items()},
        con_intervalo_real=int(a.get("lo_mM",pd.Series(dtype=float)).count()),
        version=diseno["versiones"],semilla=diseno["semilla"])
    guardar_json(PUB/"ejecucion.json",resumen_ejecucion)
    print(json.dumps(resumen_ejecucion,ensure_ascii=False,indent=2))


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--solo-validacion",action="store_true")
    if parser.parse_args().solo_validacion:
        metricas_validacion()
    else:
        informe()
