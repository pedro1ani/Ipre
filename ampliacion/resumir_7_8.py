"""Comparaciones emparejadas: un voto por celda aun cuando aparece bajo varios patrones."""
import sqlite3
from contextlib import closing

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .medir_bloques import medir
from .preparar import RAIZ

ID = ["bloque","grupo","metabolito","celda","evento"]
TOL = 1e-12


def contrastar(t):
    """Delta negativo significa menor error de la alternativa. El grupo se retira entero."""
    resumen, detalle_grupo, influencia = [], [], []
    for (f,esc,tipo), q in t.groupby(["familia","escenario","tipo"]):
        ancho = q.assign(error=(q.pred_mM-q.verdad_mM).abs()).pivot(
            index=ID,columns=["metodo","alternativa"],values="error").dropna()
        if ancho.empty:
            continue
        celdas = ancho.groupby(level=ID[:-1]).mean()
        comparadores=[(a,"control") for a in q.alternativa.unique() if a!="control"]
        if "sin_red" in set(q.alternativa):
            comparadores.append(("red","sin_red"))
        if "solo_origen" in set(q.alternativa):
            comparadores.append(("q_C","solo_origen"))
        for metodo in q.metodo.unique():
            for alt,base in comparadores:
                d = celdas[metodo,alt]-celdas[metodo,base]
                por_grupo = d.groupby(level="grupo").mean()
                por_met = d.groupby(level="metabolito").mean()
                etiquetas = dict(familia=f,escenario=esc,tipo=tipo,metodo=metodo,alternativa=alt,
                                 comparador=base)
                cambios = []
                for g, val in por_grupo.items():
                    fuera = d[d.index.get_level_values("grupo")!=g]
                    delta = fuera.groupby(level="metabolito").mean().mean()
                    cambios.append(delta)
                    detalle_grupo.append(dict(**etiquetas,grupo=g,delta_MAE_mM=val,
                                             delta_macro_sin_grupo=delta))
                for m, val in por_met.items():
                    influencia.append(dict(**etiquetas,metabolito=m,delta_MAE_mM=val,
                                           aporte_delta_macro=val/len(por_met)))
                s = pd.Series(cambios).dropna()
                bloques = d.groupby(level="bloque").mean()
                resumen.append(dict(**etiquetas,celdas=len(d),grupos=len(por_grupo),bloques=len(bloques),
                    delta_macro_mM=por_met.mean(),delta_micro_mM=d.mean(),
                    grupos_mejor=(por_grupo < -TOL).sum(),grupos_peor=(por_grupo>TOL).sum(),
                    grupos_igual=(por_grupo.abs()<=TOL).sum(),bloques_mejor=(bloques < -TOL).sum(),
                    menor_delta_sin_grupo=s.min(),mayor_delta_sin_grupo=s.max(),
                    suficiente_descriptivo=len(d)>=10 and len(por_grupo)>=5 and len(bloques)>=2))
    return pd.DataFrame(resumen),pd.DataFrame(detalle_grupo),pd.DataFrame(influencia)


def calibrar(intervalos, pred):
    claves = ["familia","tipo","sorteos","nivel","grupo","metabolito","celda","evento"]
    # Sólo Iterative produce intervalos. Se emparejan los tres brazos antes de puntuar.
    pi=intervalos.assign(cubierto=intervalos.verdad_mM.between(intervalos.inferior_mM,intervalos.superior_mM),
                          ancho_mM=intervalos.superior_mM-intervalos.inferior_mM,
                          sesgo_media_mM=intervalos.media_posterior_mM-intervalos.verdad_mM)
    wide=pi.pivot(index=claves,columns="alternativa",values="inferior_mM")
    assert set(wide.columns)=={"control","sin_red","red"}
    pares=wide.dropna().index
    intentos=pred[pred.familia.isin(pi.familia.unique()) & pred.metodo.eq("Iterative")].groupby(
        ["familia","tipo","alternativa"]).agg(intentos=("evento","size"),celdas_intentadas=("celda","nunique"))
    partes=[]
    for nombre,mask in [("disponibles",np.ones(len(pi),dtype=bool)),
                        ("comunes",pd.MultiIndex.from_frame(pi[claves]).isin(pares))]:
        v=pi[mask]
        idx=["familia","tipo","sorteos","nivel","alternativa"]
        g=v.groupby(idx+["grupo","celda"])
        c=g[["cubierto","ancho_mM","sesgo_media_mM"]].mean().assign(eventos=g.size(),aciertos=g.cubierto.sum()).reset_index()
        sub=c.groupby(idx)
        s=sub.agg(eventos_con_intervalo=("eventos","sum"),celdas=("celda","size"),grupos=("grupo","nunique"),
                  aciertos_evento=("aciertos","sum"),cobertura=("cubierto","mean"),ancho_mM=("ancho_mM","mean"),
                  sesgo_posterior_mM=("sesgo_media_mM","mean")).reset_index()
        # Mantener también combinaciones con cero intervalos.
        combinaciones=pd.MultiIndex.from_product([sorted(pi.sorteos.unique()),sorted(pi.nivel.unique())],
                                                 names=["sorteos","nivel"]).to_frame(index=False)
        s=intentos.reset_index().merge(combinaciones,how="cross").merge(s,on=idx,how="left")
        for col in ["eventos_con_intervalo","celdas","grupos","aciertos_evento"]:
            s[col]=s[col].fillna(0).astype(int)
        s["disponibilidad"]=s.eventos_con_intervalo/s.intentos
        partes.append(s.assign(alcance=nombre))
    return pd.concat(partes,ignore_index=True)


def anotar_ciegas(destino):
    t=pd.read_csv(destino/"detalle/proteinas_elegidas.csv")
    with closing(sqlite3.connect(f"file:{RAIZ/'sources/ecoli_auditada_v3.db'}?mode=ro",uri=True)) as db:
        red=pd.read_sql_query("""SELECT m.protein_id_original AS proteina, r.reaction_id AS reaccion, r.formula
            FROM ishii_mmc4_mapping m JOIN reactions r ON m.functional_annotation_reaction_id=r.reaction_id""",db)
    def relacionada(r):
        sub=red[red.proteina.eq(r.sin_red)]
        return "+".join(s.reaccion for s in sub.itertuples() if r.objetivo+"[c]" in s.formula.split())
    t["coincidencia_directa_documentada"]=[relacionada(r) for r in t.itertuples()]
    t.groupby(["par","sin_red","coincidencia_directa_documentada"],dropna=False).agg(
        bloques=("bloque","nunique"),max_diferencia_disponibilidad=("diferencia_patron_train","max"),
        minimo_cultivos_comunes_train=("n_train","min")).reset_index().to_csv(destino/"controles_sin_red.csv",index=False)


def tablas_puntuales(pred,destino):
    detalle=destino/"detalle"
    glob,met,gr=medir(pred)
    glob.to_csv(detalle/"metricas_globales.csv.gz",index=False)
    met.to_csv(detalle/"error_por_metabolito.csv.gz",index=False)
    gr.to_csv(detalle/"error_por_grupo_y_metabolito.csv.gz",index=False)
    a,b,c=contrastar(pred)
    a.to_csv(detalle/"estabilidad_contrastes.csv.gz",index=False)
    b.to_csv(detalle/"efecto_y_exclusion_por_grupo.csv.gz",index=False)
    c.to_csv(detalle/"aporte_por_metabolito.csv.gz",index=False)
    base=glob[glob.alcance.eq("comunes")].copy()
    base[["familia","fuente","tipo"]]=base.experimento.str.split("__",expand=True)
    idx=["familia","escenario","tipo","metodo","intentos","eventos_predichos","celdas_predichas","grupos_evaluados"]
    fis=base[~base.alternativa.isin(["red","sin_red"]) & base.familia.isin(["captacion","flujo","D_principal","D_con_WT"])]
    fis=fis.pivot(index=idx,columns="alternativa",values="MAE_mM_macro").reset_index()
    fis.to_csv(destino/"tasas_resultados.csv",index=False)
    prot=base[base.familia.isin(pred[pred.alternativa.eq("red")].familia.unique())]
    prot=prot.pivot(index=idx,columns="alternativa",values="MAE_mM_macro").reset_index()
    prot.to_csv(destino/"proteinas_resultados.csv",index=False)
    a=a[a.tipo.eq("individual") & a.metodo.ne("Mediana")]
    a.to_csv(destino/"estabilidad_individual.csv",index=False)
    pred.groupby(["familia","escenario","tipo","alternativa","metodo","motivo"]).size().rename(
        "eventos").to_csv(detalle/"motivos.csv")
    anotar_ciegas(destino)
    familias=pred.loc[pred.alternativa.eq("red"),"familia"].unique()
    lectura=base[base.familia.isin(familias) & base.tipo.eq("individual")]
    campos=["familia","metodo","alternativa","intentos","celdas_predichas","grupos_evaluados",
            "MAE_mM_micro","RMSE_mM_micro","sesgo_mM_micro","MAE_IQR_macro"]
    lectura[campos].to_csv(destino/"proteinas_errores.csv",index=False)
    foco=met.query("experimento=='captacion__B__individual' and metodo=='Iterative' and alcance=='comunes' and celdas_predichas>0")
    foco[["alternativa","metabolito","celdas_predichas","grupos_evaluados","MAE_mM","RMSE_mM",
          "sesgo_mM"]].to_csv(destino/"captacion_por_metabolito.csv",index=False)


def figura(destino):
    fig,(e1,e2)=plt.subplots(1,2,figsize=(10.6,4.5),layout="constrained")
    colores={"control":"#404a4c","sin_red":"#236877","red":"#bc552b"}
    prot=pd.read_csv(destino/"proteinas_resultados.csv").query("tipo=='individual' and metodo=='Iterative'")
    for j,r in enumerate(prot.itertuples()):
        cambios=[100*(getattr(r,s)/r.control-1) for s in ("sin_red","red")]
        e1.plot(cambios,[j,j],color="0.8",zorder=0)
        for k,s in enumerate(("sin_red","red")):
            e1.scatter(cambios[k],j,marker=["o","D"][k],color=colores[s],s=35,
                       label=["Por disponibilidad","Por reacción"][k] if j==0 else None)
    e1.axvline(0,linestyle="--",color=".5",linewidth=1)
    e1.set(yticks=range(len(prot)),yticklabels=prot.familia.str.replace("_"," / "),
           xlabel="Cambio de MAE respecto al control (%)",title="Proteínas · iterativo, solo Ishii")
    e1.invert_yaxis()
    e1.legend(fontsize=8)
    cap=pd.read_csv(destino/"tasas_resultados.csv").query("tipo=='individual' and familia=='captacion' and metodo=='Iterative'").set_index("escenario")
    for i,s in enumerate(["A","B"]):
        e2.plot([cap.at[s,"control"],cap.at[s,"q_C"]],[i,i],color="0.75",lw=2,zorder=0)
        for k,nom in enumerate(["control","q_C"]):
            e2.scatter(cap.at[s,nom],i,color=["#404a4c","#236877"][k],marker=["o","D"][k],
                       s=48,label=["Solo metabolitos","+ Captación de carbono"][k] if i==0 else None)
    e2.set(yticks=[0,1],yticklabels=["Ishii","Ishii + MMC2"],xlabel="MAE macro (mM)",
           title="Captación · iterativo, 375 celdas",ylim=(-.6,1.65))
    e2.legend(fontsize=8,loc="upper right")
    for ax in (e1,e2):
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="x",alpha=.2)
    fig.suptitle("Ocultamiento individual · menor error es mejor",fontsize=13,weight="bold")
    with plt.rc_context({"svg.fonttype":"none","svg.hashsalt":"IPre78"}):
        fig.savefig(destino/"tasas_y_proteinas.svg",metadata={"Date":None})
    fig.savefig(destino/"tasas_y_proteinas.png",dpi=160)
    plt.close(fig)


def resumir(pred,intervalos,destino):
    tablas_puntuales(pred,destino)
    ca=calibrar(intervalos,pred)
    ca.to_csv(destino/"detalle/calibracion_todos.csv.gz",index=False)
    ca.query("alcance=='comunes' and sorteos==100 and tipo=='individual'").to_csv(
        destino/"intervalos_proteinas.csv",index=False)
    from .proteinas import sensibilidad_rf03
    sensibilidad_rf03(destino)
    figura(destino)
