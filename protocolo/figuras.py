"""Figuras de informe. La función por_metabolito permite consultar cualquiera de las 135 señales."""
import io

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from PIL import Image
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .datos import ANCLAS, PRIV, PUB, leer_matriz
from .distribuciones import describir, matrices
from .estrategias import DET, RES, FIJOS, PRINCIPALES

FIG = RES / "figuras"
COLORES = dict(zip(["Observado", *FIJOS, "Original", "Jerarquica"],
                   ["#263746","#a76b19","#18768a","#3953a5","#a14673","#5e6871","#28814f"]))
PANEL = ["fum","pep","s7p","succ","cit","ade","atp","ish_mal","ish_f6p","ish_g6p"]


def lienzo(filas, cols, titulo, pie, ancho=12, alto=7):
    fig, ax = plt.subplots(filas, cols, figsize=(ancho,alto), squeeze=False)
    fig.subplots_adjust(left=.085,right=.965,top=.84,bottom=.15,hspace=.65,wspace=.45)
    fig.text(.065,.955,titulo,fontsize=17,weight="bold",va="top",color="#193444")
    fig.text(.065,.035,pie,fontsize=8.5,color="#40545c")
    return fig, ax.ravel()


def salvar(fig, nombre):
    FIG.mkdir(parents=True,exist_ok=True)
    stream = io.BytesIO()
    fig.savefig(stream,format="png",dpi=150,facecolor="white")
    stream.seek(0)
    with Image.open(stream) as img:
        img.convert("RGB").quantize(colors=256).save(FIG/(nombre+".png"),optimize=True)
    plt.close(fig)


def leer(nombre):
    carpeta = PUB if nombre == "inventario" else RES
    return pd.read_csv(carpeta/(nombre+".csv"))


def ecdf(ax, serie, nombre, color):
    t = np.sort(pd.Series(serie).dropna())
    if len(t):
        ax.step(t,np.arange(1,len(t)+1)/len(t),where="post",label=nombre,color=color,lw=1.6)


def xy(ax,x,y,identidad=False):
    x,y=np.asarray(x),np.asarray(y)
    ok=np.isfinite(x)&np.isfinite(y)
    ax.scatter(x[ok],y[ok],s=12,alpha=.43,color="#18768a",edgecolors="none")
    if identidad and ok.any():
        lo,hi=min(x[ok].min(),y[ok].min()),max(x[ok].max(),y[ok].max())
        ax.plot([lo,hi],[lo,hi],"--",color="#483b35",lw=1)
    return int(ok.sum())


def firma(y):
    nombres = leer("inventario").set_index("variable").original_label
    return str(nombres.get(y,y))+" ("+y+")"


def por_metabolito(y="fum", estrategias=FIJOS, ish=None, completadas=None):
    """Original, completado y estimado se distinguen por color, símbolo y denominador."""
    ish = leer_matriz("ishii") if ish is None else ish
    if y not in ish:
        raise KeyError(f"{y} no está en la tabla de Ishii. Consulte inventario.csv.")
    completadas = dict(matrices(ish)) if completadas is None else completadas
    obs = ish[y].dropna()
    v = {m:completadas[m][y].dropna() for m in estrategias}
    imp = {m:completadas[m].loc[ish[y].isna(),y].dropna() for m in estrategias}
    fig,a=lienzo(2,3,"Distribución de "+firma(y),
        f"{len(obs)} mediciones originales · {ish[y].isna().sum()} huecos · todos los ejes de concentración en mM. "
        "Los totales incluyen las mediciones; los rombos muestran solo estimaciones.",13,8)
    bins=np.histogram_bin_edges(pd.concat([obs,*v.values()]),bins=12)
    for k,(name,s) in enumerate([("Observado",obs),*v.items()]):
        c=COLORES.get(name,f"C{k%10}")
        a[0].hist(s,bins=bins,histtype="step",color=c,lw=1.4,label=f"{name} ({len(s)})")
        ecdf(a[1],s,name,c)
        a[2].boxplot(s,positions=[k],widths=.55,showfliers=False)
        if len(s)>=8 and s.nunique()>1:
            for body in a[3].violinplot(s,positions=[k],showextrema=False)["bodies"]:
                body.set_facecolor(c)
        a[3].scatter(k+np.linspace(-.11,.11,len(obs)),obs,s=7,c=COLORES["Observado"],alpha=.55)
        if name in imp:
            a[3].scatter(k+np.linspace(-.09,.09,len(imp[name])),imp[name],s=14,c=c,marker="D")
            if len(imp[name]):
                a[4].hist(imp[name],bins=bins,histtype="step",lw=1.7,color=c,label=name)
                a[5].scatter([k]*len(imp[name]),imp[name],marker="D",color=c,s=18)
    labs=["Original",*estrategias]
    for ax in (a[2],a[3],a[5]):
        ax.set_xticks(range(len(labs)),labs,rotation=24,ha="right",fontsize=8); ax.set_ylabel("Concentración (mM)")
    for ax in (a[0],a[1],a[4]): ax.set_xlabel("Concentración (mM)")
    for ax in (a[0],a[4]): ax.set_ylabel("Recuento")
    a[1].set_ylabel("Fracción acumulada");a[0].legend(fontsize=7)
    if not any(len(s) for s in imp.values()):
        texto = "Sin huecos en Ishii" if ish[y].notna().all() else "Sin estimaciones aplicables"
        for ax in [a[4],a[5]]:
            ax.text(.5,.55,texto,transform=ax.transAxes,ha="center",fontsize=10,color="#485d63")
    for ax,t in zip(a,["Histograma: originales y totales","Distribución acumulada","Cajas del conjunto disponible",
                       "Violines y puntos por cultivo","Histograma: solo imputados","Estimaciones individuales"]):ax.set_title(t,fontsize=10)
    return fig


def por_metabolito_original(y="fum", ish=None):
    """Consulta previa a la imputación, sin depender de ninguna salida de los modelos."""
    ish=leer_matriz("ishii") if ish is None else ish
    s=ish[y].dropna()
    fig,a=lienzo(1,3,"Mediciones originales de "+firma(y),
                 f"Concentraciones en mM · {len(s)} medidas · {len(ish)-len(s)} ausentes entre {len(ish)} cultivos. "
                 "Cada punto es un cultivo. Ningún valor extremo ha sido eliminado.",11,4.6)
    a[0].hist(s,bins="auto",color="#337a86");a[0].set(xlabel="Concentración (mM)",ylabel="Cultivos",title="Histograma")
    ecdf(a[1],s,"Ishii","#18768a");a[1].set(xlabel="Concentración (mM)",ylabel="Fracción acumulada",title="ECDF")
    if len(s)>=8 and s.nunique()>1:
        a[2].violinplot(s,positions=[1],showextrema=False)
    a[2].boxplot(s,positions=[1],widths=.2,showfliers=False)
    a[2].scatter(1+np.linspace(-.12,.12,len(s)),s,s=15,color="#245161",zorder=3)
    a[2].set(ylabel="Concentración (mM)",xticks=[],title="Caja y observaciones")
    return fig



def original(ish):
    d=describir(ish)
    fig,a=lienzo(2,3,"Qué midió Ishii antes de imputar",
        f"{len(ish)} cultivos · {len(d)} metabolitos · {int(d.n.sum())} mediciones. "
        "Los ejes logarítmicos cambian únicamente la representación.",13,7.8)
    for ax,k,u in [(a[0],"n","Mediciones por metabolito"),(a[1],"CV_pct","CV (%)"),(a[2],"SD","SD (mM)")]:
        t=d[k].dropna()
        bins=np.geomspace(t[t>0].min(),t.max(),20) if k=="SD" else 18
        ax.hist(t,bins=bins,color="#347d87");ax.set(xlabel=u,ylabel="Metabolitos",title=u)
        if k=="SD":ax.set_xscale("log")
    xy(a[3],d.media,d.SD);a[3].set(xscale="log",yscale="log",xlabel="Media (mM)",ylabel="SD (mM)",title="Media y dispersión")
    a[4].imshow(ish.isna(),cmap=ListedColormap(["#cce3df","#183e54"]),aspect="auto",interpolation="nearest")
    a[4].set(xlabel="Metabolito (orden de la matriz)",ylabel="Cultivo",title="Oscuro = ausencia")
    for j,y in enumerate(["fum","cit","atp"]):
        t=ish[y].dropna()
        a[5].scatter(j+np.linspace(-.14,.14,len(t)),t,s=13,label=f"{y}: n={len(t)}")
    a[5].set(yscale="log",ylabel="Concentración (mM)",xticks=range(3),xticklabels=["fum","cit","atp"],title="Mediciones individuales");a[5].legend(fontsize=7)
    salvar(fig,"01_mediciones_y_ausencias")


def estadisticos(ish, comps):
    o=describir(ish)
    fig,a=lienzo(2,4,"Qué cambia al completar los huecos reales",
        "135 metabolitos · 30 cultivos. Diagonal: sin cambio. Las líneas de abajo conservan el orden original de los metabolitos. "
        "Cocientes solo si la SD original es interpretable.",15,8)
    for m in FIJOS:
        d=describir(comps[m]);prod=d.n.gt(o.n)
        for ax,k in zip(a,["media","mediana","SD","CV_pct","IQR"]):
            col=COLORES[m]
            ax.scatter(o[k],d[k],s=7,alpha=.5,color=col,label=m)
            vals=pd.concat([o[k],d[k]]);lo=vals[vals>0].min();hi=vals.max()
            if m==FIJOS[-1]:ax.plot([lo,hi],[lo,hi],"k--",lw=.8)
            ax.set(xscale="symlog",yscale="symlog",xlabel="Original",ylabel="Posterior",title=k+(" (%)" if k=="CV_pct" else " (mM)"))
        ratio=(d.SD/o.SD.where(o.SD>np.maximum(1e-12,ish.abs().max()*1e-6))).where(prod)
        for ax,s,base,title in [(a[5],100*(ratio-1),0,"Cambio de SD (%)"),(a[6],ratio,1,"SD posterior / original")]:
            ax.plot(np.arange(len(s)),s,".",ms=3,color=COLORES[m]);ax.axhline(base,color="k",lw=.6)
            ax.set(xlabel="Metabolito",title=title)
    a[7].axis("off");a[7].legend(*a[0].get_legend_handles_labels(),loc="center",frameon=False)
    salvar(fig,"03_estadisticos_antes_despues")


def relaciones(ish, comps, m="Iterativo"):
    x=comps[m]
    fig,a=lienzo(2,3,"Relaciones entre metabolitos · "+m,
       "Pearson y Spearman sobre pares de cultivos con ambas señales; mínimo 8 pares. Gris = sin soporte. "
       "Después pueden incorporarse cultivos nuevos. Δ = posterior − original.",13,8)
    for i,tipo in enumerate(["pearson","spearman"]):
        o=ish[PANEL].corr(method=tipo,min_periods=8);p=x[PANEL].corr(method=tipo,min_periods=8)
        for ax,z,t in zip(a[3*i:3*i+3],[o,p,p-o],["Medido","Completado","Diferencia"]):
            im=ax.imshow(z,cmap=plt.get_cmap("RdBu_r").with_extremes(bad="#d6d6d6"),vmin=-1 if t!="Diferencia" else -2,
                         vmax=1 if t!="Diferencia" else 2)
            ax.set_xticks(range(len(PANEL)),PANEL,rotation=70,fontsize=6);ax.set_yticks(range(len(PANEL)),PANEL,fontsize=6)
            ax.set_title(t+" · "+tipo.capitalize(),fontsize=10);fig.colorbar(im,ax=ax,shrink=.72)
    salvar(fig,"04_correlaciones")
    fig,a=lienzo(1,3,"Pares de metabolitos y componentes principales",
       "PCA descriptivo sobre cuatro anclas. La media, la SD y los ejes se ajustan una vez, solo con cultivos originalmente completos. "
       "Los símbolos abiertos son cultivos añadidos.",13,5)
    for ax,(u,v) in zip(a,[("fum","succ"),("pep","s7p")]):
        fin=ish[[u,v]].notna().all(axis=1);nuevo=x[[u,v]].notna().all(axis=1)&~fin
        ax.scatter(ish.loc[fin,u],ish.loc[fin,v],s=22,c="#1e697d",label=f"Medidos ({fin.sum()})")
        ax.scatter(x.loc[nuevo,u],x.loc[nuevo,v],marker="D",facecolors="none",edgecolors="#ad5637",label=f"Con imputación ({nuevo.sum()})")
        ax.set(xlabel=u+" (mM)",ylabel=v+" (mM)");ax.legend(fontsize=7)
    ref=ish[list(ANCLAS)].dropna(); scaler=StandardScaler().fit(ref)
    modelo=PCA(n_components=2,svd_solver="full").fit(scaler.transform(ref))
    for tab,lab,mark in [(ref,"Original","o"),(x.loc[x.index.difference(ref.index),list(ANCLAS)].dropna(),"Añadido","D")]:
        if len(tab):
            pcs=modelo.transform(scaler.transform(tab));a[2].scatter(*pcs.T,label=f"{lab} ({len(tab)})",marker=mark,facecolors="none" if mark=="D" else "#1e697d",edgecolors="#ad5637" if mark=="D" else "#1e697d")
    v=modelo.explained_variance_ratio_*100
    a[2].set(xlabel=f"PC1 ({v[0]:.0f} %)",ylabel=f"PC2 ({v[1]:.0f} %)",title="Proyección en los mismos ejes");a[2].legend(fontsize=7)
    salvar(fig,"05_pares_y_pca")


def errores(z):
    comparables=z.pivot(index="evento",columns="estrategia",values="pred")[PRINCIPALES].dropna().index
    t=z[z.evento.isin(comparables)].copy()
    t["residuo"]=t.pred-t.verdad;t["ae"]=t.residuo.abs()
    n=len(comparables)
    fig,a=lienzo(2,3,"Errores sobre concentraciones ocultas",
        f"Perturbación, máscara individual · {n} mediciones comunes. Residuo = predicción − medición. "
        "Escala logarítmica simétrica para mostrar también los extremos.",13,8)
    for k,m in enumerate(FIJOS):
        g=t[t.estrategia.eq(m)];c=COLORES[m]
        for ax,key in [(a[0],"ae"),(a[1],"residuo")]:
            lim=max(t.ae.max(),.001);ext=np.geomspace(.0001,lim,36)
            bins=np.r_[0,ext] if key=="ae" else np.r_[-ext[::-1],0,ext]
            ax.hist(g[key],bins=bins,color=c,label=m,histtype="step",lw=1.3)
            ax.set_xscale("symlog",linthresh=.004);ax.set_ylabel("Mediciones")
        ecdf(a[2],g.ae,m,c)
        a[3].boxplot(g.ae,positions=[k],widths=.5,flierprops=dict(markersize=1))
    a[0].set_title("Error absoluto (mM)");a[0].legend(fontsize=7);a[1].set_title("Error con signo (mM)")
    a[2].set(xlabel="Error absoluto (mM)",ylabel="Fracción acumulada",title="ECDF");a[2].set_xscale("symlog",linthresh=.004)
    a[3].set(xticks=range(4),xticklabels=FIJOS,title="Error absoluto por método",ylabel="mM");a[3].set_yscale("symlog",linthresh=.004)
    g=t[t.estrategia.eq("KNN3")]
    for i,y in enumerate(PANEL[:7]):
        e=g.loc[g.objetivo.eq(y),"ae"]
        if len(e):a[4].boxplot(e,positions=[i],widths=.5,flierprops=dict(markersize=2))
    a[4].set(xticks=range(7),xticklabels=PANEL[:7],title="Error por metabolito · KNN3",ylabel="mM");a[4].set_yscale("symlog",linthresh=.005)
    xy(a[5],g.verdad,g.pred,True);a[5].set(xlabel="Medición (mM)",ylabel="Estimación (mM)",title="KNN3 · predicción frente a medición")
    for ax in [a[5]]:ax.set_xscale("symlog",linthresh=.01);ax.set_yscale("symlog",linthresh=.01)
    salvar(fig,"06_errores_externos")
    fig,a=lienzo(2,4,"Dónde cuesta más reconstruir una medición · KNN3",
        f"{n} casillas de Ishii · grupos de perturbación excluidos por bloque. "
        "SD, faltantes y donantes se calculan solo con entrenamiento; no se descartan extremos.",15,8)
    for ax,k,v,lab in zip(a,["verdad","pred","faltante_train_pct","SD_train","predictores","vecinos_disponibles"],
            ["residuo","residuo","ae","ae","ae","ae"],["Medición (mM)","Predicción (mM)","Faltantes del objetivo (%)","SD del objetivo (mM)","Predictores visibles","Donantes disponibles"]):
        xy(ax,g[k],g[v]);ax.set(xlabel=lab,ylabel=("Residuo (mM)" if v=="residuo" else "Error absoluto (mM)"))
        ax.set_yscale("symlog",linthresh=.007)
        if v=="residuo":ax.axhline(0,color="k",lw=.7)
        if k in ("pred","verdad","SD_train"):ax.set_xscale("symlog",linthresh=.01)
    por=g.groupby("grupo").ae.mean().sort_values();a[6].bar(range(len(por)),por,color="#18768a");a[6].set_xticks(range(len(por)),por.index,fontsize=5,rotation=80)
    a[6].set(title="Grupos experimentales",ylabel="MAE (mM)")
    for m in FIJOS:
        h=t[t.estrategia.eq(m)];bins=pd.qcut(h.verdad,5,duplicates="drop");k=h.groupby(bins,observed=True).agg(x=("verdad","median"),y=("ae","mean"))
        a[7].plot(k.x,k.y,"o-",ms=3,color=COLORES[m],label=m)
    a[7].set(xscale="log",yscale="log",xlabel="Medición por quintil (mM)",ylabel="MAE (mM)");a[7].legend(fontsize=6)
    salvar(fig,"07_residuos_y_soporte")


def decisiones(z):
    g=z[z.estrategia.eq("Original")]
    q=leer("alternativas_mismas_casillas").query("prueba=='perturbacion' and tipo=='individual'")
    c=leer("comparacion_estrategias").query("prueba=='perturbacion' and tipo=='individual'")
    fig,a=lienzo(2,3,"Qué decide la regla y cuánto cuesta esa decisión",
         "Ocultamiento individual por perturbación. Paneles de error condicionados a las mismas mediciones. "
         "El ideal por casilla conoce la verdad: solo es un límite retrospectivo.",14,8.8)
    freq=g.groupby(["metodo","estado"]).size().unstack(fill_value=0)
    freq.plot.bar(stacked=True,ax=a[0],rot=15,fontsize=8);a[0].legend(fontsize=7)
    a[0].set(title=f"Elección original · {len(g)} intentos",ylabel="Mediciones",xlabel="")
    h=z[z.estrategia.isin(PRINCIPALES[4:]) & z.estado.isin(["respaldo","respaldo_general"])]
    h.groupby(["estrategia","metodo"]).size().unstack(fill_value=0).plot.bar(stacked=True,ax=a[1],rot=20,fontsize=6)
    a[1].set(title="Método usado como respaldo",xlabel="",ylabel="Mediciones");a[1].legend(fontsize=6)
    h=q.assign(decision=q.estado.str[:3]+"/"+q.seleccion).pivot(index="decision",columns="alternativa",values="MAE_mM")
    h.plot.bar(ax=a[2],fontsize=6,rot=25);a[2].legend(fontsize=6);a[2].set(title="Dentro de cada decisión original",ylabel="MAE (mM)",xlabel="")
    h=g[g.objetivo.isin(PANEL) & g.metodo.notna()]
    tabla=pd.crosstab([h.objetivo,h.perfil],h.metodo)
    if len(tabla):
        cod=tabla.idxmax(axis=1).map({m:i for i,m in enumerate(FIJOS)}).unstack()
        a[3].imshow(cod,aspect="auto",vmin=-.5,vmax=3.5,cmap=ListedColormap([COLORES[k] for k in FIJOS]).with_extremes(bad="#eeeeee"))
        a[3].set_xticks(range(len(cod.columns)),[s.replace("anclas_presentes","completo") for s in cod.columns],rotation=28,ha="right",fontsize=6)
        a[3].set_yticks(range(len(cod)),cod.index,fontsize=7)
        a[3].legend(handles=[Line2D([],[],color=COLORES[k],lw=4,label=k) for k in FIJOS],fontsize=5,ncols=2)
    a[3].set_title("Método más frecuente por perfil",fontsize=10)
    p=z[z.estrategia.isin(["Original","KNN3"])].pivot(index="evento",columns="estrategia",values="pred").dropna()
    verdad=g.set_index("evento").verdad.reindex(p.index)
    dif=(p.Original-verdad).abs()-(p.KNN3-verdad).abs()
    ecdf(a[4],dif,"Regla − KNN3","#18768a");a[4].axvline(0,ls="--",lw=1,color="k")
    a[4].set_xscale("symlog",linthresh=.002);a[4].set(xlabel="Diferencia de error (mM)",ylabel="Fracción acumulada",title=f"Regla frente a KNN3 · n={len(p)}")
    for i,m in enumerate([*PRINCIPALES,"Ideal_no_aplicable"]):
        ct=c[c.estrategia.eq(m)].set_index("alcance")
        a[5].scatter(ct.at["disponibles","cobertura"]*100,ct.at["comunes","MAE_micro_mM"],s=24,label=f"{i+1}: {m}")
        a[5].annotate(str(i+1),(ct.at["disponibles","cobertura"]*100,ct.at["comunes","MAE_micro_mM"]),fontsize=6)
    a[5].set(xlabel="Cobertura sobre todos los intentos (%)",ylabel="MAE común (mM)",title="Error y cobertura");a[5].legend(fontsize=5,loc="best",ncols=2)
    salvar(fig,"08_eleccion_y_respaldo")


def incertidumbre():
    z=pd.read_csv(DET/"incertidumbre_externa.csv.gz").query("prueba=='perturbacion' and tipo=='individual'")
    fig,a=lienzo(2,4,"Cuánta incertidumbre refleja el intervalo iterativo",
        f"Regla original · {len(z)} intervalos con verdad oculta · 95 % nominal; el ancho se mide en mM. "
        "Los sorteos describen una casilla, no cultivos adicionales ni un IC del MAE.",15,8)
    # Se ilustran tres objetivos prefijados y el primer evento por ID, sin elegirlos por su error.
    ejemplos=[]
    for y in ["fum","pep","atp"]:
        g=z[z.objetivo.eq(y)].sort_values(["fold","evento"])
        if len(g):ejemplos.append(g.iloc[0])
    for i,r in enumerate(ejemplos):
        with archivo_posterior(r) as datos:
            v=datos["valores"][:,np.flatnonzero(datos["eventos"]==r.evento)[0]]
        a[i].hist(v,bins=17,color="#78acba",label=f"{len(v)} sorteos");a[i].axvspan(r.lo_mM,r.hi_mM,color="#91c7a6",alpha=.28,label="Intervalo 95 %")
        a[i].axvline(r.verdad,color="#a24835",ls="--",label="Medición");a[i].axvline(r.pred,color="#1d465b",label="Predicción")
        a[i].set(title=r.objetivo+" · "+r.cultivo,xlabel="Concentración (mM)",ylabel="Sorteos");a[i].legend(fontsize=6)
    for i,y in enumerate(PANEL[:7]):
        s=z.loc[z.objetivo.eq(y),"ancho"]
        if len(s):a[3].boxplot(s,positions=[i],widths=.5,flierprops=dict(markersize=2))
    a[3].set(xticks=range(7),xticklabels=PANEL[:7],ylabel="Ancho del 95 % (mM)",title="Anchura por metabolito")
    a[3].set_yscale("symlog",linthresh=.002)
    for ax,k,t in [(a[4],"error_abs","Error observado (mM)"),(a[6],"predictores","Predictores visibles¹"),(a[7],"faltante_train_pct","Ausencias en entrenamiento (%)")]:
        xy(ax,z[k],z.ancho);ax.set(xlabel=t,ylabel="Ancho del 95 % (mM)")
        ax.set_yscale("symlog",linthresh=.003)
        if k=="error_abs":ax.set_xscale("symlog",linthresh=.003)
    c=leer("calibracion_posterior").query("tipo=='individual'")
    for cv,h in c.groupby("prueba"):
        a[5].plot(h.nominal,h.observada,"o-",ms=4,label="Regla: "+cv)
    fijo=leer("calibracion_iterativo_fijo").query("objetivo=='todos'")
    a[5].plot(fijo.nominal,fijo.observada,"s--",ms=3,label="Iterativo fijo: perturbación")
    a[5].plot([.45,1],[.45,1],"k--",lw=.8);a[5].set(xlabel="Cobertura nominal",ylabel="Fracción observada",title="Calibración");a[5].legend(fontsize=6)
    a[6].set_xlabel("Metabolitos visibles del ajuste")
    salvar(fig,"09_posterior_y_calibracion")


def archivo_posterior(r):
    for path in sorted((PRIV/"posteriores"/r.fold).glob(r.objetivo+"__*.npz")):
        with np.load(path) as z:
            if r.evento in z["eventos"]:
                return np.load(path)
    raise ValueError(f"No se encontró la distribución posterior de {r.evento}")


def ejecutar():
    ish=leer_matriz("ishii");comps=dict(matrices(ish))
    original(ish)
    salvar(por_metabolito_original("fum",ish), "01_fumarato_original")
    for y in ["fum","cit","atp"]:
        salvar(por_metabolito(y,ish=ish,completadas=comps),"02_distribucion_"+y)
    estadisticos(ish,comps)
    relaciones(ish,comps)
    z=pd.read_csv(DET/"validacion.csv.gz",low_memory=False)
    z=z[z.fold.str.startswith("E") & z.tipo.eq("individual")]
    errores(z);decisiones(z);incertidumbre()


if __name__=="__main__":
    ejecutar()

