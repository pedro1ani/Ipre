"""Prepara las tablas breves y la figura que lee el notebook de evaluación."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .diseno_5_6 import DESTINO
from .medir_bloques import medir
from .selector import EMPATE_MM, METODOS


def comparar_selector(destino):
    detalle = destino / "detalle"
    pred = pd.read_csv(detalle / "selector_predicciones.csv.gz")
    estricto = pred.copy()
    descartar = estricto.metodo.eq("Selector") & estricto.politica.ne("eleccion_interna")
    estricto.loc[descartar, "pred_mM"] = np.nan
    glob, _, _ = medir(estricto)
    glob.to_csv(detalle / "solo_con_eleccion.csv.gz", index=False)
    comparaciones = []
    for (panel, tipo), t in pred.groupby(["panel", "tipo"]):
        casos = ["grupo", "celda", "evento"]
        errores = t.assign(error=(t.pred_mM-t.verdad_mM).abs()).pivot(
            index=casos, columns=["escenario","metodo"], values="error").dropna()
        por_celda = errores.groupby(level=casos[:2]).mean()
        por_grupo = por_celda.groupby(level="grupo").mean()
        for esc in ("A","B"):
            for fijo in METODOS:
                dif = por_grupo[esc, "Selector"] - por_grupo[esc, fijo]
                comparaciones.append(dict(panel=panel, tipo=tipo, escenario=esc, fijo=fijo,
                    grupos=len(dif), mejor=int((dif < -EMPATE_MM).sum()), peor=int((dif > EMPATE_MM).sum()),
                    empate=int((dif.abs() <= EMPATE_MM).sum()), mediana_delta_MAE_mM=dif.median(),
                    q25_delta=dif.quantile(.25), q75_delta=dif.quantile(.75)))
    pd.DataFrame(comparaciones).to_csv(destino / "comparacion_grupos.csv", index=False)
    g = pd.read_csv(destino / "selector_global.csv").query("alcance == 'comunes'").copy()
    g[["panel","tipo"]] = g.experimento.str.split("_",n=1,expand=True)
    index = ["panel", "tipo", "escenario", "intentos", "celdas_predichas", "grupos_evaluados"]
    g.pivot(index=index, columns="metodo", values="MAE_mM_macro").reset_index().to_csv(
        destino / "selector_lectura.csv", index=False)
    met = pd.read_csv(detalle / "selector_por_metabolito.csv.gz")
    met = met[met.experimento.str.endswith("_individual") & met.alcance.eq("comunes") & met.celdas_predichas.gt(0)]
    met.pivot(index=["experimento", "escenario", "metabolito", "celdas_predichas", "grupos_evaluados"],
              columns="metodo", values="MAE_mM").to_csv(destino / "selector_por_metabolito.csv")


def figura(a, b, destino):
    plt.rcParams.update({"font.size":9, "svg.fonttype":"none", "svg.hashsalt":"IPre56",
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.6, 4.4), layout="constrained")
    a = a.query("tipo == 'individual' and panel == 'comun'").set_index("escenario")
    nombres = ["Mediana","KNN","Iterative","Selector"]
    for i, nombre in enumerate(nombres):
        ax.plot([a.at["A",nombre], a.at["B",nombre]], [i,i], color="#b7c2bf", lw=2, zorder=1)
    ax.scatter(a.loc["A",nombres], range(4), color="#175f69", s=55, label="A · Ishii", zorder=2)
    ax.scatter(a.loc["B",nombres], range(4), color="#ae4723", marker="D", s=39, label="B · Ishii + MMC2", zorder=2)
    ax.set(yticks=range(4), yticklabels=["Mediana","KNN","Iterativo","Selector"],
           xlabel="MAE macro (mM) · menor es mejor", title="1. Error por método · panel común")
    ax.invert_yaxis()
    ax.grid(axis="x",alpha=.25)
    ax.legend(loc="lower right",fontsize=8)
    bx.plot([45,100],[45,100],ls="--",color="0.55",label="Cobertura nominal")
    q=b.query("tipo == 'individual' and sorteos == 200")
    for (panel,esc), v in q.groupby(["panel","escenario"]):
        color = "#175f69" if esc=="A" else "#ae4723"
        estilo = "-" if panel=="comun" else ":"
        bx.plot(v.nivel,v.cobertura*100,color=color,ls=estilo,marker="o" if panel=="comun" else "s",
                markersize=4,label=f"{esc} · {'común' if panel=='comun' else 'ampliado'}")
    bx.set(xlim=(45,100),ylim=(40,100),xticks=[50,80,95],yticks=[40,50,60,70,80,90,100],
           xlabel="Nivel nominal (%)",ylabel="Mediciones dentro del intervalo (%)",
           title="2. ¿Cuántas mediciones cubre el intervalo?")
    bx.grid(alpha=.2)
    bx.legend(fontsize=7,loc="upper left")
    fig.suptitle("Validación externa · las mismas 416 mediciones de Ishii",weight="bold",fontsize=13)
    fig.savefig(destino / "selector_e_intervalos.svg",metadata={"Date":None})
    fig.savefig(destino / "selector_e_intervalos.png",dpi=160)
    plt.close(fig)


def ejecutar(destino=DESTINO):
    comparar_selector(destino)
    g = pd.read_csv(destino / "intervalos_global.csv").query("alcance == 'comunes' and serie == 'prefijo'")
    cols = ["panel","tipo","escenario","sorteos","nivel","intentos","eventos_con_intervalo",
            "aciertos_evento","celdas_con_intervalo","grupos_evaluados","cobertura","cobertura_macro",
            "ancho_mM","mediana_ancho_mM","sesgo_mM","bajo","sobre"]
    g[cols].to_csv(destino / "calibracion_lectura.csv", index=False)
    m = pd.read_csv(destino / "intervalos_por_metabolito_200.csv")
    m = m.query("tipo == 'individual' and nivel == 95 and celdas_con_intervalo > 0")
    m[["panel","escenario","metabolito","celdas_con_intervalo","grupos_evaluados","aciertos_evento",
       "cobertura","ancho_mM","sesgo_mM","bajo","sobre","SD_condicional_mM","SD_multiple_mM",
       "n_objetivo_train","n_predictores"]].to_csv(destino / "calibracion_por_metabolito.csv", index=False)
    figura(pd.read_csv(destino / "selector_lectura.csv"), g, destino)


if __name__ == "__main__":
    ejecutar()
