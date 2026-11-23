"""Contrastes por metabolito y figura de los dos experimentos nuevos."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

from .medir_bloques import medir


def contrastes_y_figura(pred, glob, met, por_grupo, destino):
    contrastes = []
    for alt in ("opciones_A", "reajuste"):
        contrastes.append(("fuentes", alt, alt, "A", "B", "MMC2_"+alt))
    contrastes.append(("curva", "50pct", "100pct", "A", "A", "tamano_50_a_100"))
    contrastes.extend([("ocho_filas", "base", "mas_8_Ishii", "A", "A", "sumar_ocho_Ishii"),
                       ("ocho_filas", "base", "mas_8_MMC2", "A", "B", "sumar_ocho_MMC2"),
                       ("ocho_filas", "mas_8_Ishii", "mas_8_MMC2", "A", "B", "igual_numero_de_filas")])
    for esc in ("A", "B"):
        for tipo in ("documental", "elegible"):
            contrastes.extend([("patron_"+tipo, "individual", "conjunto", esc, esc, "patron_"+tipo),
                               ("ajuste_"+tipo, "fijo", "reajuste", esc, esc, "reajuste_"+tipo)])
    resultados, detalle = [], []
    for exp, alt1, alt2, a, b, etiqueta in contrastes:
        def pareja(tabla, claves):
            t = tabla[(tabla.experimento.eq(exp)) & tabla.alcance.eq("comunes")]
            u = t[t.alternativa.eq(alt1) & t.escenario.eq(a)]
            v = t[t.alternativa.eq(alt2) & t.escenario.eq(b)]
            return u.merge(v, on=claves, suffixes=("_1","_2"), validate="one_to_one")
        clave = ["metodo", "metabolito"]
        g = pareja(por_grupo, clave+["grupo"])
        g["delta_MAE"] = g.MAE_mM_2-g.MAE_mM_1
        g = g[g.delta_MAE.notna()].copy()
        g["mejora"] = g.delta_MAE < -1e-12
        g["empeora"] = g.delta_MAE > 1e-12
        resumen = g.groupby(clave).agg(grupos_mejora=("mejora","sum"), grupos_empeora=("empeora","sum"),
            grupos_con_comparacion=("grupo","nunique"), mediana_delta_grupos=("delta_MAE","median"))
        m = pareja(met, clave)
        assert m.celdas_predichas_1.equals(m.celdas_predichas_2)
        m["delta_MAE_mM"] = m.MAE_mM_2 - m.MAE_mM_1
        m["delta_sesgo_mM"] = m.sesgo_mM_2 - m.sesgo_mM_1
        cols = clave+["celdas_predichas_1","eventos_predichos_1","grupos_evaluados_1","MAE_mM_1",
                      "MAE_mM_2","delta_MAE_mM","sesgo_mM_1","sesgo_mM_2","delta_sesgo_mM"]
        resultados.append(m[cols].join(resumen, on=clave).assign(contraste=etiqueta, escenario_1=a, escenario_2=b))
        detalle.append(g[clave+["grupo","celdas_predichas_1","MAE_mM_1","MAE_mM_2","delta_MAE",
                                "sesgo_mM_1","sesgo_mM_2"]].assign(contraste=etiqueta, escenario_1=a, escenario_2=b))
    cambios = pd.concat(resultados, ignore_index=True)
    cambios.to_csv(destino / "cambios_por_metabolito.csv", index=False)
    foco = cambios.contraste.eq("MMC2_opciones_A") & cambios.metabolito.isin(
        ["fum","cit","ade","atp","pep","s7p"])
    cambios.loc[foco].to_csv(destino / "objetivos_MMC2.csv", index=False, float_format="%.6g")
    pd.concat(detalle, ignore_index=True).to_csv(destino / "detalle/cambios_por_grupo.csv.gz", index=False,
                                                compression={"method":"gzip","mtime":0})
    curva = []
    for rep, v in pred[pred.experimento.eq("curva")].groupby("rep"):
        r, _, _ = medir(v)
        curva.append(r[r.alcance.eq("comunes")].assign(rep=rep))
    curva = pd.concat(curva, ignore_index=True)
    curva.to_csv(destino / "curva_repeticiones.csv", index=False)
    cols = ["experimento","alternativa","escenario","metodo","intentos","eventos_predichos",
            "celdas_predichas","grupos_evaluados","metabolitos_evaluados","MAE_mM_macro","MAE_mM_micro",
            "RMSE_mM_micro","sesgo_mM_micro","MAE_IQR_macro","MAE_R_macro"]
    glob.loc[glob.alcance.eq("comunes"),cols].to_csv(destino/"comparacion.csv", index=False,
                                                     float_format="%.7g")
    matplotlib.rcParams.update({"svg.fonttype":"none","svg.hashsalt":"IPre-bloques-1-a-4"})
    fig, (izq, der) = plt.subplots(1, 2, figsize=(10, 4.7), constrained_layout=True)
    escala = ["50pct", "75pct", "100pct"]
    colores = {"Mediana":"#62676b","KNN":"#007c79","Iterative":"#a24f19"}
    x = np.arange(3)
    for metodo, color in colores.items():
        sub = curva[curva.metodo.eq(metodo)].groupby("alternativa").MAE_mM_macro
        izq.plot(x, sub.median().reindex(escala), marker="o", color=color, label=metodo)
        izq.fill_between(x, sub.quantile(.25).reindex(escala), sub.quantile(.75).reindex(escala),
                         color=color, alpha=.12)
    izq.set(xticks=x, xticklabels=["50 %", "75 %", "100 %"], xlabel="Grupos Ishii de entrenamiento",
            ylabel="MAE macro (mM)", title="Añadir cultivos de Ishii")
    izq.legend(frameon=False)
    datos = glob[(glob.experimento=="patron_elegible") & (glob.alcance=="comunes")].set_index(
        ["escenario","metodo","alternativa"])
    ticks = []
    for i,(esc, metodo) in enumerate([(s,m) for s in ("A","B") for m in colores]):
        a,b = [datos.loc[(esc,metodo,k), "MAE_mM_macro"] for k in ["individual","conjunto"]]
        der.plot([a,b],[i,i],lw=2,color=colores[metodo])
        # El círculo abierto representa el ocultamiento individual en las mismas celdas.
        der.plot(a,i,"o",mfc="white",mec=colores[metodo],ms=8)
        der.plot(b,i,"x",color=colores[metodo],mew=2,ms=7)
        ticks.append(esc+" · "+metodo)
    der.set(yticks=range(6), yticklabels=ticks, xlabel="MAE macro (mM)",
            title="Ocultar uno (○) o varios (×)")
    der.invert_yaxis()
    for ax, axis in [(izq,izq.yaxis),(der,der.xaxis)]:
        ax.spines[["top","right"]].set_visible(False)
        axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.2f}".replace(".",",")))
        ax.grid(alpha=.19)
    fig.savefig(destino/"cultivos_y_patrones.png", dpi=160)
    fig.savefig(destino/"cultivos_y_patrones.svg", metadata={"Date":None})
    plt.close(fig)
