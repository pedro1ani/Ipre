"""Gráficos de resultados existentes"""

import argparse
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from imputadores import METHODS

RAIZ = Path(__file__).resolve().parent


def main():
    argumentos = argparse.ArgumentParser()
    argumentos.add_argument(
        "--results", type=Path, default=RAIZ / "resultados/ishii"
    )
    salida = argumentos.parse_args().results
    predicciones = pd.read_csv(salida / "validacion/predicciones.csv")
    carpeta_figuras = salida / "figuras"
    carpeta_figuras.mkdir(exist_ok=True)
    for prueba, tabla_prueba in predicciones.groupby("experiment"):
        metabolitos = list(tabla_prueba.metabolite_id.unique())
        figura, ejes = plt.subplots(
            len(metabolitos),
            3,
            figsize=(12, 3 * len(metabolitos)),
            squeeze=False,
        )
        for i, metabolito in enumerate(metabolitos):
            for j, metodo in enumerate(METHODS):
                eje = ejes[i, j]
                valores = tabla_prueba[
                    (tabla_prueba.metabolite_id == metabolito)
                    & (tabla_prueba.method == metodo)
                ].dropna(subset=["estimate_mM"])
                con_soporte = valores.support == "supported"
                for seleccionados, color, etiqueta in [
                    (con_soporte, "#147d92", "Con soporte"),
                    (~con_soporte, "#d78b35", "General o limitado"),
                ]:
                    grupo = valores[seleccionados]
                    eje.scatter(
                        grupo.truth_mM,
                        grupo.estimate_mM,
                        color=color,
                        s=23,
                        alpha=0.7,
                        label=etiqueta,
                    )
                maximo = (
                    max(
                        valores.truth_mM.max(),
                        valores.estimate_mM.max(),
                        0.001,
                    )
                    if len(valores)
                    else 1
                )
                eje.plot([0, maximo], [0, maximo], "--", color="gray", lw=1)
                eje.set_title(
                    f"{metodo}\n{metabolito}; n={len(valores)}", fontsize=10
                )
                eje.set_xlabel("Medición ocultada (mM)")
                eje.set_ylabel("Estimación (mM)")
                eje.grid(alpha=0.15)
        ejes[0, 0].legend(fontsize=8)
        figura.suptitle(
            f"{prueba}: valor conocido vs imputado\n"
            "Todos los candidatos numéricos; diagonal = predicción exacta"
        )
        figura.tight_layout(rect=(0, 0, 1, 0.95))
        figura.savefig(
            carpeta_figuras / f"{prueba}_real_vs_imputado.png", dpi=150
        )
        plt.close(figura)
    print(f"Gráficos: {carpeta_figuras}")


if __name__ == "__main__":
    main()
