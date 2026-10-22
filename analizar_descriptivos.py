"""Descriptivos y gráficos exploratorios de valores observados; sin
imputación."""

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

from apoyo import describir

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

RAIZ = Path(__file__).resolve().parent
ETIQUETAS = [
    "Fumarato",
    "PEP",
    "s7p",
    "Succinato",
    "Piruvato",
    "α-cetoglutarato",
    "Citrato",
]


def main():
    archivo = argparse.ArgumentParser()
    archivo.add_argument("--output", type=Path, default=RAIZ / "descriptivos")
    opciones = archivo.parse_args()
    salida = opciones.output.resolve()
    salida.mkdir(parents=True, exist_ok=False)
    entradas = []
    tablas = {}
    estadisticas = []
    for estudio, etiqueta in [("ishii", "Ishii"), ("mmc2", "MMC2")]:
        for unidad in ["mM", "R"]:
            archivo = RAIZ / "normalizacion" / f"matriz_{estudio}_{unidad}.csv"
            entradas.append(archivo)
            tabla = pd.read_csv(archivo, index_col=0)
            tablas[estudio, unidad] = tabla
            for metabolito in tabla:
                valores = tabla[metabolito].dropna()
                resumen = describir(valores)
                estadisticas.append(
                    dict(
                        group=etiqueta,
                        unit=unidad,
                        metabolite=metabolito,
                        total=len(tabla),
                        n=len(valores),
                        missing=len(tabla) - len(valores),
                        coverage_pct=100 * len(valores) / len(tabla),
                        mean=resumen["mean_mM"],
                        median=resumen["median_mM"],
                        SD=resumen["SD_mM"],
                        Q1=valores.quantile(0.25),
                        Q3=valores.quantile(0.75),
                        IQR=resumen["IQR_mM"],
                        CV_pct=resumen["CV_pct"],
                        min=valores.min(),
                        max=valores.max(),
                    )
                )
    estadisticas = pd.DataFrame(estadisticas)
    orden_unidad = pd.Categorical(
        estadisticas.unit, categories=["mM", "R"], ordered=True
    )
    estadisticas = (
        estadisticas.assign(orden=orden_unidad)
        .sort_values("orden", kind="stable")
        .drop(columns="orden")
    )
    estadisticas.to_csv(
        salida / "estadisticos.csv", index=False, encoding="utf-8-sig"
    )
    ishii = tablas["ishii", "R"]
    mmc2 = tablas["mmc2", "R"]
    ratios = pd.concat([ishii, mmc2])
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.facecolor": "white",
        }
    )
    figuras = []
    # Distribución de concentraciones por estudio.
    figura, ejes = plt.subplots(2, 4, figsize=(14, 8))
    aleatorio = np.random.default_rng(20261015)
    for j, (metabolito, etiqueta) in enumerate(zip(ratios.columns, ETIQUETAS)):
        eje = ejes.flat[j]
        for posicion, (nombre, datos, color) in enumerate(
            [("Ishii", ishii, "#196b9a"), ("MMC2", mmc2, "#d47929")], 1
        ):
            valores = datos[metabolito].dropna()
            eje.boxplot(
                valores,
                positions=[posicion],
                widths=0.45,
                showfliers=False,
                patch_artist=True,
                boxprops={"facecolor": color, "alpha": 0.22},
                medianprops={"color": "black"},
            )
            eje.scatter(
                posicion + aleatorio.uniform(-0.13, 0.13, len(valores)),
                valores,
                s=25,
                color=color,
                alpha=0.75,
            )
        eje.axhline(1, color="#999999", ls="--", lw=1)
        eje.set_yscale("log")
        eje.set_xticks(
            [1, 2],
            [
                "Ishii\nn=" + str(ishii[metabolito].count()),
                "MMC2\nn=" + str(mmc2[metabolito].count()),
            ],
        )
        eje.set_title(etiqueta)
        eje.set_ylabel("Ratio respecto a glucosa")
        eje.grid(axis="y", alpha=0.15)
    ejes.flat[7].axis("off")
    ejes.flat[7].text(
        0,
        0.85,
        (
            "Solo valores observados.\nEje vertical lo"
            "garítmico.\nLínea punteada: referencia = "
            "1.\nPuntos: muestras/condiciones.\nCajas: "
            "Q1–Q3 y mediana.\nBigotes: hasta 1,5 IQR."
            "\nEscalas distintas entre metabolitos.\nMM"
            "C2 incluye la fila de referencia."
        ),
        va="top",
        linespacing=1.6,
    )
    figura.suptitle("Distribuciones antes de imputar", fontsize=18)
    figura.tight_layout(rect=(0, 0, 1, 0.95))
    figuras.append(("01_distribuciones", figura))
    # El logaritmo se usa solo en esta visualización.
    figura, eje = plt.subplots(figsize=(12, 11))
    valores_ratio = ratios.to_numpy()
    mapa = np.ma.masked_invalid(
        np.log2(np.where(valores_ratio > 0, valores_ratio, np.nan))
    )
    colores = plt.get_cmap("RdBu_r").copy()
    colores.set_bad("#d6d6d6")
    limite = max(abs(mapa.min()), abs(mapa.max()))
    imagen = eje.imshow(
        mapa, aspect="auto", cmap=colores, vmin=-limite, vmax=limite
    )
    eje.set_xticks(range(7), ETIQUETAS, rotation=25, ha="right")
    eje.set_yticks(range(len(ratios)), ratios.index, fontsize=8)
    eje.axhline(29.5, color="black", lw=2)
    figura.colorbar(
        imagen,
        ax=eje,
        label="log₂(ratio): 0 = referencia; +1 = doble; −1 = mitad",
    )
    eje.set_title(
        (
            "Patrones por muestra y metabolito\nGris ="
            " faltante; línea negra separa Ishii y MM"
            "C2"
        ),
        fontsize=16,
    )
    figura.tight_layout()
    figuras.append(("02_mapa_muestras", figura))
    figura, ejes = plt.subplots(1, 2, figsize=(13, 5))
    valores = np.arange(7)
    ancho = 0.35
    for k, (nombre, tabla, color) in enumerate(
        [("Ishii", ishii, "#196b9a"), ("MMC2", mmc2, "#d47929")]
    ):
        cobertura = tabla.notna().mean() * 100
        barras = ejes[0].bar(
            valores + (k - 0.5) * ancho,
            cobertura,
            ancho,
            label=nombre,
            color=color,
        )
        for barra, pares in zip(barras, tabla.count()):
            ejes[0].text(
                barra.get_x() + barra.get_width() / 2,
                barra.get_height() + 2,
                str(pares),
                ha="center",
                fontsize=8,
            )
    ejes[0].set_xticks(valores, ETIQUETAS, rotation=35, ha="right")
    ejes[0].set_ylim(0, 115)
    ejes[0].set_ylabel("Cobertura (%)")
    ejes[0].legend()
    ejes[0].set_title(
        "Cobertura por metabolito\nEtiquetas = número de valores medidos"
    )
    conteo = ishii.notna().sum(axis=1)
    frecuencias = conteo.value_counts().reindex(range(8), fill_value=0)
    ejes[1].bar(frecuencias.index, frecuencias, color="#196b9a")
    ejes[1].set_xticks(range(8))
    ejes[1].set_xlabel("Metabolitos observados de 7")
    ejes[1].set_ylabel("Número de muestras Ishii")
    ejes[1].set_title("Información disponible por muestra")
    figura.tight_layout()
    figuras.append(("03_cobertura", figura))
    # Correlaciones por pares medidos, sin completar huecos.
    pares = ishii.notna().astype(int).T.dot(ishii.notna().astype(int))
    correlaciones = ishii.corr(method="spearman", min_periods=8)
    correlaciones.to_csv(
        salida / "correlacion_spearman_ishii.csv", encoding="utf-8-sig"
    )
    pares.to_csv(salida / "n_pares_ishii.csv", encoding="utf-8-sig")
    figura, eje = plt.subplots(figsize=(10, 8))
    imagen = eje.imshow(
        np.ma.masked_invalid(correlaciones), vmin=-1, vmax=1, cmap="RdBu_r"
    )
    eje.set_xticks(range(7), ETIQUETAS, rotation=35, ha="right")
    eje.set_yticks(range(7), ETIQUETAS)
    for i in range(7):
        for j in range(7):
            valor = correlaciones.iloc[i, j]
            eje.text(
                j,
                i,
                ("—" if pd.isna(valor) else f"{valor:.2f}")
                + f"\nn={pares.iloc[i, j]}",
                ha="center",
                va="center",
                color=(
                    "white"
                    if pd.notna(valor) and abs(valor) > 0.65
                    else "black"
                ),
                fontsize=10,
            )
    figura.colorbar(imagen, ax=eje, label="ρ de Spearman")
    eje.set_title(
        (
            "Asociaciones exploratorias dentro de Ish"
            "ii\nPares observados; mínimo 8 pares; sin"
            " imputación"
        ),
        fontsize=15,
    )
    figura.tight_layout()
    figuras.append(("04_correlaciones_ishii", figura))
    with PdfPages(salida / "graficos_exploratorios.pdf") as pdf:
        for nombre, figura in figuras:
            figura.savefig(salida / (nombre + ".png"), dpi=180)
            pdf.savefig(figura)
            plt.close(figura)
    conteo.rename("observed_metabolites").to_csv(
        salida / "cobertura_muestras_ishii.csv", encoding="utf-8-sig"
    )
    manifest = {
        "stage": "exploratory_observed_only",
        "input_sha256": {
            str(archivo.relative_to(RAIZ)): hashlib.sha256(
                archivo.read_bytes()
            ).hexdigest()
            for archivo in entradas
        },
        "script_sha256": hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        "versions": {
            biblioteca: importlib.metadata.version(biblioteca)
            for biblioteca in ["numpy", "pandas", "matplotlib", "scipy"]
        },
        "SD_ddof": 1,
        "quantiles": "linear",
        "correlation": "Spearman pairwise complete, min 8",
        "jitter_seed": 20261015,
        "outputs_sha256": {
            archivo.name: hashlib.sha256(archivo.read_bytes()).hexdigest()
            for archivo in salida.iterdir()
        },
    }
    (salida / "manifest_descriptivos.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(
        estadisticas.query("unit=='R' and group!='Combinado'")[
            ["group", "metabolite", "n", "median", "CV_pct", "min", "max"]
        ].to_string(index=False)
    )
    print("\nCORRELACIONES\n", correlaciones.round(3).to_string())
    print("\nMUESTRAS CON <=2 VALORES\n", conteo[conteo <= 2].to_string())


if __name__ == "__main__":
    main()
