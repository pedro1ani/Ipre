"""Tablas, histogramas y diagnóstico de intervalos del estudio ampliado."""

import argparse
import itertools

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from estudio_ampliado import METHODS, ROOT, SEED, summarize

COLORS = {"Mediana": "#6e7781", "KNN": "#007c91", "Iterative": "#a84c89"}


def ecdf_distance(x, y):
    grid = np.sort(np.r_[x, y])
    return np.max(
        np.abs(
            np.searchsorted(np.sort(x), grid, side="right") / len(x)
            - np.searchsorted(np.sort(y), grid, side="right") / len(y)
        )
    )


def paired_distribution_test(table, groups):
    """Intercambia medido/estimado por grupo; inferencia exploratoria."""
    x = table.truth_mM.to_numpy()
    y = table.estimate_mM.to_numpy()
    labels = groups.loc[table["sample"]].to_numpy()
    unique = np.unique(labels)
    observed = ecdf_distance(x, y)
    if len(unique) <= 12:
        switches = np.array(
            list(itertools.product([False, True], repeat=len(unique)))
        )
        exact = True
    else:
        switches = np.random.default_rng(SEED).integers(
            0, 2, (1999, len(unique))
        )
        exact = False
    indices = np.searchsorted(unique, labels)
    statistics = []
    for switch in switches:
        swap = switch[indices].astype(bool)
        statistics.append(
            ecdf_distance(np.where(swap, y, x), np.where(swap, x, y))
        )
    exceed = np.sum(np.asarray(statistics) >= observed - 1e-12)
    pvalue = (
        exceed / len(statistics)
        if exact
        else (exceed + 1) / (len(statistics) + 1)
    )
    return (
        observed,
        pvalue,
        len(unique),
        "exact" if exact else "1999_permutations",
    )


def analyze(dataset):
    out = ROOT / "resultados_ampliados" / dataset
    original = pd.read_csv(out / f"matriz_{dataset}_mM.csv", index_col=0)
    pred = pd.read_csv(out / "predicciones.csv")
    candidates = pd.read_csv(out / "candidatos.csv")
    macro = summarize(original, pred, candidates, out)
    groups = (
        pd.read_csv(ROOT / "normalizacion/muestras_y_exclusiones.csv")
        .set_index("sample_id")
        .perturbation_group
    )
    tests = []
    for (method, metabolite), g in pred.groupby(["method", "metabolite"]):
        g = g.dropna(subset=["estimate_mM"])
        stat, pvalue, n, mode = paired_distribution_test(g, groups)
        tests.append(
            dict(
                method=method,
                metabolite=metabolite,
                ECDF_distance=stat,
                p_value=pvalue,
                n_groups=n,
                mode=mode,
            )
        )
    tests = pd.DataFrame(tests)
    order = np.argsort(tests.p_value.to_numpy())
    tests["p_holm"] = np.nan
    tests.loc[order, "p_holm"] = np.minimum(
        1,
        np.maximum.accumulate(
            tests.p_value.to_numpy()[order]
            * (len(tests) - np.arange(len(tests)))
        ),
    )
    tests.to_csv(out / "distribuciones_validacion_test.csv", index=False)
    figs = out / "figuras"
    figs.mkdir(exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 120,
        }
    )
    for method in METHODS:
        completed = pd.read_csv(out / f"imputada_{method}.csv", index_col=0)
        fig, axes = plt.subplots(
            5, 2, figsize=(11, 15), constrained_layout=True
        )
        for ax, metabolite in zip(axes.flat, original):
            before = original[metabolite].dropna()
            after = completed[metabolite].dropna()
            bins = np.histogram_bin_edges(np.r_[before, after], bins=7)
            ax.hist(
                before,
                bins=bins,
                density=True,
                histtype="step",
                linewidth=2,
                color="#202d3b",
                label=f"Medido (n={len(before)})",
            )
            ax.hist(
                after,
                bins=bins,
                density=True,
                alpha=0.28,
                color=COLORS[method],
                label=f"Completado (n={len(after)})",
            )
            ax.set(
                title=metabolite,
                xlabel="Concentración (mM)",
                ylabel="Densidad",
            )
            ax.legend(fontsize=8)
        fig.suptitle(
            f"{dataset.upper()} · {method} · Antes y después\n"
            "Histogramas descriptivos; comparten las mediciones originales",
            fontsize=15,
        )
        fig.savefig(figs / f"distribuciones_{method}.png")
        plt.close(fig)
    metrics = pd.read_csv(out / "metricas.csv")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    pivot = metrics.pivot(index="metabolite", columns="method", values="MAE_R")
    pivot[METHODS].plot.bar(ax=axes[0], color=[COLORS[m] for m in METHODS])
    axes[0].set(
        title="Error por metabolito · evaluación externa",
        ylabel="MAE en ratios",
        xlabel="",
    )
    macro.loc[METHODS, "MAE_R"].plot.bar(
        ax=axes[1], color=[COLORS[m] for m in METHODS]
    )
    axes[1].set(
        title="Promedio con igual peso por metabolito",
        ylabel="Macro MAE en ratios",
        xlabel="",
    )
    fig.savefig(figs / "comparacion_metodos.png")
    plt.close(fig)
    for metabolite in ["fum[c]", "atp[c]"]:
        g = pred[
            (pred.method == "Iterative") & (pred.metabolite == metabolite)
        ].sort_values("truth_mM")
        fig, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
        positions = np.arange(len(g))
        ax.vlines(
            positions,
            g.posterior_low_mM,
            g.posterior_high_mM,
            color=COLORS["Iterative"],
            alpha=0.65,
            label="Cuantiles posteriores 2,5–97,5 %",
        )
        ax.scatter(
            positions,
            g.truth_mM,
            color="#202d3b",
            label="Medición ocultada",
            zorder=3,
        )
        ax.scatter(
            positions,
            g.posterior_mean_mM,
            marker="x",
            color=COLORS["Iterative"],
            label="Media de 30 imputaciones",
            zorder=4,
        )
        ax.set_xticks(positions, g["sample"], rotation=90)
        ax.set(
            title=(
                f"{dataset.upper()} · {metabolite} · "
                "Incertidumbre contrastada con mediciones"
            ),
            ylabel="Concentración (mM)",
            xlabel="Muestras ordenadas por medición",
        )
        ax.legend(fontsize=8)
        fig.savefig(figs / f'intervalos_{metabolite.split("[")[0]}.png')
        plt.close(fig)
    coverage = pd.read_csv(out / "intervalos_evaluacion.csv")
    winner = macro.MAE_R.idxmin()
    iterative = pred[pred.method == "Iterative"]
    covered = (iterative.truth_mM >= iterative.posterior_low_mM) & (
        iterative.truth_mM <= iterative.posterior_high_mM
    )
    convergence = iterative.converged.mean()
    diagnostic = pred.groupby("method").agg(
        predictions=("estimate_mM", "count"),
        converged_fraction=("converged", "mean"),
        minimum_target_count=("train_target_n", "min"),
    )
    diagnostic.to_csv(out / "diagnostico_ajustes.csv")
    report = [
        f"# Resultados ampliados: {dataset.upper()}",
        "",
        f"El menor macro MAE_R externo corresponde a **{winner}**.",
        "",
        "```text\n" + macro.to_string() + "\n```",
        "",
        "Cada metabolito tiene el mismo peso. El ranking es exploratorio; "
        "no demuestra "
        "superioridad universal ni valida los huecos reales. La selección "
        "de parámetros "
        "se hizo dentro del entrenamiento de cada partición.",
        "",
        f"Iterative: {100*convergence:.1f} % de las predicciones externas "
        "proceden de ajustes que cumplieron el criterio de parada. "
        "El ranking incluye todas las predicciones finitas, también las "
        "no convergentes; este diagnóstico debe acompañar su interpretación.",
        "",
        "## Incertidumbre",
        "",
        f"Los cuantiles posteriores contienen {int(covered.sum())} de "
        f"{len(covered)} mediciones ocultadas ({100*covered.mean():.1f} %). "
        "La cobertura observada no alcanza el 95 % nominal de los cuantiles.",
        "",
        "```text\n" + coverage.to_string(index=False) + "\n```",
        "",
        "Los intervalos residuales son rangos empíricos de errores "
        "internos, no intervalos "
        "con cobertura garantizada. Con menos de cinco muestras de "
        "calibración se omiten. "
        "Los intervalos posteriores de Iterative provienen de 30 "
        "reajustes con semillas "
        "distintas. Sus cuantiles son aproximados; se reporta su "
        "cobertura externa real. "
        "La SD es de imputación bajo el modelo, no la SD experimental. No "
        "incluye la "
        "incertidumbre del volumen celular ni de la referencia Gerosa.",
        "",
        "## Distribuciones y pruebas",
        "",
        "Los histogramas usan los mismos intervalos de clase "
        "antes/después. Se reportan "
        "media, mediana, SD, IQR, CV y distancia de Wasserstein. No se "
        "calcula t-test, "
        "Mann–Whitney ni KS independiente sobre antes/después porque las "
        "muestras se "
        "solapan. Para valores ocultados se reporta distancia entre ECDF "
        "y permutación "
        "pareada por grupo, junto con prueba de signos de sesgo y "
        "corrección de Holm. "
        "Estas pruebas son exploratorias: requieren intercambiabilidad o "
        "simetría de "
        "errores por grupo y los modelos de validación cruzada comparten "
        "entrenamiento. "
        "Un p alto no demuestra equivalencia de distribuciones ni "
        "imputaciones correctas.",
        "",
        "## Límites del estudio",
        "",
        "- Diez metabolitos, 30 muestras Ishii u ocho condiciones MMC2, "
        "por separado.",
        "- Glucosa MMC2 se excluye de evaluación: es la referencia conocida.",
        "- La evaluación nueva oculta una celda por consulta; las "
        "ausencias reales se "
        "imputan conjuntamente. No sustituye las pruebas de patrones "
        "reales y series "
        "analíticas conservadas en el flujo anterior.",
        "- Cada candidato real se obtiene sin su bloque externo de "
        "grupos. Se usa la "
        "misma selección interna; no se reajusta usando la muestra "
        "consultada.",
        "- ATP, ADP y AMP están completos en la cohorte Ishii: amplían la "
        "evaluación "
        "y los predictores, pero no generan huecos reales artificiales.",
        "- La reducción de dispersión no se interpreta como una mejora.",
        "- Los resultados y parámetros nuevos no tienen que coincidir con "
        "los antiguos.",
        "",
    ]
    (out / "CONCLUSIONES.md").write_text("\n".join(report), encoding="utf8")
    print(dataset, macro.to_string(), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=["ishii", "mmc2"], required=True)
    analyze(parser.parse_args().dataset)
