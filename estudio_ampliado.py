"""Diez metabolitos: selección interna, evaluación externa e incertidumbre.

python estudio_ampliado.py --dataset ishii
python estudio_ampliado.py --dataset mmc2
"""

import argparse
import hashlib
import json
import sqlite3
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, wasserstein_distance
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer, KNNImputer, SimpleImputer
from sklearn.linear_model import BayesianRidge
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent
SEED = 20261022
DRAWS = 30
METHODS = ["Mediana", "KNN", "Iterative"]
# Rejilla pequeña y predefinida; nunca se amplía mirando el error externo.
CONFIGS = {
    "Mediana": [{"id": "mediana"}],
    "KNN": [
        {"id": "k3_uniform", "k": 3, "weights": "uniform"},
        {"id": "k5_uniform", "k": 5, "weights": "uniform"},
        {"id": "k3_distance", "k": 3, "weights": "distance"},
    ],
    "Iterative": [
        {
            "id": "base",
            "rounds": 50,
            "tol": 0.001,
            "reg_iter": 300,
            "reg_tol": 0.001,
            "prior": 1e-6,
        },
        {
            "id": "convergencia",
            "rounds": 100,
            "tol": 0.0001,
            "reg_iter": 1000,
            "reg_tol": 0.0001,
            "prior": 1e-6,
        },
        {
            "id": "prior_0.01",
            "rounds": 50,
            "tol": 0.001,
            "reg_iter": 300,
            "reg_tol": 0.001,
            "prior": 0.01,
        },
    ],
}
ADDED = {
    "amp[c]": ("ISHII:metabolites:538", "AMP"),
    "adp[c]": ("ISHII:metabolites:554", "ADP"),
    "atp[c]": ("ISHII:metabolites:566", "ATP"),
}


def save_json(path, value):
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False), encoding="utf8"
    )


def prepare(out):
    """Añade equivalencias explícitas sin modificar la base original."""
    with sqlite3.connect(ROOT / "data/concentraciones_mM.db") as con:
        records = pd.read_sql("SELECT * FROM metabolite_measurements_mM", con)
    config = json.loads((ROOT / "config_normalizacion.json").read_text())
    ids = dict(config["ishii_analyte_ids"])
    ids.update({m: item[0] for m, item in ADDED.items()})
    selected = []
    audit = []
    references = []
    for metabolite, analyte in ids.items():
        ish = records[
            (records.dataset_id == "ISHII2007")
            & (records.analyte_id == analyte)
        ].copy()
        ref = records[
            (records.dataset_id == "MMC2")
            & (records.condition_id == "MMC2_GLUCOSE")
            & (records.metabolite_id == metabolite)
        ]
        if len(ref) != 1 or ish.empty:
            raise ValueError(f"Identidad o referencia ausente: {metabolite}")
        r = ref.iloc[0]
        assert r.value_mM > 0 and r.absence_code == "reported"
        assert r.uncertainty_type == "SD" and pd.notna(r.uncertainty_mM)
        assert (
            r.source_uncertainty_cell and not r.composite and not r.is_imputed
        )
        assert (ish.compartment == "intracellular").all()
        assert (ish.unit_original == "mM").all()
        assert not ish.composite.any() and not ish.is_imputed.any()
        if metabolite in ADDED:
            assert set(ish.original_label) == {ADDED[metabolite][1]}
        audit.append(
            dict(
                metabolite=metabolite,
                ishii_analyte=analyte,
                source_label=ish.original_label.iloc[0],
                mapping=(
                    "explicit_label_review"
                    if metabolite in ADDED
                    else "previous_audited_mapping"
                ),
                gerosa_label=r.original_label,
                reference_mM=r.value_mM,
                reference_SD_mM=r.uncertainty_mM,
                reference_uid=r.measurement_uid,
            )
        )
        ish["metabolite_id"] = metabolite
        ish["sample_id"] = ish.sample_key
        selected.append(ish)
        mmc = records[
            (records.dataset_id == "MMC2")
            & (records.metabolite_id == metabolite)
            & records.condition_id.isin(config["mmc2_conditions"])
        ].copy()
        assert (mmc.compartment == "intracellular").all()
        assert not mmc.composite.any() and not mmc.is_imputed.any()
        mmc["sample_id"] = mmc.condition_id
        selected.append(mmc)
        references.append(r)
    pd.DataFrame(audit).to_csv(out / "homologacion.csv", index=False)
    trace = pd.concat(selected, ignore_index=True)
    meta = pd.read_csv(ROOT / "normalizacion/muestras_y_exclusiones.csv")
    cohort = meta[meta.included].sample_id
    trace = trace[trace.sample_id.isin(cohort)]
    assert not trace.duplicated(
        ["dataset_id", "sample_id", "metabolite_id"]
    ).any()
    trace.to_csv(out / "trazabilidad.csv", index=False)
    pd.DataFrame(references).to_csv(out / "referencias.csv", index=False)
    refs = pd.Series({r.metabolite_id: r.value_mM for r in references})
    matrices = {}
    for name, dataset in [("ishii", "ISHII2007"), ("mmc2", "MMC2")]:
        table = trace[trace.dataset_id == dataset].pivot(
            index="sample_id", columns="metabolite_id", values="value_mM"
        )
        table = table.reindex(columns=list(ids)).sort_index()
        assert table.shape == ((30, 10) if name == "ishii" else (8, 10))
        table.to_csv(out / f"matriz_{name}_mM.csv")
        ratios = table / refs
        ratios.to_csv(out / f"matriz_{name}_R.csv")
        groups = meta.set_index("sample_id").loc[
            table.index, "perturbation_group"
        ]
        matrices[name] = (ratios, refs[table.columns], groups)
    return matrices


def cases(matrix, real=False):
    """Una consulta por celda observada; huecos reales juntos por muestra."""
    rows, info = [], []
    for sample, row in matrix.iterrows():
        if sample == "MMC2_GLUCOSE" and not real:
            continue
        if real:
            missing = row.index[row.isna()]
            if len(missing):
                rows.append(row.to_numpy().copy())
                for metabolite in missing:
                    info.append(
                        dict(
                            query=len(rows) - 1,
                            sample=sample,
                            metabolite=metabolite,
                            truth_R=np.nan,
                        )
                    )
        else:
            for metabolite in row.index[row.notna()]:
                query = row.copy()
                query[metabolite] = np.nan
                rows.append(query.to_numpy())
                info.append(
                    dict(
                        query=len(rows) - 1,
                        sample=sample,
                        metabolite=metabolite,
                        truth_R=float(row[metabolite]),
                    )
                )
    return pd.DataFrame(rows, columns=matrix.columns), pd.DataFrame(info)


def predict(train, query, method, config, posterior=False, seed=SEED):
    """Aprende todas las transformaciones solo con train."""
    active = train.columns[train.notna().any()]
    output = query.copy()
    log = dict(converged=True, rounds=0, warnings="")
    if len(active) == 0 or len(query) == 0:
        return output, log
    x, q = train[active], query[active]
    scaler = None
    if method == "Mediana":
        model = SimpleImputer(strategy="median")
    else:
        scaler = StandardScaler().fit(x)
        x, q = scaler.transform(x), scaler.transform(q)
        if method == "KNN":
            model = KNNImputer(
                n_neighbors=config["k"], weights=config["weights"]
            )
        else:
            prior = config["prior"]
            regressor = BayesianRidge(
                max_iter=config["reg_iter"],
                tol=config["reg_tol"],
                alpha_1=prior,
                alpha_2=prior,
                lambda_1=prior,
                lambda_2=prior,
            )
            model = IterativeImputer(
                estimator=regressor,
                max_iter=config["rounds"],
                tol=config["tol"],
                initial_strategy="median",
                skip_complete=False,
                sample_posterior=posterior,
                random_state=seed,
                min_value=-scaler.mean_ / scaler.scale_,
            )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model.fit(x)
        y = model.transform(q)
    if scaler is not None:
        y = scaler.inverse_transform(y)
    log["warnings"] = "; ".join(str(w.message) for w in caught)
    log["converged"] = not any(
        "Early stopping" in str(w.message) for w in caught
    )
    log["rounds"] = int(getattr(model, "n_iter_", 0))
    for j, metabolite in enumerate(active):
        missing = query[metabolite].isna()
        output.loc[missing, metabolite] = np.maximum(0, y[missing, j])
    return output, log


def extract(prediction, info, method, config, refs, train, log):
    rows = []
    for item in info.itertuples():
        value = prediction.loc[item.query, item.metabolite]
        rows.append(
            dict(
                sample=item.sample,
                metabolite=item.metabolite,
                method=method,
                config=config["id"],
                truth_R=item.truth_R,
                estimate_R=value,
                truth_mM=item.truth_R * refs[item.metabolite],
                estimate_mM=value * refs[item.metabolite],
                train_target_n=int(train[item.metabolite].notna().sum()),
                converged=log["converged"],
            )
        )
    return pd.DataFrame(rows)


def inner_selection(train, groups, refs, fold):
    """Valida por grupos dentro de train; no recibe las verdades externas."""
    results, splits = [], []
    splitter = GroupKFold(n_splits=min(3, groups.nunique()))
    for inner, (a, b) in enumerate(splitter.split(train, groups=groups)):
        fit, validation = train.iloc[a], train.iloc[b]
        assert not set(groups.iloc[a]) & set(groups.iloc[b])
        splits.append(
            dict(
                outer=fold,
                inner=inner,
                train=list(fit.index),
                validation=list(validation.index),
            )
        )
        query, info = cases(validation)
        for method in METHODS:
            for config in CONFIGS[method]:
                pred, log = predict(fit, query, method, config)
                part = extract(pred, info, method, config, refs, fit, log)
                part["inner_fold"] = inner
                results.append(part)
    all_predictions = pd.concat(results, ignore_index=True)
    all_predictions["absolute_error_R"] = abs(
        all_predictions.estimate_R - all_predictions.truth_R
    )
    # Intersección de casos: una opción no puede ganar absteniéndose más.
    all_predictions["case"] = (
        all_predictions["sample"] + "|" + all_predictions.metabolite
    )
    sets = [
        set(g.loc[g.estimate_R.notna(), "case"])
        for _, g in all_predictions.groupby(["method", "config"])
    ]
    common = set.intersection(*sets)
    scored = all_predictions[all_predictions.case.isin(common)]
    macro = scored.groupby(
        ["method", "config", "metabolite"]
    ).absolute_error_R.mean()
    scores = (
        macro.groupby(["method", "config"])
        .mean()
        .rename("macro_MAE_R")
        .reset_index()
    )
    scores["outer_fold"] = fold
    chosen = {}
    for method in METHODS:
        ordered = scores[scores.method == method].sort_values(
            ["macro_MAE_R", "config"]
        )
        if ordered.empty:
            raise ValueError(
                f"Sin validación interna comparable para {method}"
            )
        winner = ordered.iloc[0]["config"]
        chosen[method] = next(c for c in CONFIGS[method] if c["id"] == winner)
    all_predictions["outer_fold"] = fold
    return chosen, all_predictions, scores, splits


def residual_intervals(rows, inner, refs):
    """Rango empírico de errores internos; no promete cobertura nominal."""
    rows = rows.copy()
    lows, highs, counts = [], [], []
    for row in rows.itertuples():
        pool = inner[
            (inner.method == row.method)
            & (inner.config == row.config)
            & (inner.metabolite == row.metabolite)
        ].dropna(subset=["estimate_R"])
        # Promediar por muestra evita contar varias veces la misma unidad.
        errors = (
            (pool.truth_R - pool.estimate_R).groupby(pool["sample"]).mean()
        )
        counts.append(len(errors))
        if len(errors) < 5 or not np.isfinite(row.estimate_R):
            lows.append(np.nan)
            highs.append(np.nan)
            continue
        low, high = np.quantile(errors, [0.025, 0.975])
        lows.append(max(0, row.estimate_R + low) * refs[row.metabolite])
        highs.append(max(0, row.estimate_R + high) * refs[row.metabolite])
    rows["residual_low_mM"], rows["residual_high_mM"] = lows, highs
    rows["residual_n"] = counts
    return rows


def posterior_intervals(train, query, info, config, refs, seed):
    """Reajusta con semillas distintas para propagar también huecos de train."""
    draws = []
    for draw in range(DRAWS):
        pred, _ = predict(train, query, "Iterative", config, True, seed + draw)
        draws.append(
            [
                pred.loc[r.query, r.metabolite] * refs[r.metabolite]
                for r in info.itertuples()
            ]
        )
    values = np.asarray(draws)
    return values


def run(dataset, out, matrices):
    matrix, refs, groups = matrices[dataset]
    final, real, tuning, logs, inner_tables, draw_tables = (
        [],
        [],
        [],
        [],
        [],
        [],
    )
    outer = GroupKFold(n_splits=min(5, groups.nunique()))
    for fold, (a, b) in enumerate(outer.split(matrix, groups=groups)):
        train, test = matrix.iloc[a], matrix.iloc[b]
        assert not set(groups.iloc[a]) & set(groups.iloc[b])
        chosen, inner, scores, inner_splits = inner_selection(
            train, groups.iloc[a], refs, fold
        )
        tuning.append(scores)
        inner_tables.append(inner)
        logs.append(
            dict(
                fold=fold,
                train=list(train.index),
                test=list(test.index),
                chosen=chosen,
                inner_splits=inner_splits,
            )
        )
        print(f"{dataset}: fold {fold+1}, selección {chosen}", flush=True)
        for is_real in [False, True]:
            query, info = cases(test, real=is_real)
            if info.empty:
                continue
            for method in METHODS:
                config = chosen[method]
                pred, fit_log = predict(train, query, method, config)
                part = extract(
                    pred, info, method, config, refs, train, fit_log
                )
                part["fold"] = fold
                part = residual_intervals(part, inner, refs)
                if method == "Iterative":
                    draws = posterior_intervals(
                        train, query, info, config, refs, SEED + 1000 * fold
                    )
                    part["posterior_mean_mM"] = draws.mean(axis=0)
                    part["posterior_SD_mM"] = draws.std(axis=0, ddof=1)
                    part["posterior_low_mM"] = np.quantile(
                        draws, 0.025, axis=0
                    )
                    part["posterior_high_mM"] = np.quantile(
                        draws, 0.975, axis=0
                    )
                    for draw, values in enumerate(draws):
                        draw_table = info[["sample", "metabolite"]].copy()
                        draw_table["draw"] = draw
                        draw_table["value_mM"] = values
                        draw_table["fold"] = fold
                        draw_table["real_missing"] = is_real
                        draw_tables.append(draw_table)
                (real if is_real else final).append(part)
    predictions = pd.concat(final, ignore_index=True)
    candidates = pd.concat(real, ignore_index=True) if real else pd.DataFrame()
    predictions["error_mM"] = predictions.estimate_mM - predictions.truth_mM
    predictions["absolute_error_mM"] = predictions.error_mM.abs()
    predictions["absolute_error_R"] = (
        predictions.estimate_R - predictions.truth_R
    ).abs()
    for table, name in [
        (predictions, "predicciones"),
        (candidates, "candidatos"),
        (pd.concat(tuning), "hiperparametros"),
        (pd.concat(inner_tables), "validacion_interna"),
        (pd.concat(draw_tables), "muestras_posteriores"),
    ]:
        table.to_csv(out / f"{name}.csv", index=False)
    save_json(out / "particiones.json", logs)
    return matrix * refs, predictions, candidates


def summarize(original, predictions, candidates, out):
    p = predictions.copy()
    p["case"] = p["sample"] + "|" + p.metabolite
    sets = [
        set(g.loc[g.estimate_mM.notna(), "case"])
        for _, g in p.groupby("method")
    ]
    common = set.intersection(*sets)
    p = p[p.case.isin(common)]
    metrics = p.groupby(["method", "metabolite"]).agg(
        n=("estimate_mM", "count"),
        MAE_mM=("absolute_error_mM", "mean"),
        MAE_R=("absolute_error_R", "mean"),
        bias_mM=("error_mM", "mean"),
    )
    metrics.to_csv(out / "metricas.csv")
    macro = metrics.groupby("method")[["MAE_R", "MAE_mM"]].mean()
    macro.to_csv(out / "comparacion_global.csv")
    coverage, distributions, signed = [], [], []
    for (method, metabolite), g in p.groupby(["method", "metabolite"]):
        for interval in ["residual", "posterior"]:
            lo, hi = interval + "_low_mM", interval + "_high_mM"
            if lo not in g or (
                interval == "posterior" and method != "Iterative"
            ):
                continue
            valid = g.dropna(subset=[lo, hi])
            coverage.append(
                dict(
                    method=method,
                    metabolite=metabolite,
                    interval=interval,
                    n=len(valid),
                    coverage=(
                        (valid.truth_mM >= valid[lo])
                        & (valid.truth_mM <= valid[hi])
                    ).mean(),
                    mean_width_mM=(valid[hi] - valid[lo]).mean(),
                )
            )
        # Prueba de signos por grupo; descriptiva/exploratoria por CV compartida.
        diff = g.groupby("sample").error_mM.mean()
        meta = pd.read_csv(ROOT / "normalizacion/muestras_y_exclusiones.csv")
        group_map = meta.set_index("sample_id").perturbation_group
        diff.index = group_map.loc[diff.index].to_numpy()
        diff = diff.groupby(level=0).mean()
        nonzero = diff[np.abs(diff) > 1e-12]
        pv = (
            binomtest(int((nonzero > 0).sum()), len(nonzero), 0.5).pvalue
            if len(nonzero)
            else 1.0
        )
        signed.append(
            dict(
                method=method,
                metabolite=metabolite,
                test="group_sign_exploratory",
                n_groups=len(diff),
                p_value=pv,
                median_error_mM=diff.median(),
            )
        )
    pd.DataFrame(coverage).to_csv(
        out / "intervalos_evaluacion.csv", index=False
    )
    tests = pd.DataFrame(signed)
    # Holm en toda la familia de pruebas, no por metabolito a posteriori.
    order = np.argsort(tests.p_value.to_numpy())
    adjusted = np.minimum(
        1,
        np.maximum.accumulate(
            tests.p_value.to_numpy()[order]
            * (len(tests) - np.arange(len(tests)))
        ),
    )
    tests["p_holm"] = np.nan
    tests.loc[order, "p_holm"] = adjusted
    tests.to_csv(out / "pruebas_exploratorias.csv", index=False)
    for method in METHODS:
        completed = original.copy()
        if not candidates.empty:
            for row in candidates[candidates.method == method].itertuples():
                assert pd.isna(original.loc[row.sample, row.metabolite])
                completed.loc[row.sample, row.metabolite] = row.estimate_mM
        assert np.array_equal(
            completed.to_numpy()[original.notna()],
            original.to_numpy()[original.notna()],
        )
        completed.to_csv(out / f"imputada_{method}.csv")
        for metabolite in original:
            before = original[metabolite].dropna()
            after = completed[metabolite].dropna()
            for phase, values in [("antes", before), ("despues", after)]:
                distributions.append(
                    dict(
                        method=method,
                        metabolite=metabolite,
                        phase=phase,
                        n=len(values),
                        mean=values.mean(),
                        median=values.median(),
                        SD=values.std(),
                        IQR=values.quantile(0.75) - values.quantile(0.25),
                        CV=(
                            values.std() / values.mean()
                            if values.mean()
                            else np.nan
                        ),
                        wasserstein_from_before=wasserstein_distance(
                            before, values
                        ),
                    )
                )
    pd.DataFrame(distributions).to_csv(out / "distribuciones.csv", index=False)
    return macro


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=["ishii", "mmc2"], required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = args.output or ROOT / "resultados_ampliados" / args.dataset
    out.mkdir(parents=True, exist_ok=False)
    matrices = prepare(out)
    original, predictions, candidates = run(args.dataset, out, matrices)
    macro = summarize(original, predictions, candidates, out)
    import importlib.metadata

    save_json(
        out / "protocolo.json",
        dict(
            seed=SEED,
            posterior_draws=DRAWS,
            grid=CONFIGS,
            dataset=args.dataset,
            columns=list(original.columns),
            validation=(
                "5 outer grouped folds; 3 inner grouped folds; "
                "one hidden cell per query"
            ),
            versions={
                p: importlib.metadata.version(p)
                for p in ["numpy", "pandas", "scipy", "scikit-learn"]
            },
            hashes={
                p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
                for p in [
                    "estudio_ampliado.py",
                    "data/concentraciones_mM.db",
                    "config_normalizacion.json",
                ]
            },
        ),
    )
    print(macro, flush=True)


if __name__ == "__main__":
    main()
