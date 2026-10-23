"""Verifica las salidas guardadas sin volver a entrenar."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def main():
    metadata = pd.read_csv(ROOT / "normalizacion/muestras_y_exclusiones.csv")
    groups = metadata.set_index("sample_id").perturbation_group
    summary = {}
    for dataset in ["ishii", "mmc2"]:
        path = ROOT / "resultados_ampliados" / dataset
        original = pd.read_csv(path / f"matriz_{dataset}_mM.csv", index_col=0)
        predictions = pd.read_csv(path / "predicciones.csv")
        candidates = pd.read_csv(path / "candidatos.csv")
        draws = pd.read_csv(path / "muestras_posteriores.csv")
        partitions = json.loads((path / "particiones.json").read_text())
        assert original.shape[1] == 10
        expected = int(original.notna().sum().sum())
        if dataset == "mmc2":
            expected -= int(original.loc["MMC2_GLUCOSE"].notna().sum())
        assert len(predictions) == expected * 3
        assert not predictions.duplicated(
            ["sample", "metabolite", "method"]
        ).any()
        assert len(candidates) == int(original.isna().sum().sum()) * 3
        assert not candidates.duplicated(
            ["sample", "metabolite", "method"]
        ).any()
        for split in partitions:
            train, test = split["train"], split["test"]
            assert not set(train) & set(test)
            assert not set(groups.loc[train]) & set(groups.loc[test])
            for inner in split["inner_splits"]:
                assert set(inner["train"] + inner["validation"]) <= set(train)
                assert not set(groups.loc[inner["train"]]) & set(
                    groups.loc[inner["validation"]]
                )
            for table in [predictions, candidates]:
                assert set(
                    table[table.fold == split["fold"]]["sample"]
                ) <= set(test)
        for row in predictions.itertuples():
            assert np.isclose(
                row.truth_mM, original.loc[row.sample, row.metabolite]
            )
            assert np.isclose(
                row.absolute_error_mM, abs(row.estimate_mM - row.truth_mM)
            )
        for method in ["Mediana", "KNN", "Iterative"]:
            completed = pd.read_csv(
                path / f"imputada_{method}.csv", index_col=0
            )
            np.testing.assert_array_equal(
                completed.to_numpy()[original.notna()],
                original.to_numpy()[original.notna()],
            )
            assert (completed.stack() >= 0).all()
        assert (draws.value_mM >= 0).all()
        assert (
            draws.groupby(["sample", "metabolite", "real_missing"]).size()
            == 30
        ).all()
        for real, table in [(False, predictions), (True, candidates)]:
            table = table[table.method == "Iterative"]
            for row in table.itertuples():
                values = draws[
                    (draws["sample"] == row.sample)
                    & (draws.metabolite == row.metabolite)
                    & (draws.real_missing == real)
                ].value_mM
                assert np.isclose(values.std(), row.posterior_SD_mM)
                assert np.isclose(values.quantile(0.025), row.posterior_low_mM)
                assert np.isclose(
                    values.quantile(0.975), row.posterior_high_mM
                )
        summary[dataset] = dict(
            metabolites=10,
            masked_cells=expected,
            predictions=len(predictions),
            candidates=len(candidates),
            posterior_draws=len(draws),
            disjoint_groups=True,
            observed_unchanged=True,
            posterior_verified=True,
        )
    (ROOT / "VERIFICACION_AMPLIACION.json").write_text(
        json.dumps(summary, indent=2)
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
