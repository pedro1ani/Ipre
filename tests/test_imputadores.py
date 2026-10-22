import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from imputadores import METHODS, Imputador

CFG = json.loads(
    (
        Path(__file__).resolve().parents[1] / "config_imputacion.json"
    ).read_text()
)


class BehaviorTests(unittest.TestCase):
    def frame(self, rows):
        return pd.DataFrame(rows, columns=["a", "b", "c"])

    def test_median_is_raw_training_median(self):
        train = self.frame([[1, 2, 3], [5, 4, 6], [20, 8, 9], [30, 10, 12]])
        out, d = (
            Imputador("SimpleImputer", CFG)
            .fit(train)
            .transform(self.frame([[np.nan, 123, 456]]))
        )
        self.assertEqual(out.iloc[0, 0], 12.5)
        self.assertEqual(d.iloc[0].support, "training_only")
        self.assertEqual(out.iloc[0, 1], 123)

    def test_knn_no_overlap_mean_not_median(self):
        train = self.frame([[1, 2, 3], [2, 4, 6], [12, 8, 9]])
        out, d = (
            Imputador("KNNImputer", CFG)
            .fit(train)
            .transform(self.frame([[np.nan, np.nan, np.nan]]))
        )
        self.assertAlmostEqual(out.iloc[0, 0], 5)
        self.assertTrue((d.support == "training_only").all())
        self.assertTrue(d.mechanism.str.contains("training_mean").all())

    def test_knn_two_neighbors_flagged(self):
        train = self.frame([[1, 2, 3], [5, 4, 6]])
        out, d = (
            Imputador("KNNImputer", CFG)
            .fit(train)
            .transform(self.frame([[np.nan, 3, 4]]))
        )
        self.assertAlmostEqual(out.iloc[0, 0], 3)
        self.assertEqual(d.iloc[0].n_usable_neighbors, 2)
        self.assertEqual(d.iloc[0].support, "limited_support")

    def test_empty_column_never_zero(self):
        train = self.frame([[np.nan, 2, 3], [np.nan, 4, 6], [np.nan, 8, 9]])
        for m in METHODS:
            out, d = (
                Imputador(m, CFG)
                .fit(train)
                .transform(self.frame([[np.nan, 5, 7]]))
            )
            self.assertTrue(np.isnan(out.iloc[0, 0]))
            self.assertEqual(d.iloc[0].support, "abstention")

    def test_all_empty_training(self):
        for m in METHODS:
            out, d = (
                Imputador(m, CFG)
                .fit(self.frame([[np.nan] * 3]))
                .transform(self.frame([[np.nan] * 3]))
            )
            self.assertTrue(out.isna().all().all())
            self.assertTrue((d.support == "abstention").all())

    def test_training_only_scaling_and_preserved_observations(self):
        train = self.frame([[1, 2, 3], [2, 4, 6], [3, 8, 9]])
        for m in METHODS:
            model = Imputador(m, CFG).fit(train)
            before = train.copy()
            query = self.frame([[np.nan, 1e8, 3]])
            out, d = model.transform(query)
            pd.testing.assert_frame_equal(train, before)
            self.assertEqual(out.iloc[0, 1], 1e8)
            self.assertEqual(out.iloc[0, 2], 3)
            if model.escalador is not None:
                np.testing.assert_allclose(model.escalador.mean_, train.mean())

    def test_iterative_no_observations_not_supported(self):
        train = self.frame([[1, 2, 3], [2, 4, 6], [3, 8, 9]])
        out, d = (
            Imputador("IterativeImputer", CFG)
            .fit(train)
            .transform(self.frame([[np.nan] * 3]))
        )
        self.assertTrue((d.support == "training_only").all())
        self.assertTrue(np.isfinite(out).all().all())

    def test_nonconvergence_is_visible(self):
        cfg = json.loads(json.dumps(CFG))
        cfg["rondas_maximas"] = 1
        cfg["tolerancia"] = 1e-15
        train = self.frame(
            [[1, 2, np.nan], [2, np.nan, 5], [3, 7, 8], [4, 8, 9], [5, 12, 15]]
        )
        model = Imputador("IterativeImputer", cfg).fit(train)
        self.assertFalse(model.convergio)
        self.assertTrue(model.advertencias)
        _, d = model.transform(self.frame([[np.nan, 6, 7]]))
        self.assertEqual(d.iloc[0].support, "limited_support")

    def test_order_guard(self):
        model = Imputador("SimpleImputer", CFG).fit(self.frame([[1, 2, 3]]))
        with self.assertRaises(ValueError):
            model.transform(pd.DataFrame([[1, 2, 3]], columns=["b", "a", "c"]))

    def test_deterministic(self):
        train = self.frame(
            [[1, 2, np.nan], [2, 4, 5], [3, 7, 8], [4, np.nan, 9]]
        )
        for m in METHODS:
            a, _ = (
                Imputador(m, CFG)
                .fit(train)
                .transform(self.frame([[np.nan, 6, 7]]))
            )
            b, _ = (
                Imputador(m, CFG)
                .fit(train)
                .transform(self.frame([[np.nan, 6, 7]]))
            )
            pd.testing.assert_frame_equal(a, b)


class InputChecksTests(unittest.TestCase):
    def test_mixed_origin_not_blocked(self):
        from ejecutar_imputacion import validar_entrada

        data = pd.DataFrame({"a": [1, 2]}, index=["KO01", "MMC2_GLUCOSE"])
        meta = pd.DataFrame(
            {"dataset_id": ["ISHII2007", "MMC2"]}, index=data.index
        )
        validar_entrada(data, meta)

    def test_each_dataset_allowed_separately(self):
        from ejecutar_imputacion import validar_entrada

        meta = pd.DataFrame(
            {"dataset_id": ["ISHII2007", "MMC2"]},
            index=["KO01", "MMC2_GLUCOSE"],
        )
        for dataset, key in [("ishii", "KO01"), ("mmc2", "MMC2_GLUCOSE")]:
            validar_entrada(pd.DataFrame({"a": [1]}, index=[key]), meta)


if __name__ == "__main__":
    unittest.main()
