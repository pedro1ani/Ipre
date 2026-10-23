"""Controles de fuga de datos, observados e incertidumbre."""

import unittest

import numpy as np
import pandas as pd
from estudio_ampliado import CONFIGS, METHODS, cases, inner_selection, predict


class AmpliacionTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(5)
        self.train = pd.DataFrame(
            rng.uniform(0.2, 2, (12, 4)),
            columns=list("abcd"),
            index=[f"S{i}" for i in range(12)],
        )
        self.train.iloc[0, 1] = np.nan
        self.query = pd.DataFrame(
            [[np.nan, 0.7, 1.0, 1.4]], columns=list("abcd")
        )

    def test_hidden_truth_is_outside_query(self):
        query, info = cases(self.train)
        self.assertEqual(len(info), int(self.train.notna().sum().sum()))
        for row in info.itertuples():
            self.assertTrue(pd.isna(query.loc[row.query, row.metabolite]))
            self.assertEqual(
                row.truth_R, self.train.loc[row.sample, row.metabolite]
            )

    def test_preserves_observed_and_nonnegative(self):
        for method in METHODS:
            pred, _ = predict(
                self.train, self.query, method, CONFIGS[method][0]
            )
            np.testing.assert_array_equal(
                pred[["b", "c", "d"]], self.query[["b", "c", "d"]]
            )
            self.assertTrue((pred.to_numpy() >= 0).all())

    def test_unobserved_target_abstains(self):
        self.train["a"] = np.nan
        for method in METHODS:
            pred, _ = predict(
                self.train, self.query, method, CONFIGS[method][0]
            )
            self.assertTrue(pd.isna(pred.loc[0, "a"]))

    def test_query_batch_does_not_retrain(self):
        for method in METHODS:
            config = CONFIGS[method][0]
            first, _ = predict(self.train, self.query, method, config)
            batch = pd.concat(
                [self.query, self.query * 100], ignore_index=True
            )
            second, _ = predict(self.train, batch, method, config)
            np.testing.assert_allclose(first.iloc[0], second.iloc[0])

    def test_posterior_reproducible_but_not_constant(self):
        config = CONFIGS["Iterative"][0]
        a, _ = predict(self.train, self.query, "Iterative", config, True, 11)
        b, _ = predict(self.train, self.query, "Iterative", config, True, 11)
        c, _ = predict(self.train, self.query, "Iterative", config, True, 12)
        np.testing.assert_array_equal(a, b)
        self.assertNotEqual(a.loc[0, "a"], c.loc[0, "a"])

    def test_inner_groups_disjoint(self):
        groups = pd.Series(np.repeat(np.arange(6), 2), index=self.train.index)
        refs = pd.Series(1.0, index=self.train.columns)
        _, _, _, splits = inner_selection(self.train, groups, refs, 0)
        for split in splits:
            self.assertFalse(
                set(groups.loc[split["train"]])
                & set(groups.loc[split["validation"]])
            )


if __name__ == "__main__":
    unittest.main()
