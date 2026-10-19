"""Verifica resultados de una ejecución independiente sin reentrenar modelos."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', type=Path, default=ROOT/'resultados/ishii')
    out = parser.parse_args().results
    manifest = json.loads((out/'manifest.json').read_text())
    dataset = manifest['dataset']
    for name, expected in manifest['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected, name
    for name, expected in manifest['output_sha256'].items():
        assert hashlib.sha256((out/name).read_bytes()).hexdigest() == expected, name
    data = pd.read_csv(ROOT/'normalizacion'/f'matriz_{dataset}_R.csv', index_col=0)
    mm = pd.read_csv(ROOT/'normalizacion'/f'matriz_{dataset}_mM.csv', index_col=0)
    meta = pd.read_csv(ROOT/'normalizacion/muestras_y_exclusiones.csv', index_col=0)
    refs = pd.read_csv(ROOT/'normalizacion/referencias_glucosa.csv', index_col=0).value_mM
    fits = json.loads((out/'ajustes.json').read_text())
    for fit in fits:
        train = data.drop(index=fit['excluded'])
        assert list(train.index) == fit['train_ids']
        assert not set(fit['excluded']) & set(train.index)
        if fit['scaler_mean'] is not None:
            np.testing.assert_allclose(fit['scaler_mean'], train[fit['active_features']].mean(), rtol=1e-12, atol=1e-14)
    masks = json.loads((out/'validacion/mascaras.json').read_text())
    trials = {trial['trial_id']:trial for trial in masks['trials']}
    pred = pd.read_csv(out/'validacion/predicciones.csv')
    truth = pd.read_csv(out/'validacion/verdades.csv')
    for _, group in pred.groupby('method'):
        assert set(group.case_id) == set(truth.case_id)
    if dataset == 'mmc2':
        assert 'MMC2_GLUCOSE' not in set(pred.sample_id)
        assert 'analytical_series' not in set(pred.experiment)
    for row in pred.itertuples():
        trial, fit = trials[row.trial_id], fits[row.fit_id]
        assert row.sample_id in trial['excluded']
        assert row.metabolite_id in trial['hidden'] and row.metabolite_id in trial['score']
        assert fit['excluded'] == sorted(trial['excluded']) and fit['method'] == row.method
        assert np.isclose(row.truth_mM, mm.loc[row.sample_id, row.metabolite_id])
        assert np.isclose(row.truth_R, data.loc[row.sample_id, row.metabolite_id])
        assert np.isclose(row.absolute_error_mM, abs(row.estimate_mM-row.truth_mM), equal_nan=True)
        assert np.isclose(row.estimate_mM, row.estimate_R*refs[row.metabolite_id], equal_nan=True)
    candidates = pd.read_csv(out/'imputacion/candidatos_y_soporte.csv')
    for row in candidates.itertuples():
        assert pd.isna(data.loc[row.sample_id, row.metabolite_id])
        assert row.metabolite_id in manifest['config']['targets']
        excluded = sorted(data.index[meta.loc[data.index, 'perturbation_group'] == meta.loc[row.sample_id, 'perturbation_group']])
        assert fits[row.fit_id]['excluded'] == excluded
    for path in (out/'imputacion').glob('*_mM.csv'):
        matrix = pd.read_csv(path, index_col=0)
        assert matrix.index.equals(mm.index) and matrix.columns.equals(mm.columns)
        observed = mm.notna().to_numpy()
        np.testing.assert_allclose(matrix.to_numpy()[observed], mm.to_numpy()[observed], rtol=1e-12, atol=1e-14)
    if not candidates.empty:
        stats = pd.read_csv(out/'imputacion/estadisticas_antes_despues.csv')
        assert ((stats.observed+stats.imputed+stats.pending) == len(data)).all()
        common = pd.read_csv(out/'imputacion/casillas_comunes.csv')
        for _, part in common.groupby('policy'):
            sets = [set(zip(g.loc[g.common_methods, 'sample_id'], g.loc[g.common_methods, 'metabolite_id'])) for _, g in part.groupby('method')]
            assert all(s == sets[0] for s in sets)
    print(f'OK {dataset}: una sola base, máscaras comunes, train exclusivo, errores y observados intactos.')

if __name__ == '__main__':
    main()
