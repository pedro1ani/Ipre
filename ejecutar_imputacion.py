<<<<<<< HEAD
"""Un estudio por ejecución: 1) evaluar con máscaras; 2) imputar candidatos reales.

Uso: python ejecutar_imputacion.py --dataset ishii --output resultados/ishii
MMC2 se procesa por separado cambiando únicamente --dataset y --output.
"""
=======
"""Imputa un único estudio por ejecución"""
>>>>>>> 95ef1e4dd97b7e57a5354f4487c6e4a147b94d28
from pathlib import Path
import argparse
import hashlib
import importlib.metadata
import json
import platform
import sqlite3
import numpy as np
import pandas as pd
from imputadores import AuditedImputer, METHODS

ROOT = Path(__file__).resolve().parent
POLICIES = ("available", "supported_plus_baseline")

def require(condition, message):
    if not condition:
        raise ValueError(message)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def save(path, frame):
    frame.to_csv(path, index=False, encoding="utf-8-sig", float_format="%.17g")

def dump(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")

def validate_input(data, samples):
    require(data.index.is_unique and data.columns.is_unique, "IDs duplicados")
    require(data.index.isin(samples.index).all(), "Muestra sin metadatos")

def accepted(frame, policy):
    finite = np.isfinite(frame.estimate_R)
    if policy == "available":
        return finite
    return finite & ((frame.method == "SimpleImputer") | (frame.support == "supported"))

def make_masks(ish,mm,metadata,series,targets,fingerprint, dataset="ishii"):
 trials=[];truth=[];pattern_rows=[]
 def add(experiment,sample,hidden,score,excluded,pattern=''):
  tid=f'{experiment}:{sample}:{pattern or score[0]}'
  for met in score:
   require(pd.notna(ish.loc[sample,met]),'Verdad no observada')
   truth.append(dict(case_id=tid+'|'+met,trial_id=tid,sample_id=sample,metabolite_id=met,experiment=experiment,truth_R=float(ish.loc[sample,met]),truth_mM=float(mm.loc[sample,met])))
  trials.append(dict(trial_id=tid,experiment=experiment,sample_id=sample,hidden=hidden,score=score,excluded=sorted(excluded),pattern_id=pattern))
 for sample in ish.index:
  group=metadata.loc[sample,'perturbation_group'];excluded=list(ish.index[metadata.loc[ish.index,'perturbation_group']==group])
  if dataset=='mmc2' and sample=='MMC2_GLUCOSE':continue
  if series is not None:
   family=series.loc[sample,'series_family'];series_excluded=list(ish.index[(series.loc[ish.index,'series_family']==family)|(metadata.loc[ish.index,'perturbation_group']==group)])
  for met in targets:
   if pd.notna(ish.loc[sample,met]):
    add('perturbation',sample,[met],[met],excluded)
    if series is not None:add('analytical_series',sample,[met],[met],series_excluded)
 patterns={}
 for sample in ish.index[ish[targets].isna().any(axis=1)]:
  key=tuple(ish.columns[ish.loc[sample].isna()]);patterns.setdefault(key,[]).append(sample)
 for k,(hidden,origins) in enumerate(sorted(patterns.items()),1):
  pid=f'P{k:02d}';hidden=list(hidden);score=[m for m in targets if m in hidden];visible=[m for m in ish if m not in hidden]
  eligible=ish.index[ish[score].notna().all(axis=1)&ish[visible].notna().all(axis=1)]
  pattern_rows.append(dict(pattern_id=pid,source_samples=list(origins),hidden=hidden,score=score,eligible_samples=list(eligible)))
  for sample in eligible:
   if dataset=='mmc2' and sample=='MMC2_GLUCOSE':continue
   group=metadata.loc[sample,'perturbation_group'];excluded=list(ish.index[metadata.loc[ish.index,'perturbation_group']==group])
   add('real_patterns',sample,hidden,score,excluded,pid)
 trials.sort(key=lambda x:x['trial_id']);truth.sort(key=lambda x:x['case_id'])
 require(len({t['trial_id'] for t in trials})==len(trials),'Trials duplicados')
 return dict(schema=1,fingerprint=fingerprint,feature_order=list(ish.columns),sample_order=list(ish.index),trials=trials,patterns=pattern_rows),pd.DataFrame(truth)


def evaluate_metrics(predictions):
    """Compara errores propios y sobre las mismas celdas para los tres métodos."""
    rows, reasons, memberships = [], [], []
    for policy in POLICIES:
        sets = {method: set(part.loc[accepted(part, policy), 'case_id'])
                for method, part in predictions.groupby('method')}
        common = set.intersection(*(sets[method] for method in METHODS))
        for method, part in predictions.groupby('method'):
            for row in part.itertuples():
                responded = row.case_id in sets[method]
                memberships.append(dict(policy=policy, method=method, case_id=row.case_id,
                                        responded=responded, common_methods=row.case_id in common))
                reason = 'responded' if responded else f'{row.support}: {row.mechanism}'
                if not responded and row.method == 'IterativeImputer' and row.converged == False:
                    reason += '; not_converged'
                reasons.append(dict(policy=policy, method=method, experiment=row.experiment,
                                    metabolite_id=row.metabolite_id, reason=reason))
            for scope, selected in [('own', sets[method]), ('common_methods', common)]:
                for (experiment, met), group in part.groupby(['experiment', 'metabolite_id']):
                    keep = group[group.case_id.isin(selected)]
                    rows.append(dict(policy=policy, scope=scope, method=method,
                        experiment=experiment, metabolite_id=met, attempts=len(group),
                        predictions=len(keep), coverage_pct=100 * len(keep) / len(group),
                        MAE_mM=keep.absolute_error_mM.mean(), MAE_R=keep.absolute_error_R.mean()))
    reasons = pd.DataFrame(reasons).value_counts(sort=False).rename('cells').reset_index()
    return pd.DataFrame(rows), reasons, pd.DataFrame(memberships)


def summarize_candidates(original_mM, candidates, directory):
    """Conserva los observados; completa solo objetivos candidatos y calcula estadísticas."""
    stats, membership = [], []
    for policy in POLICIES:
        sets = {method: set(zip(part.loc[accepted(part, policy), 'sample_id'],
                               part.loc[accepted(part, policy), 'metabolite_id']))
                for method, part in candidates.groupby('method')}
        common = set.intersection(*(sets[method] for method in METHODS))
        for method in METHODS:
            part = candidates[candidates.method == method]
            for row in part.itertuples():
                key = (row.sample_id, row.metabolite_id)
                membership.append(dict(policy=policy, method=method, sample_id=row.sample_id,
                    metabolite_id=row.metabolite_id, selected=key in sets[method], common_methods=key in common))
            for scope, selected in [('own', sets[method]), ('common_methods', common)]:
                completed = original_mM.copy()
                for row in part.itertuples():
                    if (row.sample_id, row.metabolite_id) in selected:
                        require(pd.isna(completed.loc[row.sample_id, row.metabolite_id]), 'Observado modificado')
                        completed.loc[row.sample_id, row.metabolite_id] = row.estimate_mM
                if scope == 'own':
                    completed.to_csv(directory / f'{policy}_{method}_mM.csv', encoding='utf-8-sig', float_format='%.17g')
                for stage, matrix in [('before', original_mM), ('after', completed)]:
                    for met in matrix:
                        values = matrix[met].dropna()
                        mean, sd = values.mean(), values.std(ddof=1)
                        status = ('n_less_than_2' if len(values) < 2 else
                                  'mean_nonpositive' if mean <= 0 else
                                  'mean_near_zero' if abs(mean) <= 1e-12 else 'interpretable')
                        observed = int(original_mM[met].count())
                        stats.append(dict(policy=policy, scope=scope, method=method, stage=stage,
                            metabolite_id=met, total_samples=len(matrix), observed=observed,
                            imputed=len(values)-observed, pending=len(matrix)-len(values),
                            mean_mM=mean, median_mM=values.median(), SD_mM=sd,
                            IQR_mM=values.quantile(.75)-values.quantile(.25),
                            CV_pct=100*sd/mean if status == 'interpretable' else np.nan, CV_status=status))
    return pd.DataFrame(stats), pd.DataFrame(membership)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset', choices=['ishii', 'mmc2'], default='ishii')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    dataset = args.dataset
    out = (args.output or ROOT / 'resultados' / dataset).resolve()
    require(not out.exists(), 'Utilice una carpeta de salida nueva')

    # PREPARACIÓN: una sola base seleccionada; nunca se concatenan Ishii y MMC2.
    config_path = ROOT / f'config_imputacion_{dataset}.json'
    cfg = json.loads(config_path.read_text())
    targets = cfg['targets']
    matrix_path = ROOT / 'normalizacion' / f'matriz_{dataset}_R.csv'
    mm_path = ROOT / 'normalizacion' / f'matriz_{dataset}_mM.csv'
    samples_path = ROOT / 'normalizacion/muestras_y_exclusiones.csv'
    refs_path = ROOT / 'normalizacion/referencias_glucosa.csv'
    data = pd.read_csv(matrix_path, index_col=0)
    original_mM = pd.read_csv(mm_path, index_col=0)
    samples = pd.read_csv(samples_path, index_col=0)
    refs = pd.read_csv(refs_path, index_col=0).value_mM
    validate_input(data, samples)
    require(cfg['dataset'] == dataset and tuple(cfg['methods']) == METHODS, 'Configuración incorrecta')
    require(bool(targets) and len(set(targets)) == len(targets) and set(targets) <= set(data.columns), 'Objetivos inválidos')
    require(data.index.equals(original_mM.index) and data.columns.equals(original_mM.columns), 'Orden inconsistente')
    require(np.allclose(data.mul(refs).to_numpy(), original_mM.to_numpy(), equal_nan=True), 'mM y ratios inconsistentes')
    paths = [config_path, matrix_path, mm_path, samples_path, refs_path, ROOT/'imputadores.py', Path(__file__)]
    series = None
    if dataset == 'ishii':
        source = ROOT / 'sources/ecoli_auditada_v3.db'
        paths.append(source)
        with sqlite3.connect(source.resolve().as_uri() + '?mode=ro', uri=True) as connection:
            series = pd.read_sql_query("SELECT sample_id_original AS sample_id, source_series_id AS series_raw, source_column FROM ishii_metabolite_source_columns WHERE source_sheet='metabolites'", connection).set_index('sample_id').loc[data.index]
        series['series_family'] = series.series_raw.astype(str).str.replace('*', '', regex=False)
    hashes = {str(path.relative_to(ROOT)): digest(path) for path in paths}
    collection, truth = make_masks(data, original_mM, samples, series, targets, hashes, dataset)
    require(not truth.empty, 'No hay valores conocidos elegibles para evaluación')
    out.mkdir(parents=True)
    validation, imputation = out/'validacion', out/'imputacion'
    validation.mkdir()
    imputation.mkdir()
    dump(validation/'mascaras.json', collection)
    save(validation/'verdades.csv', truth)
    if series is not None:
        save(validation/'series_analiticas.csv', series.reset_index())

    models, fits = {}, []
    def get_model(method, excluded):
        key = (method, tuple(sorted(excluded)))
        if key not in models:
            train = data.drop(index=excluded)
            require(not set(train.index) & set(excluded), 'Filtración de muestras')
            model = AuditedImputer(method, cfg).fit(train)
            model.audit_fit_id = len(fits)
            fits.append(dict(fit_id=len(fits), excluded=sorted(excluded), **model.metadata()))
            models[key] = model
        return models[key]

    # SECCIÓN 1 — EVALUAR: ocultar valores conocidos y comparar estimado vs medido.
    predictions = []
    truths = truth.set_index('case_id')
    for trial in collection['trials']:
        query = data.loc[[trial['sample_id']]].copy()
        query.loc[:, trial['hidden']] = np.nan
        for method in METHODS:
            model = get_model(method, trial['excluded'])
            _, diagnostics = model.transform(query)
            for row in diagnostics[diagnostics.metabolite_id.isin(trial['score'])].to_dict('records'):
                case_id = trial['trial_id'] + '|' + row['metabolite_id']
                known = truths.loc[case_id]
                row.update(fit_id=model.audit_fit_id, case_id=case_id, trial_id=trial['trial_id'],
                    experiment=trial['experiment'], truth_R=float(known.truth_R), truth_mM=float(known.truth_mM),
                    estimate_mM=row['estimate_R'] * refs[row['metabolite_id']])
                row['absolute_error_R'] = abs(row['estimate_R'] - row['truth_R'])
                row['absolute_error_mM'] = abs(row['estimate_mM'] - row['truth_mM'])
                predictions.append(row)
    predictions = pd.DataFrame(predictions)
    require(len(predictions) == len(truth)*len(METHODS), 'Evaluación incompleta')
    save(validation/'predicciones.csv', predictions)
    metrics, reasons, common = evaluate_metrics(predictions)
    save(validation/'metricas.csv', metrics)
    save(validation/'cobertura_y_motivos.csv', reasons)
    save(validation/'casos_comunes.csv', common)
    dump(validation/'estado.json', dict(status='COMPLETE', dataset=dataset, scored_cells=len(truth)))
    print(f'{dataset}: evaluación completa, {len(truth)} celdas por método.', flush=True)

    # SECCIÓN 2 — IMPUTAR: aplicar los métodos a los huecos reales identificados.
    # Son candidatos solo los NaN de los cuatro objetivos; no se modifican observados.
    candidate_rows = []
    for sample in data.index[data[targets].isna().any(axis=1)]:
        group = samples.loc[sample, 'perturbation_group']
        excluded = list(data.index[samples.loc[data.index, 'perturbation_group'] == group])
        query = data.loc[[sample]].copy()
        for method in METHODS:
            model = get_model(method, excluded)
            _, diagnostics = model.transform(query)
            for row in diagnostics[diagnostics.metabolite_id.isin(targets)].to_dict('records'):
                row.update(fit_id=model.audit_fit_id, excluded_group=group,
                           estimate_mM=row['estimate_R'] * refs[row['metabolite_id']])
                candidate_rows.append(row)
    # Columnas explícitas permiten también ejecutar una base sin candidatos reales.
    candidates = pd.DataFrame(candidate_rows, columns=list(dict.fromkeys(
        list(predictions.columns) + ['excluded_group'])))
    save(imputation/'candidatos_y_soporte.csv', candidates)
    if not candidates.empty:
        stats, membership = summarize_candidates(original_mM, candidates, imputation)
        save(imputation/'estadisticas_antes_despues.csv', stats)
        save(imputation/'casillas_comunes.csv', membership)
    dump(out/'ajustes.json', fits)
    require(hashes == {str(path.relative_to(ROOT)): digest(path) for path in paths}, 'Entradas modificadas')
    dump(out/'manifest.json', dict(dataset=dataset, status='SUCCESS', config=cfg,
        input_sha256=hashes, python=platform.python_version(),
        versions={name:importlib.metadata.version(name) for name in ['numpy','pandas','scipy','scikit-learn']},
        validation_attempts=len(predictions), real_attempts=len(candidates),
        unavailable_tests=[] if series is not None else ['analytical_series: no metadata homologous to Ishii'],
        excluded_evaluation_samples=[] if dataset=='ishii' else ['MMC2_GLUCOSE: reference reveals truth'],
        output_sha256={str(path.relative_to(out)):digest(path) for path in sorted(out.rglob('*')) if path.is_file()}))
    print(f'{dataset}: {len(candidates)} candidatos (todos los métodos). Salida: {out}', flush=True)

if __name__ == '__main__':
    main()
