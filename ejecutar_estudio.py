"""Máscaras únicas -> validación común -> candidatos reales -> estadísticas. Sin log."""
from pathlib import Path
import argparse,csv,hashlib,json,sqlite3,platform,importlib.metadata
import numpy as np
import pandas as pd
from imputadores import AuditedImputer,METHODS
ROOT=Path(__file__).resolve().parent
SCENARIOS=('ishii','ishii_mmc2')
POLICIES=('available','supported_plus_baseline')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def save(p,d):d.to_csv(p,index=False,encoding='utf-8-sig',float_format='%.17g')
def require(ok,message):
 if not ok:raise ValueError(message)

def make_masks(ish,mm,metadata,series,targets,fingerprint):
 trials=[];truth=[];pattern_rows=[]
 def add(experiment,sample,hidden,score,excluded,pattern=''):
  tid=f'{experiment}:{sample}:{pattern or score[0]}'
  for met in score:
   require(pd.notna(ish.loc[sample,met]),'Verdad no observada')
   truth.append(dict(case_id=tid+'|'+met,trial_id=tid,sample_id=sample,metabolite_id=met,experiment=experiment,truth_R=float(ish.loc[sample,met]),truth_mM=float(mm.loc[sample,met])))
  trials.append(dict(trial_id=tid,experiment=experiment,sample_id=sample,hidden=hidden,score=score,excluded=sorted(excluded),pattern_id=pattern))
 for sample in ish.index:
  group=metadata.loc[sample,'perturbation_group'];excluded=list(ish.index[metadata.loc[ish.index,'perturbation_group']==group])
  family=series.loc[sample,'series_family'];series_excluded=list(ish.index[(series.loc[ish.index,'series_family']==family)|(metadata.loc[ish.index,'perturbation_group']==group)])
  for met in targets:
   if pd.notna(ish.loc[sample,met]):
    add('perturbation',sample,[met],[met],excluded)
    add('analytical_series',sample,[met],[met],series_excluded)
 patterns={}
 for sample in ish.index[ish[targets].isna().any(axis=1)]:
  key=tuple(ish.columns[ish.loc[sample].isna()]);patterns.setdefault(key,[]).append(sample)
 for k,(hidden,origins) in enumerate(sorted(patterns.items()),1):
  pid=f'P{k:02d}';hidden=list(hidden);score=[m for m in targets if m in hidden];visible=[m for m in ish if m not in hidden]
  eligible=ish.index[ish[score].notna().all(axis=1)&ish[visible].notna().all(axis=1)]
  pattern_rows.append(dict(pattern_id=pid,source_samples=list(origins),hidden=hidden,score=score,eligible_samples=list(eligible)))
  for sample in eligible:
   group=metadata.loc[sample,'perturbation_group'];excluded=list(ish.index[metadata.loc[ish.index,'perturbation_group']==group])
   add('real_patterns',sample,hidden,score,excluded,pid)
 trials.sort(key=lambda x:x['trial_id']);truth.sort(key=lambda x:x['case_id'])
 require(len({t['trial_id'] for t in trials})==len(trials),'Trials duplicados')
 return dict(schema=1,fingerprint=fingerprint,feature_order=list(ish.columns),sample_order=list(ish.index),trials=trials,patterns=pattern_rows),pd.DataFrame(truth)

def eligible(d,policy):
 finite=np.isfinite(d.estimate_R.to_numpy(float))
 return pd.Series(finite,index=d.index) if policy=='available' else pd.Series(finite,index=d.index)&((d.method=='SimpleImputer')|(d.support=='supported'))

def rejection(row,policy):
 if not np.isfinite(row.estimate_R):return row.mechanism
 if policy=='available' or row.method=='SimpleImputer' or row.support=='supported':return 'responded'
 reasons=[row.support,row.mechanism]
 if row.method=='IterativeImputer' and row.converged is False:reasons.append('not_converged')
 return ';'.join(reasons)

def metrics(pred):
 rows=[];membership=[];delta=[];coverage=[]
 for policy in POLICIES:
  d=pred.copy();d['responded']=eligible(d,policy)
  sets={(sc,m):set(g.loc[g.responded,'case_id']) for (sc,m),g in d.groupby(['scenario','method'])}
  all6=set.intersection(*(sets[(sc,m)] for sc in SCENARIOS for m in METHODS))
  same_sc={sc:set.intersection(*(sets[(sc,m)] for m in METHODS)) for sc in SCENARIOS}
  same_method={m:sets[('ishii',m)]&sets[('ishii_mmc2',m)] for m in METHODS}
  for (sc,m),g in d.groupby(['scenario','method']):
   for scope,ids in [('own',sets[(sc,m)]),('common_methods',same_sc[sc]),('common_all_six',all6)]:
    for (experiment,met),part in g.groupby(['experiment','metabolite_id']):
     keep=part[part.case_id.isin(ids)];n=len(keep)
     rows.append(dict(policy=policy,scope=scope,scenario=sc,method=m,experiment=experiment,metabolite_id=met,attempts=len(part),predictions=n,coverage_pct=100*n/len(part),MAE_mM=float(keep.absolute_error_mM.mean()) if n else np.nan,MAE_R=float(keep.absolute_error_R.mean()) if n else np.nan,bias_mM=float((keep.estimate_mM-keep.truth_mM).mean()) if n else np.nan))
   for row in g.itertuples():
    membership.append(dict(policy=policy,case_id=row.case_id,scenario=sc,method=m,responded=bool(row.responded),common_methods=row.case_id in same_sc[sc],common_all_six=row.case_id in all6,common_scenarios_same_method=row.case_id in same_method[m]))
   g['reason']=[rejection(x,policy) for x in g.itertuples()]
   for (experiment,met,reason),part in g.groupby(['experiment','metabolite_id','reason']):coverage.append(dict(policy=policy,scenario=sc,method=m,experiment=experiment,metabolite_id=met,reason=reason,cells=len(part)))
  for m,ids in same_method.items():
   a=d[(d.scenario=='ishii')&(d.method==m)&d.case_id.isin(ids)]
   b=d[(d.scenario=='ishii_mmc2')&(d.method==m)&d.case_id.isin(ids)]
   joined=a.merge(b,on='case_id',suffixes=('_ishii','_mixed'),validate='one_to_one')
   for (experiment,met),g in d[(d.scenario=='ishii')&(d.method==m)].groupby(['experiment','metabolite_id']):
    sub=joined[(joined.experiment_ishii==experiment)&(joined.metabolite_id_ishii==met)]
    delta.append(dict(policy=policy,method=m,experiment=experiment,metabolite_id=met,paired_cells=len(sub),MAE_ishii_mM=sub.absolute_error_mM_ishii.mean(),MAE_ishii_mmc2_mM=sub.absolute_error_mM_mixed.mean(),delta_MAE_mM=(sub.absolute_error_mM_mixed-sub.absolute_error_mM_ishii).mean(),delta_MAE_R=(sub.absolute_error_R_mixed-sub.absolute_error_R_ishii).mean()))
 return pd.DataFrame(rows),pd.DataFrame(coverage),pd.DataFrame(membership),pd.DataFrame(delta)

def statistics(ish_mM,real,targets,epsilon):
 rows=[];members=[];matrices={};long=[]
 baseline={}
 def desc(x):
  x=x.dropna();n=len(x);mean=x.mean();sd=x.std(ddof=1)
  reason='n_less_than_2' if n<2 else ('mean_nonpositive' if mean<=0 else ('mean_near_zero' if abs(mean)<=epsilon else 'interpretable'))
  return dict(mean_mM=mean,median_mM=x.median(),SD_mM=sd,Q1_mM=x.quantile(.25),Q3_mM=x.quantile(.75),IQR_mM=x.quantile(.75)-x.quantile(.25),CV_pct=100*sd/mean if reason=='interpretable' else np.nan,CV_status=reason)
 for met in ish_mM:baseline[met]=desc(ish_mM[met])
 for policy in POLICIES:
  sets={(sc,m):set(zip(g.loc[eligible(g,policy),'sample_id'],g.loc[eligible(g,policy),'metabolite_id'])) for (sc,m),g in real.groupby(['scenario','method'])}
  common_sc={sc:set.intersection(*(sets[(sc,m)] for m in METHODS)) for sc in SCENARIOS}
  common6=set.intersection(*(sets[(sc,m)] for sc in SCENARIOS for m in METHODS))
  for (sc,method),g in real.groupby(['scenario','method']):
   for scope,chosen in [('own',sets[(sc,method)]),('common_methods',common_sc[sc]),('common_all_six',common6)]:
    matrix=ish_mM.copy()
    for v in g.itertuples():
     if (v.sample_id,v.metabolite_id) in chosen:
      require(pd.isna(matrix.loc[v.sample_id,v.metabolite_id]),'Intento de sustituir observado');matrix.loc[v.sample_id,v.metabolite_id]=v.estimate_mM
    observed=ish_mM.notna();require(np.array_equal(matrix.to_numpy()[observed],ish_mM.to_numpy()[observed]),'Original modificado')
    if scope=='own':
     matrices[(policy,sc,method)]=matrix
     for v in g.itertuples():members.append(dict(policy=policy,scenario=sc,method=method,sample_id=v.sample_id,metabolite_id=v.metabolite_id,selected=(v.sample_id,v.metabolite_id) in chosen,common_methods=(v.sample_id,v.metabolite_id) in common_sc[sc],common_all_six=(v.sample_id,v.metabolite_id) in common6))
    for met in matrix:
     n_original=int(ish_mM[met].notna().sum());n_final=int(matrix[met].notna().sum())
     for stage,values,imputed,pending in [('before',baseline[met],0,30-n_original),('after',desc(matrix[met]),n_final-n_original,30-n_final)]:
      rows.append(dict(policy=policy,scope=scope,scenario=sc,method=method,metabolite_id=met,stage=stage,total_samples=30,observed=n_original,imputed=imputed,pending=pending,**values))
 return pd.DataFrame(rows),pd.DataFrame(members),matrices

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=ROOT/'resultados_estudio');parser.add_argument('--masks',type=Path,default=ROOT/'validacion_comun');args=parser.parse_args()
 out=args.output.resolve();require(not out.exists(),'La salida debe ser nueva')
 cfg=json.loads((ROOT/'config_imputacion_ishii.json').read_text());vcfg=json.loads((ROOT/'config_validacion.json').read_text());targets=cfg['targets']
 require(vcfg['scenarios']==list(SCENARIOS) and vcfg['methods']==list(METHODS),'Escenarios/métodos fuera del protocolo');require(cfg['seed']==vcfg['seed'],'Semillas incoherentes')
 ish=pd.read_csv(ROOT/'normalizacion/matriz_ishii_R.csv',index_col=0);mm=pd.read_csv(ROOT/'normalizacion/matriz_ishii_mM.csv',index_col=0);ger=pd.read_csv(ROOT/'normalizacion/matriz_mmc2_R.csv',index_col=0)
 metadata=pd.read_csv(ROOT/'normalizacion/muestras_y_exclusiones.csv',index_col=0);refs=pd.read_csv(ROOT/'normalizacion/referencias_glucosa.csv',index_col=0).value_mM
 require(len(ish)==30 and list(ish.columns)==list(ger.columns)==list(mm.columns) and list(ish.index)==list(mm.index),'Orden/matriz incorrectos');require(np.allclose(ish.mul(refs).to_numpy(),mm.to_numpy(),equal_nan=True),'mM/ratio incoherente')
 source=ROOT/'sources/ecoli_auditada_v3.db'
 with sqlite3.connect(source.resolve().as_uri()+'?mode=ro',uri=True) as c:
  series=pd.read_sql_query("SELECT sample_id_original AS sample_id,source_series_id AS series_raw,source_column FROM ishii_metabolite_source_columns WHERE source_sheet='metabolites'",c).set_index('sample_id').loc[ish.index].copy()
 require(series.series_raw.notna().all(),'Serie desconocida');series['series_family']=series.series_raw.str.replace('*','',regex=False)
 input_paths=[ROOT/'normalizacion'/f for f in ['matriz_ishii_R.csv','matriz_ishii_mM.csv','matriz_mmc2_R.csv','muestras_y_exclusiones.csv','referencias_glucosa.csv']]+[source,ROOT/'config_validacion.json',ROOT/'config_imputacion_ishii.json',ROOT/'imputadores.py',Path(__file__).resolve()]
 before={str(p.relative_to(ROOT)):sha(p) for p in input_paths}
 fingerprint={k:v for k,v in before.items() if not k.endswith('.py')}
 collection,truth=make_masks(ish,mm,metadata,series,targets,fingerprint)
 masks=args.masks.resolve()
 if masks.exists():
  require(json.loads((masks/'masks.json').read_text())==collection,'Máscaras incompatibles con matriz/configuración actuales; use nueva carpeta de máscaras')
  stored=pd.read_csv(masks/'verdades_ishii.csv');pd.testing.assert_frame_equal(stored,truth,check_exact=False,rtol=1e-14,atol=1e-15)
 else:
  masks.mkdir(parents=True);dump(masks/'masks.json',collection);save(masks/'verdades_ishii.csv',truth);save(masks/'series_analiticas.csv',series.reset_index())
 out.mkdir(parents=True);vdir=out/'validacion';vdir.mkdir();idir=out/'imputacion';sdir=out/'estadisticas'
 for p in [idir,sdir]:p.mkdir()
 cache={};fit_records=[]
 def get_model(sc,method,excluded):
  key=(sc,method,tuple(sorted(excluded)))
  if key not in cache:
   train=ish.drop(index=excluded)
   if sc=='ishii_mmc2':train=pd.concat([train,ger])
   require(not set(train.index)&set(excluded),'Filtración de filas')
   estimator=AuditedImputer(method,cfg).fit(train);estimator.audit_fit_id=len(fit_records);cache[key]=estimator
   fit_records.append(dict(fit_id=len(fit_records),scenario=sc,excluded=list(excluded),**estimator.metadata()))
  return cache[key]
 truth_lookup=truth.set_index('case_id');rows=[];hidden_cells=[]
 for i,t in enumerate(collection['trials']):
  q=ish.loc[[t['sample_id']]].copy();q.loc[:,t['hidden']]=np.nan
  require(q[t['hidden']].isna().all().all(),'Máscara no aplicada')
  if t['experiment']=='real_patterns':require(set(q.columns[q.iloc[0].isna()])==set(t['hidden']),'Patrón no exacto')
  for met in t['hidden']:hidden_cells.append(dict(trial_id=t['trial_id'],sample_id=t['sample_id'],metabolite_id=met,originally_observed=pd.notna(ish.loc[t['sample_id'],met]),scored=met in t['score']))
  for sc in SCENARIOS:
   for method in METHODS:
    model=get_model(sc,method,t['excluded']);_,d=model.transform(q)
    for x in d[d.metabolite_id.isin(t['score'])].to_dict('records'):
     cid=t['trial_id']+'|'+x['metabolite_id'];tt=truth_lookup.loc[cid]
     x.update(fit_id=model.audit_fit_id,case_id=cid,trial_id=t['trial_id'],experiment=t['experiment'],scenario=sc,truth_R=float(tt.truth_R),truth_mM=float(tt.truth_mM),estimate_mM=x['estimate_R']*refs[x['metabolite_id']])
     x['absolute_error_mM']=abs(x['estimate_mM']-x['truth_mM']);x['absolute_error_R']=abs(x['estimate_R']-x['truth_R']);rows.append(x)
  if (i+1)%60==0:print(f'Validación: {i+1}/{len(collection["trials"])} ensayos',flush=True)
 pred=pd.DataFrame(rows).sort_values(['case_id','scenario','method']).reset_index(drop=True)
 require(len(pred)==len(truth)*6,'Intentos incompletos')
 for _,g in pred.groupby(['scenario','method']):require(set(g.case_id)==set(truth.case_id),'Celdas distintas entre métodos/escenarios')
 save(vdir/'predicciones.csv',pred);save(vdir/'celdas_ocultadas.csv',pd.DataFrame(hidden_cells))
 met,cov,common,delta=metrics(pred)
 for name,d in [('metricas',met),('cobertura_y_motivos',cov),('casos_comunes',common),('efecto_MMC2_pareado',delta)]:save(vdir/(name+'.csv'),d)
 # Gate explícito: las métricas completas se guardan antes de estimar faltantes reales.
 dump(vdir/'estado.json',dict(status='COMPLETE',truths='original Ishii only',trials=len(collection['trials']),scored_cases=len(truth),attempts=len(pred)))
 print('Validación completa. Generando candidatos reales.',flush=True)
 real_rows=[]
 missing_samples=ish.index[ish[targets].isna().any(axis=1)]
 for sample in missing_samples:
  group=metadata.loc[sample,'perturbation_group'];excluded=list(ish.index[metadata.loc[ish.index,'perturbation_group']==group])
  for sc in SCENARIOS:
   for method in METHODS:
    model=get_model(sc,method,excluded);_,d=model.transform(ish.loc[[sample]])
    for x in d[d.metabolite_id.isin(targets)].to_dict('records'):
     x.update(fit_id=model.audit_fit_id,scenario=sc,excluded_group=group,estimate_mM=x['estimate_R']*refs[x['metabolite_id']]);real_rows.append(x)
 real=pd.DataFrame(real_rows).sort_values(['scenario','method','sample_id','metabolite_id']);require(len(real)==114,'Huecos reales incorrectos')
 save(idir/'candidatos_y_soporte.csv',real)
 stats,members,matrices=statistics(mm,real,targets,vcfg['cv_mean_epsilon_mM'])
 save(sdir/'antes_despues.csv',stats);save(sdir/'casillas_comunes.csv',members)
 for (policy,sc,method),matrix in matrices.items():matrix.to_csv(idir/f'{policy}_{sc}_{method}_mM.csv',encoding='utf-8-sig',float_format='%.17g')
 dump(out/'ajustes_y_particiones.json',fit_records)
 require(before=={str(p.relative_to(ROOT)):sha(p) for p in input_paths},'Entradas modificadas')
 summary=dict(validation_cases_by_experiment=truth.groupby('experiment').size().to_dict(),trials=len(collection['trials']),evaluation_attempts=len(pred),real_attempts=len(real),fit_count=len(cache),checks=dict(same_cases_all_methods_scenarios=True,hidden_targets_absent_from_training=True,originals_unchanged=True,validation_precedes_real_imputation=True,no_log=True))
 versions={x:importlib.metadata.version(x) for x in ['numpy','pandas','scipy','scikit-learn']}
 dump(out/'manifest_estudio.json',dict(status='SUCCESS',summary=summary,input_sha256=before,mask_sha256=sha(masks/'masks.json'),truth_sha256=sha(masks/'verdades_ishii.csv'),python=platform.python_version(),versions=versions,method_config=cfg,validation_config=vcfg,output_sha256={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}))
 print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
