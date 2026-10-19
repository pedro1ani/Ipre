"""Imputa un único estudio por ejecución"""
from pathlib import Path
import argparse,hashlib,json,importlib.metadata,platform
import pandas as pd
import numpy as np
from imputadores import AuditedImputer,METHODS
ROOT=Path(__file__).resolve().parent

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def validate_input(data,samples):
    if not data.index.is_unique or not data.columns.is_unique:raise ValueError('IDs duplicados')
    if not data.index.isin(samples.index).all():raise ValueError('Muestra sin metadatos')

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--dataset',choices=['ishii','mmc2'],default='ishii')
 p.add_argument('--output',type=Path)
 args=p.parse_args();dataset=args.dataset
 out=(args.output or ROOT/'imputacion'/dataset).resolve()
 if out.exists():raise FileExistsError('Utilice una carpeta de salida nueva')
 config_path=ROOT/f'config_imputacion_{dataset}.json';cfg=json.loads(config_path.read_text());targets=cfg['targets']
 if cfg['dataset']!=dataset or tuple(cfg['methods'])!=METHODS:raise ValueError('Configuración no corresponde al estudio')
 matrix_path=ROOT/'normalizacion'/f'matriz_{dataset}_R.csv';data=pd.read_csv(matrix_path,index_col=0)
 samples_path=ROOT/'normalizacion/muestras_y_exclusiones.csv';samples=pd.read_csv(samples_path,index_col=0)
 validate_input(data,samples)
 if not targets or len(set(targets))!=len(targets) or not set(targets)<=set(data.columns):raise ValueError('Objetivos inválidos')
 ref_path=ROOT/'normalizacion/referencias_glucosa.csv';refs=pd.read_csv(ref_path,index_col=0).value_mM
 paths=[config_path,matrix_path,samples_path,ref_path,ROOT/'imputadores.py',Path(__file__)]
 hashes={str(f.relative_to(ROOT)):digest(f) for f in paths}
 results=[];fits=[];missing_samples=data.index[data[targets].isna().any(axis=1)]
 for group in sorted(set(samples.loc[missing_samples,'perturbation_group'])):
  excluded=data.index[samples.loc[data.index,'perturbation_group']==group];test_ids=missing_samples.intersection(excluded,sort=False)
  train=data.drop(index=excluded);query=data.loc[test_ids].copy()
  validate_input(train,samples)
  if set(train.index)&set(excluded):raise AssertionError('Filtración de grupo')
  for method in METHODS:
   estimator=AuditedImputer(method,cfg).fit(train);completed,diag=estimator.transform(query)
   observed=query.notna()
   if not np.array_equal(completed.to_numpy()[observed],query.to_numpy()[observed]):raise AssertionError('Medidas alteradas')
   diag=diag[diag.metabolite_id.isin(targets)].copy();diag['dataset']=dataset;diag['excluded_group']=group
   diag['estimate_mM']=diag.estimate_R*diag.metabolite_id.map(refs)
   diag['accepted_supported_mM']=diag.estimate_mM.where(diag.support=='supported')
   results.extend(diag.to_dict('records'));fits.append(dict(dataset=dataset,excluded_group=group,excluded_ids=list(excluded),query_ids=list(test_ids),**estimator.metadata()))
 if not results:raise ValueError('No hay faltantes en los objetivos seleccionados; no se requiere imputar')
 out.mkdir(parents=True,exist_ok=False)
 d=pd.DataFrame(results).sort_values(['method','sample_id','metabolite_id']).reset_index(drop=True)
 expected=int(data[targets].isna().sum().sum())*len(METHODS)
 if len(d)!=expected:raise AssertionError('Número de intentos incorrecto')
 d.to_csv(out/'candidatos_y_soporte.csv',index=False,encoding='utf-8-sig')
 summary=d.groupby(['dataset','method','support']).size().rename('cells').reset_index();summary.to_csv(out/'resumen_soporte.csv',index=False,encoding='utf-8-sig')
 (out/'ajustes_y_particiones.json').write_text(json.dumps(fits,indent=2,ensure_ascii=False,default=lambda x:int(x) if isinstance(x,np.integer) else x)+'\n')
 for method in METHODS:
  selected=d[d.method==method];matrix=data[targets].mul(refs[targets],axis=1)
  for row in selected.itertuples():
   if row.support=='supported':matrix.loc[row.sample_id,row.metabolite_id]=row.estimate_mM
  matrix.to_csv(out/f'objetivos_soportados_{method}_mM.csv',encoding='utf-8-sig')
 if hashes!={str(f.relative_to(ROOT)):digest(f) for f in paths}:raise AssertionError('Entradas modificadas')
 manifest=dict(stage='single_dataset_candidates_not_accuracy_validated',dataset=dataset,input_sha256=hashes,config=cfg,python=platform.python_version(),versions={x:importlib.metadata.version(x) for x in ['numpy','pandas','scipy','scikit-learn']},attempts=len(d),originals_unchanged=True,group_disjoint=True,training_studies=sorted(samples.loc[data.index,'dataset_id'].unique().tolist()),output_sha256={f.name:digest(f) for f in sorted(out.iterdir())})
 (out/'manifest_imputacion.json').write_text(json.dumps(manifest,indent=2)+'\n');print(summary.to_string(index=False))
if __name__=='__main__':main()
