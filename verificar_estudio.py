"""Verifica separación, máscaras, unidades, métricas y conservación sin reentrenar."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent

def main():
 p=argparse.ArgumentParser();p.add_argument('--results',type=Path,default=ROOT/'resultados_estudio');p.add_argument('--masks',type=Path,default=ROOT/'validacion_comun');args=p.parse_args();out=args.results
 manifest=json.loads((out/'manifest_estudio.json').read_text());collection=json.loads((args.masks/'masks.json').read_text());fits=json.loads((out/'ajustes_y_particiones.json').read_text());trials={t['trial_id']:t for t in collection['trials']}
 for name,h in manifest['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
 for name,h in manifest['output_sha256'].items():assert hashlib.sha256((out/name).read_bytes()).hexdigest()==h,name
 assert hashlib.sha256((args.masks/'masks.json').read_bytes()).hexdigest()==manifest['mask_sha256']
 ish=pd.read_csv(ROOT/'normalizacion/matriz_ishii_R.csv',index_col=0);mm=pd.read_csv(ROOT/'normalizacion/matriz_ishii_mM.csv',index_col=0);ger=pd.read_csv(ROOT/'normalizacion/matriz_mmc2_R.csv',index_col=0)
 references=pd.read_csv(ROOT/'normalizacion/referencias_glucosa.csv',index_col=0).value_mM
 for fit in fits:
  train=ish.drop(index=fit['excluded'])
  if fit['scenario']=='ishii_mmc2':train=pd.concat([train,ger])
  assert list(train.index)==fit['train_ids'];assert not set(fit['excluded'])&set(train.index)
  if fit['scaler_mean'] is not None:np.testing.assert_allclose(fit['scaler_mean'],train[fit['active_features']].mean().to_numpy(),rtol=1e-12,atol=1e-14)
 pred=pd.read_csv(out/'validacion/predicciones.csv');truth=pd.read_csv(args.masks/'verdades_ishii.csv')
 for _,g in pred.groupby(['scenario','method']):assert set(g.case_id)==set(truth.case_id)
 for row in pred.itertuples():
  t=trials[row.trial_id];fit=fits[row.fit_id]
  assert row.sample_id in t['excluded'];assert row.metabolite_id in t['hidden'];assert row.metabolite_id in t['score']
  assert fit['excluded']==sorted(t['excluded']);assert fit['scenario']==row.scenario and fit['method']==row.method
  assert np.isclose(mm.loc[row.sample_id,row.metabolite_id],row.truth_mM,rtol=1e-12,atol=1e-14)
  assert np.isclose(ish.loc[row.sample_id,row.metabolite_id],row.truth_R,rtol=1e-12,atol=1e-14)
  assert np.isclose(row.absolute_error_mM,abs(row.estimate_mM-row.truth_mM),rtol=1e-12,atol=1e-14,equal_nan=True)
  assert np.isclose(row.absolute_error_R,abs(row.estimate_R-row.truth_R),rtol=1e-12,atol=1e-14,equal_nan=True)
  assert np.isclose(row.estimate_mM,row.estimate_R*references[row.metabolite_id],rtol=1e-12,atol=1e-14,equal_nan=True)
 for path in (out/'imputacion').glob('*_mM.csv'):
  d=pd.read_csv(path,index_col=0);observed=mm.notna();np.testing.assert_allclose(d.to_numpy()[observed],mm.to_numpy()[observed],rtol=1e-12,atol=1e-14)
 stats=pd.read_csv(out/'estadisticas/antes_despues.csv');assert ((stats.observed+stats.imputed+stats.pending)==30).all()
 memberships=pd.read_csv(out/'estadisticas/casillas_comunes.csv')
 for (_,sc),g in memberships.groupby(['policy','scenario']):
  sets=[set(zip(part.loc[part.common_methods,'sample_id'],part.loc[part.common_methods,'metabolite_id'])) for _,part in g.groupby('method')]
  assert all(s==sets[0] for s in sets)
 for _,g in memberships.groupby('policy'):
  sets=[set(zip(part.loc[part.common_all_six,'sample_id'],part.loc[part.common_all_six,'metabolite_id'])) for _,part in g.groupby(['scenario','method'])];assert all(s==sets[0] for s in sets)
 print('OK: hashes, train exclusivo, escalado, 2610 comparaciones, unidades, originales y casillas comunes.')
if __name__=='__main__':main()
