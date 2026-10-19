"""Pruebas de las máscaras entregadas y de ausencia de filtración."""
import json,sys,unittest
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from imputadores import AuditedImputer,METHODS

class ValidationTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.collection=json.loads((ROOT/'resultados/ishii/validacion/mascaras.json').read_text())
  cls.ish=pd.read_csv(ROOT/'normalizacion/matriz_ishii_R.csv',index_col=0)
  cls.meta=pd.read_csv(ROOT/'normalizacion/muestras_y_exclusiones.csv',index_col=0)
  cls.series=pd.read_csv(ROOT/'resultados/ishii/validacion/series_analiticas.csv',index_col=0)
  cls.cfg=json.loads((ROOT/'config_imputacion_ishii.json').read_text())
 def test_masks_score_only_originals_and_patterns_exact(self):
  for t in self.collection['trials']:
   self.assertTrue(self.ish.loc[t['sample_id'],t['score']].notna().all())
   q=self.ish.loc[[t['sample_id']]].copy();q.loc[:,t['hidden']]=np.nan
   if t['experiment']=='real_patterns':self.assertEqual(set(q.columns[q.iloc[0].isna()]),set(t['hidden']))
 def test_exclusion_complete(self):
  for t in self.collection['trials']:
   sample=t['sample_id'];group=self.meta.loc[sample,'perturbation_group']
   expected=set(self.ish.index[self.meta.loc[self.ish.index,'perturbation_group']==group])
   if t['experiment']=='analytical_series':
    family=self.series.loc[sample,'series_family'];expected|=set(self.series.index[self.series.series_family==family])
   self.assertEqual(expected,set(t['excluded']))
 def test_hidden_sentinel_does_not_change_predictions(self):
  # Un caso por experimento, todos los métodos. Valores ocultos no entran al transform.
  for experiment in ['perturbation','real_patterns','analytical_series']:
   t=next(t for t in self.collection['trials'] if t['experiment']==experiment)
   a=self.ish.loc[[t['sample_id']]].copy();b=a.copy();b.loc[:,t['hidden']]=1e12
   a.loc[:,t['hidden']]=np.nan;b.loc[:,t['hidden']]=np.nan
   train=self.ish.drop(index=t['excluded'])
   for method in METHODS:
    model=AuditedImputer(method,self.cfg).fit(train)
    pa,_=model.transform(a);pb,_=model.transform(b);pd.testing.assert_frame_equal(pa,pb)

if __name__=='__main__':unittest.main()
