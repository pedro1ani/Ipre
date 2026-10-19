"""Tres imputadores sklearn con diagnóstico por celda. Sin KNN propio ni PCA."""
import warnings
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import SimpleImputer, KNNImputer, IterativeImputer
from sklearn.linear_model import BayesianRidge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import nan_euclidean_distances
from sklearn.exceptions import ConvergenceWarning

METHODS=('SimpleImputer','KNNImputer','IterativeImputer')

class AuditedImputer:
    def __init__(self, method, config):
        if method not in METHODS:
            raise ValueError('Método no permitido: '+method)
        self.method=method;self.config=config
        fixed={'estimator':'BayesianRidge','initial_strategy':'median','skip_complete':False,'sample_posterior':False,'imputation_order':'ascending'}
        if any(config['iterative'].get(k)!=v for k,v in fixed.items()):
            raise ValueError('Parámetro Iterative fuera del protocolo documentado')
        if config['simple']['strategy']!='median' or config['knn']['weights']!='uniform' or config['knn']['metric']!='nan_euclidean':
            raise ValueError('Parámetro fuera del protocolo documentado')
        for value in [config['knn']['n_neighbors'],config['iterative']['max_iter'],config['iterative']['regressor_max_iter'],config['support']['minimum_donors'],config['support']['minimum_observed_predictors']]:
            if type(value) is not int or value<1:raise ValueError('Conteo debe ser entero positivo')
        if config['iterative']['tol']<=0 or config['iterative']['regressor_tol']<=0:raise ValueError('Tolerancia debe ser positiva')

    def fit(self, train):
        if not train.index.is_unique or not train.columns.is_unique:
            raise ValueError('Identificadores duplicados')
        self.columns=list(train.columns);self.train=train.astype(float).copy()
        a=self.train.to_numpy()
        if np.isinf(a).any() or (a[np.isfinite(a)]<0).any():
            raise ValueError('Entradas deben ser ratios no negativos o NaN')
        self.active=[m for m in self.columns if self.train[m].notna().any()]
        self.empty=[m for m in self.columns if m not in self.active]
        self.warnings=[];self.model=None;self.scaler=None;self.converged=None;self.n_iter=None
        if not self.active:
            return self
        self.fit_data=self.train[self.active].to_numpy()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            if self.method=='SimpleImputer':
                self.model=SimpleImputer(strategy='median',add_indicator=False,keep_empty_features=False)
            else:
                self.scaler=StandardScaler()
                self.fit_data=self.scaler.fit_transform(self.fit_data)
                if self.method=='KNNImputer':
                    self.model=KNNImputer(n_neighbors=self.config['knn']['n_neighbors'],weights='uniform',metric='nan_euclidean',add_indicator=False,keep_empty_features=False)
                else:
                    c=self.config['iterative']
                    self.model=IterativeImputer(estimator=BayesianRidge(max_iter=c['regressor_max_iter'],tol=c['regressor_tol']),max_iter=c['max_iter'],tol=c['tol'],initial_strategy='median',imputation_order='ascending',n_nearest_features=None,skip_complete=False,sample_posterior=False,random_state=self.config['seed'],min_value=-self.scaler.mean_/self.scaler.scale_,max_value=np.inf,add_indicator=False,keep_empty_features=False)
            self.model.fit(self.fit_data)
        self.warnings=[dict(category=w.category.__name__,message=str(w.message)) for w in caught]
        if self.method=='IterativeImputer':
            self.n_iter=int(self.model.n_iter_)
            self.converged=(not any(issubclass(w.category,ConvergenceWarning) for w in caught)) if self.n_iter>0 else None
        return self

    def transform(self, query):
        if list(query.columns)!=self.columns:
            raise ValueError('Orden de columnas diferente al entrenamiento')
        x=query.astype(float).copy();a=x.to_numpy()
        if np.isinf(a).any() or (a[np.isfinite(a)]<0).any():
            raise ValueError('Consulta debe contener ratios no negativos o NaN')
        out=x.copy();diagnostics=[]
        if self.active:
            q=x[self.active].to_numpy();z=self.scaler.transform(q) if self.scaler is not None else q
            predicted=self.model.transform(z)
            raw=self.scaler.inverse_transform(predicted) if self.scaler is not None else predicted
            # Mantener exactamente las mediciones, sin roundtrip del escalador.
            for k,m in enumerate(self.active):
                missing=x[m].isna();out.loc[missing,m]=raw[missing.to_numpy(),k]
            distances=nan_euclidean_distances(z,self.fit_data) if self.method=='KNNImputer' else None
        for i,(sample,row) in enumerate(x.iterrows()):
            for met in self.columns:
                if pd.notna(row[met]):continue
                donors=self.train[met].dropna();preds=[m for m in self.active if m!=met]
                observed=[m for m in preds if pd.notna(row[m])]
                informative=[m for m in observed if self.train[m].dropna().nunique()>1]
                value=out.loc[sample,met]
                info=dict(sample_id=sample,metabolite_id=met,method=self.method,estimate_R=value,n_target_observed=len(donors),n_observed_predictors=len(observed),n_informative_predictors=len(informative),support='training_only',mechanism='',n_usable_neighbors=None,neighbor_ids='',min_shared_predictors=None,iterations=self.n_iter,converged=self.converged,roundoff_zero_corrected=False)
                if met in self.empty:
                    info.update(support='abstention',mechanism='no_target_values_in_training');out.loc[sample,met]=np.nan
                elif self.method=='SimpleImputer':
                    info['mechanism']='training_median'
                elif self.method=='KNNImputer':
                    k=self.active.index(met);donor_index=np.flatnonzero(~np.isnan(self.fit_data[:,k]));d=distances[i,donor_index]
                    if not np.isfinite(d).any():
                        info['mechanism']='automatic_training_mean_no_defined_distance';info['n_usable_neighbors']=0
                    else:
                        count=min(self.config['knn']['n_neighbors'],len(donor_index))
                        # Reproduce únicamente la selección para auditar; sklearn calcula la estimación.
                        chosen=np.argpartition(d,count-1)[:count];chosen=chosen[np.isfinite(d[chosen])];positions=donor_index[chosen]
                        shared=np.sum(~np.isnan(self.fit_data[positions]) & ~np.isnan(z[i]),axis=1)
                        info.update(n_usable_neighbors=len(positions),neighbor_ids=';'.join(str(self.train.index[j]) for j in positions),min_shared_predictors=int(shared.min()),mechanism='neighbor_mean')
                        adequate=len(positions)>=self.config['support']['minimum_donors'] and int(shared.min())>=self.config['support']['minimum_observed_predictors'] and len(informative)>=self.config['support']['minimum_observed_predictors']
                        info['support']='supported' if adequate else 'limited_support'
                        if len(informative)==0:
                            info.update(support='training_only',mechanism='neighbor_selection_without_informative_observation')
                    if donors.nunique()<=1:
                        info.update(support='training_only',mechanism='constant_target_in_training')
                else:
                    if self.n_iter==0:
                        info['mechanism']='initial_training_median_no_regression'
                    elif len(informative)==0:
                        info['mechanism']='iterative_population_only_no_informative_observation'
                    elif donors.nunique()<=1:
                        info['mechanism']='constant_target_in_training'
                    else:
                        info['mechanism']='conditional_iterative_model_with_internal_fills'
                        adequate=len(informative)>=self.config['support']['minimum_observed_predictors'] and len(donors)>=self.config['support']['minimum_donors'] and self.converged is True
                        info['support']='supported' if adequate else 'limited_support'
                if pd.notna(value) and value<0 and value>=-1e-12:
                    out.loc[sample,met]=0.0;info['roundoff_zero_corrected']=True
                elif pd.notna(value) and (not np.isfinite(value) or value<0):
                    out.loc[sample,met]=np.nan;info.update(support='abstention',mechanism='invalid_prediction')
                info['estimate_R']=out.loc[sample,met]
                diagnostics.append(info)
        return out,pd.DataFrame(diagnostics)

    def metadata(self):
        sequence=getattr(self.model,'imputation_sequence_',[])
        return dict(method=self.method,train_ids=list(self.train.index),active_features=self.active,empty_features=self.empty,scaler_mean=None if self.scaler is None else self.scaler.mean_.tolist(),scaler_scale=None if self.scaler is None else self.scaler.scale_.tolist(),iterations=self.n_iter,converged=self.converged,warnings=self.warnings,regressor_fits=len(sequence),regressor_max_iterations_used=max((int(t.estimator.n_iter_) for t in sequence),default=0),regressor_fits_at_limit=sum(t.estimator.n_iter_>=self.config['iterative']['regressor_max_iter'] for t in sequence))
