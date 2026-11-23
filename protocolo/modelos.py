"""Tres imputadores de scikit-learn. Nunca se utilizan estimaciones como observaciones."""
import warnings

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer, KNNImputer, SimpleImputer
from sklearn.preprocessing import StandardScaler

from .datos import SEMILLA

METODOS = ("Mediana", "KNN3", "KNN5", "Iterativo")
RONDAS = 5
TOL = .001


def matriz_consultas(ish, q, cols, auxiliar=None):
    datos = ish.reindex(index=q.cultivo, columns=cols).to_numpy(copy=True)
    datos[:,0] = np.nan
    for j,c in enumerate(cols[1:],1):
        oculta = [c in mask.split("|") for mask in q.ocultas]
        datos[oculta,j] = np.nan
    if auxiliar is not None:
        datos = np.column_stack([datos, auxiliar.reindex(q.cultivo).to_numpy()])
    return datos


def iterador(z, scaler, posterior=False, semilla=SEMILLA):
    modelo = IterativeImputer(max_iter=RONDAS, tol=TOL, sample_posterior=posterior,
                              random_state=semilla, min_value=-scaler.mean_ / scaler.scale_,
                              skip_complete=False, initial_strategy="median")
    with warnings.catch_warnings(record=True) as adv:
        warnings.simplefilter("always", ConvergenceWarning)
        modelo.fit(z)
    alerta = any(issubclass(w.category, ConvergenceWarning) for w in adv)
    return modelo, alerta


def estimar(train, consulta, n_met, metodos=METODOS):
    """KNN solo ve donantes con objetivo medido y al menos dos metabolitos en común."""
    n = len(consulta)
    cumple = np.isfinite(consulta[:,1:n_met])
    visibles = cumple.sum(axis=1)
    medidas = np.isfinite(train[:,1:n_met])
    donantes = (cumple.astype(int) @ medidas.T.astype(int) >= 2) & np.isfinite(train[:,0])
    n_donantes = donantes.sum(axis=1)
    valido_aux = np.isfinite(consulta[:,-1]) if train.shape[1] > n_met else np.ones(n,dtype=bool)
    out = {}
    scaler = StandardScaler().fit(train)
    z, v = scaler.transform(train), scaler.transform(consulta)
    for met in metodos:
        pred = np.full(n,np.nan)
        motivo = np.full(n,"sin_auxiliar", dtype=object)
        it, alerta = 0, False
        if met == "Mediana":
            pred[valido_aux] = SimpleImputer(strategy="median").fit(train[:,:1]).statistics_[0]
            motivo[valido_aux] = "general_entrenamiento"
        elif met.startswith("KNN"):
            k = int(met[-1])
            motivo[valido_aux] = "menos_de_2_acompanantes"
            motivo[valido_aux & (visibles>=2)] = "menos_de_3_donantes"
            for mascara in np.unique(donantes[valido_aux & (n_donantes>=3)],axis=0):
                filas = valido_aux & (donantes==mascara).all(axis=1)
                knn = KNNImputer(n_neighbors=k, keep_empty_features=True).fit(z[mascara])
                pred[filas] = knn.transform(v[filas])[:,0]*scaler.scale_[0]+scaler.mean_[0]
                motivo[filas] = "donantes_observados"
        else:
            valido = valido_aux & (visibles>=2)
            motivo[valido_aux] = "menos_de_2_acompanantes"
            if valido.any():
                modelo, alerta = iterador(z,scaler)
                pred[valido] = modelo.transform(v[valido])[:,0]*scaler.scale_[0]+scaler.mean_[0]
                it = modelo.n_iter_
                motivo[valido] = "modelo_con_acompanantes"
        mala = np.isfinite(pred) & (pred < -1e-12)
        pred[mala], motivo[mala] = np.nan, "resultado_negativo"
        pred[np.isfinite(pred)] = np.maximum(0,pred[np.isfinite(pred)])
        out[met] = (pred, motivo, visibles, np.minimum(n_donantes,int(met[-1])) if met.startswith("KNN")
                    else np.zeros(n,dtype=int), it, alerta)
    return out
