"""Un ajuste por entrenamiento y una consulta independiente por celda oculta."""
import warnings

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer, KNNImputer, SimpleImputer
from sklearn.linear_model import BayesianRidge
from sklearn.metrics.pairwise import nan_euclidean_distances
from sklearn.preprocessing import RobustScaler, StandardScaler

SEMILLA = 20261109
CONFIG = {
    "Mediana": [{}],
    "KNN": [{"n_neighbors": k, "weights": "uniform"} for k in (3, 5, 7)]
           + [{"n_neighbors": 3, "weights": "distance"}],
    "Iterative": [{"max_iter": 50, "tol": 0.001}, {"max_iter": 100, "tol": 0.0001}],
}


def elegibilidad(train, grupos, minimo=10):
    medidas = train.notna().astype(int)
    n = medidas.sum()
    ng = medidas.groupby(grupos.loc[train.index]).max().sum()
    basicas = (n >= minimo) & (ng >= 5) & (train.nunique() >= 2)
    pares = (medidas.T @ medidas).to_numpy(copy=True)
    np.fill_diagonal(pares, 0)  # Un metabolito no puede respaldarse a sí mismo.
    respaldo = (pares[:, basicas] >= 5).sum(axis=1)
    return pd.DataFrame({"n": n, "grupos": ng, "predictores_con_5_pares": respaldo,
                         "elegible": basicas & (respaldo >= 2)})


def diagnostico_extremos(train):
    """Sensibilidad descriptiva de la SD a retirar el punto más alejado de la mediana."""
    registros = []
    for nombre, serie in train.items():
        x = serie.dropna()
        q1, q3 = x.quantile([0.25, 0.75])
        extremo = (x - x.median()).abs().idxmax()
        sd = x.std(ddof=0)
        cambio = 1 - x.drop(extremo).std(ddof=0) / sd if sd > 0 else np.nan
        n_fuera = ((x < q1 - 1.5 * (q3 - q1)) | (x > q3 + 1.5 * (q3 - q1))).sum()
        registros.append((nombre, len(x), n_fuera, cambio))
    return pd.DataFrame(registros, columns=["metabolito", "n", "fuera_cercas_IQR", "caida_SD"])


def predecir(train, ishii, consultas, columnas, metodo, numero, referencia=None, robusto=False,
             covariable=None):
    """Solo train interviene en la mediana, las escalas, los vecinos y las regresiones."""
    salida = consultas.reset_index(drop=True).copy()
    salida["pred_mM"] = np.nan
    salida["motivo"] = "objetivo_no_elegible"
    salida["n_predictores"] = 0
    salida["n_donantes"] = np.nan
    if not len(columnas):
        return salida, {"rondas": 0, "aviso_convergencia": False}
    x = train[columnas].to_numpy(dtype=float, copy=True)
    q = ishii.loc[salida.muestra, columnas].to_numpy(dtype=float, copy=True)
    j = pd.Index(columnas).get_indexer(salida.metabolito)
    i = np.flatnonzero(j >= 0)
    q[i, j[i]] = np.nan  # Mantiene los demás faltantes reales de ese mismo cultivo.
    divisor = np.ones(len(columnas)) if referencia is None else referencia.loc[columnas].to_numpy()
    x, q = x / divisor, q / divisor
    observados_x, observados_q = np.isfinite(x), np.isfinite(q)
    contexto_activo = False
    if covariable is not None and metodo != "Mediana":
        cx = covariable.loc[train.index].to_numpy(dtype=float)
        cq = covariable.loc[salida.muestra].to_numpy(dtype=float)
        if not (np.isfinite(cx).all() and np.isfinite(cq).all()):
            raise ValueError("La covariable debe estar disponible en cada fila de esta comparación.")
        # Una columna constante tampoco debe contar en la distancia de KNN.
        contexto_activo = np.ptp(cx) > 0
        if contexto_activo:
            x, q = np.column_stack([x, cx]), np.column_stack([q, cq])


    params = CONFIG[metodo][numero]
    if metodo == "Mediana":
        modelo = SimpleImputer(strategy="median")
        escalador = None
    else:
        escalador = (RobustScaler() if robusto else StandardScaler()).fit(x)
        x, q = escalador.transform(x), escalador.transform(q)
        if metodo == "KNN":
            modelo = KNNImputer(**params)
        else:
            cero = escalador.transform(np.zeros((1, x.shape[1])))[0]
            modelo = IterativeImputer(estimator=BayesianRidge(), initial_strategy="median",
                                      imputation_order="ascending", skip_complete=False,
                                      random_state=SEMILLA, min_value=cero, **params)
    with warnings.catch_warnings(record=True) as avisos:
        warnings.simplefilter("always", ConvergenceWarning)
        modelo.fit(x)
        estimados = modelo.transform(q)
    if escalador is not None:
        estimados = escalador.inverse_transform(estimados)
    estimados = estimados[:, :len(columnas)] * divisor
    dist = nan_euclidean_distances(q, x) if metodo == "KNN" else None
    for fila in i:
        col = j[fila]
        n = observados_q[fila].sum()
        salida.at[fila, "n_predictores"] = n
        if metodo != "Mediana" and n < 2:
            salida.at[fila, "motivo"] = "menos_de_dos_predictores"
            continue
        if dist is not None:
            candidatos = np.isfinite(dist[fila]) & observados_x[:, col]
            if contexto_activo:
                solapan = (observados_x & observados_q[fila]).any(axis=1)
                if (candidatos & ~solapan).any():
                    salida.at[fila, "motivo"] = "posible_vecino_solo_por_covariable"
                    continue
            donantes = candidatos.sum()
            salida.at[fila, "n_donantes"] = donantes
            if donantes < 3:
                salida.at[fila, "motivo"] = "menos_de_tres_donantes"
                continue
        valor = estimados[fila, col]
        if not np.isfinite(valor) or valor < -1e-12:
            salida.at[fila, "motivo"] = "no_finita_o_negativa"
            continue
        salida.at[fila, "pred_mM"] = max(0.0, valor)
        salida.at[fila, "motivo"] = "mediana_general" if metodo == "Mediana" else "con_soporte"
    info = {"rondas": getattr(modelo, "n_iter_", 0), "aviso_convergencia":
            any(issubclass(a.category, ConvergenceWarning) for a in avisos)}
    return salida, info
