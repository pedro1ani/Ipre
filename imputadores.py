"""Tres métodos sklearn con el mismo ajuste, predicción y diagnóstico."""

import warnings

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer, KNNImputer, SimpleImputer
from sklearn.linear_model import BayesianRidge
from sklearn.metrics.pairwise import nan_euclidean_distances
from sklearn.preprocessing import StandardScaler

METHODS = ("SimpleImputer", "KNNImputer", "IterativeImputer")


class Imputador:
    """Conserva el entrenamiento para aplicar el mismo modelo a varias
    consultas."""

    def __init__(self, metodo, parametros):
        if metodo not in METHODS:
            raise ValueError(f"Método no permitido: {metodo}")
        for nombre in (
            "vecinos",
            "rondas_maximas",
            "minimo_donantes",
            "minimo_predictores",
        ):
            if type(parametros[nombre]) is not int or parametros[nombre] < 1:
                raise ValueError(f"{nombre} debe ser un entero positivo")
        if parametros["tolerancia"] <= 0:
            raise ValueError("La tolerancia debe ser positiva")
        self.metodo = metodo
        self.parametros = parametros

    def fit(self, entrenamiento):
        """Ajusta escalador y método exclusivamente con entrenamiento."""
        if (
            not entrenamiento.index.is_unique
            or not entrenamiento.columns.is_unique
        ):
            raise ValueError("Identificadores duplicados")
        self.entrenamiento = entrenamiento.astype(float).copy()
        valores = self.entrenamiento.to_numpy()
        if (
            np.isinf(valores).any()
            or (valores[np.isfinite(valores)] < 0).any()
        ):
            raise ValueError("Se esperan ratios no negativos o NaN")

        self.columnas = list(entrenamiento.columns)
        self.activas = list(entrenamiento.columns[entrenamiento.notna().any()])
        self.vacias = [
            nombre for nombre in self.columnas if nombre not in self.activas
        ]
        self.modelo = None
        self.escalador = None
        self.advertencias = []
        self.convergio = None
        self.rondas = None
        if not self.activas:
            return self

        self.datos_ajuste = self.entrenamiento[self.activas].to_numpy()
        with warnings.catch_warnings(record=True) as advertencias:
            warnings.simplefilter("always")
            if self.metodo == "SimpleImputer":
                self.modelo = SimpleImputer(strategy="median")
            else:
                self.escalador = StandardScaler()
                self.datos_ajuste = self.escalador.fit_transform(
                    self.datos_ajuste
                )
                if self.metodo == "KNNImputer":
                    self.modelo = KNNImputer(
                        n_neighbors=self.parametros["vecinos"],
                        weights="uniform",
                        metric="nan_euclidean",
                    )
                else:
                    # El límite R=0 se traslada a la escala estandarizada.
                    limite_inferior = (
                        -self.escalador.mean_ / self.escalador.scale_
                    )
                    self.modelo = IterativeImputer(
                        estimator=BayesianRidge(max_iter=300, tol=0.001),
                        max_iter=self.parametros["rondas_maximas"],
                        tol=self.parametros["tolerancia"],
                        initial_strategy="median",
                        imputation_order="ascending",
                        n_nearest_features=None,
                        skip_complete=False,
                        sample_posterior=False,
                        random_state=self.parametros["semilla"],
                        min_value=limite_inferior,
                    )
            self.modelo.fit(self.datos_ajuste)

        self.advertencias = [
            {
                "category": aviso.category.__name__,
                "message": str(aviso.message),
            }
            for aviso in advertencias
        ]
        if self.metodo == "IterativeImputer":
            self.rondas = int(self.modelo.n_iter_)
            if self.rondas > 0:
                self.convergio = not any(
                    issubclass(aviso.category, ConvergenceWarning)
                    for aviso in advertencias
                )
        return self

    def transform(self, consulta):
        """Completa ausencias y conserva cada medición original exactamente."""
        if list(consulta.columns) != self.columnas:
            raise ValueError("Orden de columnas diferente al entrenamiento")
        consulta = consulta.astype(float).copy()
        valores = consulta.to_numpy()
        if (
            np.isinf(valores).any()
            or (valores[np.isfinite(valores)] < 0).any()
        ):
            raise ValueError(
                "La consulta debe contener ratios no negativos o NaN"
            )
        completada = consulta.copy()
        distancias = None
        consulta_escalada = None
        if self.activas:
            consulta_escalada = consulta[self.activas].to_numpy()
            if self.escalador is not None:
                consulta_escalada = self.escalador.transform(consulta_escalada)
            estimaciones = self.modelo.transform(consulta_escalada)
            if self.escalador is not None:
                estimaciones = self.escalador.inverse_transform(estimaciones)
            for posicion, metabolito in enumerate(self.activas):
                faltantes = consulta[metabolito].isna()
                completada.loc[faltantes, metabolito] = estimaciones[
                    faltantes.to_numpy(), posicion
                ]
            if self.metodo == "KNNImputer":
                distancias = nan_euclidean_distances(
                    consulta_escalada, self.datos_ajuste
                )

        diagnosticos = []
        for posicion, (muestra, fila) in enumerate(consulta.iterrows()):
            for metabolito in consulta.columns[fila.isna()]:
                donantes = self.entrenamiento[metabolito].dropna()
                observados = [
                    nombre
                    for nombre in self.activas
                    if nombre != metabolito and pd.notna(fila[nombre])
                ]
                informativos = [
                    nombre
                    for nombre in observados
                    if self.entrenamiento[nombre].dropna().nunique() > 1
                ]
                valor = completada.loc[muestra, metabolito]
                diagnostico = {
                    "sample_id": muestra,
                    "metabolite_id": metabolito,
                    "method": self.metodo,
                    "estimate_R": valor,
                    "n_target_observed": len(donantes),
                    "n_observed_predictors": len(observados),
                    "n_informative_predictors": len(informativos),
                    "support": "training_only",
                    "mechanism": "",
                    "n_usable_neighbors": None,
                    "neighbor_ids": "",
                    "min_shared_predictors": None,
                    "iterations": self.rondas,
                    "converged": self.convergio,
                    "roundoff_zero_corrected": False,
                }
                if metabolito in self.vacias:
                    diagnostico.update(
                        support="abstention",
                        mechanism="no_target_values_in_training",
                    )
                elif self.metodo == "SimpleImputer":
                    diagnostico["mechanism"] = "training_median"
                elif self.metodo == "KNNImputer":
                    self._diagnosticar_vecinos(
                        diagnostico,
                        posicion,
                        metabolito,
                        distancias,
                        consulta_escalada,
                        donantes,
                        informativos,
                    )
                else:
                    self._diagnosticar_iterativo(
                        diagnostico, donantes, informativos
                    )

                if pd.notna(valor) and -1e-12 <= valor < 0:
                    completada.loc[muestra, metabolito] = 0.0
                    diagnostico["roundoff_zero_corrected"] = True
                elif pd.notna(valor) and (not np.isfinite(valor) or valor < 0):
                    completada.loc[muestra, metabolito] = np.nan
                    diagnostico.update(
                        support="abstention", mechanism="invalid_prediction"
                    )
                diagnostico["estimate_R"] = completada.loc[muestra, metabolito]
                diagnosticos.append(diagnostico)
        return completada, pd.DataFrame(diagnosticos)

    def _diagnosticar_vecinos(
        self,
        diagnostico,
        posicion,
        metabolito,
        distancias,
        consulta_escalada,
        donantes,
        informativos,
    ):
        columna = self.activas.index(metabolito)
        indices_donantes = np.flatnonzero(
            ~np.isnan(self.datos_ajuste[:, columna])
        )
        distancias_donantes = distancias[posicion, indices_donantes]
        if not np.isfinite(distancias_donantes).any():
            diagnostico.update(
                mechanism="automatic_training_mean_no_defined_distance",
                n_usable_neighbors=0,
            )
        else:
            # Replica la selección de sklearn solo para registrar los donantes.
            cantidad = min(self.parametros["vecinos"], len(indices_donantes))
            elegidos = np.argpartition(distancias_donantes, cantidad - 1)[
                :cantidad
            ]
            elegidos = elegidos[np.isfinite(distancias_donantes[elegidos])]
            vecinos = indices_donantes[elegidos]
            compartidos = np.sum(
                ~np.isnan(self.datos_ajuste[vecinos])
                & ~np.isnan(consulta_escalada[posicion]),
                axis=1,
            )
            soporte = (
                len(vecinos) >= self.parametros["minimo_donantes"]
                and compartidos.min() >= self.parametros["minimo_predictores"]
                and len(informativos) >= self.parametros["minimo_predictores"]
            )
            diagnostico.update(
                n_usable_neighbors=len(vecinos),
                neighbor_ids=";".join(
                    str(self.entrenamiento.index[i]) for i in vecinos
                ),
                min_shared_predictors=int(compartidos.min()),
                mechanism="neighbor_mean",
                support="supported" if soporte else "limited_support",
            )
            if not informativos:
                diagnostico.update(
                    support="training_only",
                    mechanism=(
                        "neighbor_selection_without_informative_o" "bservation"
                    ),
                )
        if donantes.nunique() <= 1:
            diagnostico.update(
                support="training_only",
                mechanism="constant_target_in_training",
            )

    def _diagnosticar_iterativo(self, diagnostico, donantes, informativos):
        if self.rondas == 0:
            diagnostico["mechanism"] = "initial_training_median_no_regression"
        elif not informativos:
            diagnostico["mechanism"] = (
                "iterative_population_only_no_informative_observation"
            )
        elif donantes.nunique() <= 1:
            diagnostico["mechanism"] = "constant_target_in_training"
        else:
            soporte = (
                len(informativos) >= self.parametros["minimo_predictores"]
                and len(donantes) >= self.parametros["minimo_donantes"]
                and self.convergio is True
            )
            diagnostico.update(
                mechanism="conditional_iterative_model_with_internal_fills",
                support="supported" if soporte else "limited_support",
            )

    def metadata(self):
        secuencia = getattr(self.modelo, "imputation_sequence_", [])
        iteraciones = [int(paso.estimator.n_iter_) for paso in secuencia]
        return {
            "method": self.metodo,
            "train_ids": list(self.entrenamiento.index),
            "active_features": self.activas,
            "empty_features": self.vacias,
            "scaler_mean": (
                None
                if self.escalador is None
                else self.escalador.mean_.tolist()
            ),
            "scaler_scale": (
                None
                if self.escalador is None
                else self.escalador.scale_.tolist()
            ),
            "iterations": self.rondas,
            "converged": self.convergio,
            "warnings": self.advertencias,
            "regressor_fits": len(secuencia),
            "regressor_max_iterations_used": max(iteraciones, default=0),
            "regressor_fits_at_limit": sum(
                numero >= 300 for numero in iteraciones
            ),
        }
