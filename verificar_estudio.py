"""Verifica resultados de una ejecución independiente sin reentrenar
modelos."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parent


def main():
    argumentos = argparse.ArgumentParser()
    argumentos.add_argument(
        "--results", type=Path, default=RAIZ / "resultados/ishii"
    )
    salida = argumentos.parse_args().results
    registro_ejecucion = json.loads((salida / "manifest.json").read_text())
    estudio = registro_ejecucion["dataset"]
    for nombre, esperado in registro_ejecucion["input_sha256"].items():
        if (
            hashlib.sha256((RAIZ / nombre).read_bytes()).hexdigest()
            != esperado
        ):
            raise ValueError(nombre)
    for nombre, esperado in registro_ejecucion["output_sha256"].items():
        if (
            hashlib.sha256((salida / nombre).read_bytes()).hexdigest()
            != esperado
        ):
            raise ValueError(nombre)
    ratios = pd.read_csv(
        RAIZ / "normalizacion" / f"matriz_{estudio}_R.csv", index_col=0
    )
    concentraciones = pd.read_csv(
        RAIZ / "normalizacion" / f"matriz_{estudio}_mM.csv", index_col=0
    )
    metadatos = pd.read_csv(
        RAIZ / "normalizacion/muestras_y_exclusiones.csv", index_col=0
    )
    referencias = pd.read_csv(
        RAIZ / "normalizacion/referencias_glucosa.csv", index_col=0
    ).value_mM
    ajustes = json.loads((salida / "ajustes.json").read_text())
    for ajuste in ajustes:
        entrenamiento = ratios.drop(index=ajuste["excluded"])
        if list(entrenamiento.index) != ajuste["train_ids"]:
            raise ValueError(
                ("Las filas de entrenamiento no coinciden con el registro.")
            )
        if set(ajuste["excluded"]) & set(entrenamiento.index):
            raise ValueError(
                ("El entrenamiento incluye una muestra excluida.")
            )
        if ajuste["scaler_mean"] is not None:
            np.testing.assert_allclose(
                ajuste["scaler_mean"],
                entrenamiento[ajuste["active_features"]].mean(),
                rtol=1e-12,
                atol=1e-14,
            )
    mascaras = json.loads((salida / "validacion/mascaras.json").read_text())
    ensayos = {ensayo["trial_id"]: ensayo for ensayo in mascaras["trials"]}
    predicciones = pd.read_csv(salida / "validacion/predicciones.csv")
    verdades = pd.read_csv(salida / "validacion/verdades.csv")
    for _, grupo in predicciones.groupby("method"):
        if set(grupo.case_id) != set(verdades.case_id):
            raise ValueError(("Los métodos no evaluaron los mismos casos."))
    if estudio == "mmc2":
        if "MMC2_GLUCOSE" in set(predicciones.sample_id):
            raise ValueError(
                ("Glucosa MMC2 no debe usarse como verdad de evaluación.")
            )
        if "analytical_series" in set(predicciones.experiment):
            raise ValueError(
                ("MMC2 carece de series analíticas homologables.")
            )
    for registro in predicciones.itertuples():
        ensayo, ajuste = (ensayos[registro.trial_id], ajustes[registro.fit_id])
        if registro.sample_id not in ensayo["excluded"]:
            raise ValueError(
                ("La consulta no está excluida del entrenamiento.")
            )
        if not (
            registro.metabolite_id in ensayo["hidden"]
            and registro.metabolite_id in ensayo["score"]
        ):
            raise ValueError(
                ("La celda puntuada no corresponde a la máscara.")
            )
        if not (
            ajuste["excluded"] == sorted(ensayo["excluded"])
            and ajuste["method"] == registro.method
        ):
            raise ValueError(
                ("El ajuste no corresponde al ensayo o al método.")
            )
        if not np.isclose(
            registro.truth_mM,
            concentraciones.loc[registro.sample_id, registro.metabolite_id],
        ):
            raise ValueError(
                ("La verdad mM no coincide con la medición original.")
            )
        if not np.isclose(
            registro.truth_R,
            ratios.loc[registro.sample_id, registro.metabolite_id],
        ):
            raise ValueError(
                ("La verdad R no coincide con el ratio original.")
            )
        if not np.isclose(
            registro.absolute_error_mM,
            abs(registro.estimate_mM - registro.truth_mM),
            equal_nan=True,
        ):
            raise ValueError(
                ("El error absoluto no coincide con la predicción.")
            )
        if not np.isclose(
            registro.estimate_mM,
            registro.estimate_R * referencias[registro.metabolite_id],
            equal_nan=True,
        ):
            raise ValueError(("La estimación en mM no coincide con su ratio."))
    candidatos = pd.read_csv(salida / "imputacion/candidatos_y_soporte.csv")
    for registro in candidatos.itertuples():
        if not pd.isna(ratios.loc[registro.sample_id, registro.metabolite_id]):
            raise ValueError(
                ("Se intentó imputar una celda originalmente medida.")
            )
        if registro.metabolite_id not in registro_ejecucion["targets"]:
            raise ValueError(
                ("Se intentó completar un metabolito que no es objetivo.")
            )
        excluidas = sorted(
            ratios.index[
                metadatos.loc[ratios.index, "perturbation_group"]
                == metadatos.loc[registro.sample_id, "perturbation_group"]
            ]
        )
        if ajustes[registro.fit_id]["excluded"] != excluidas:
            raise ValueError(
                ("La exclusión del candidato no abarca todo su grupo.")
            )
    for archivo in (salida / "imputacion").glob("*_mM.csv"):
        matriz = pd.read_csv(archivo, index_col=0)
        if not (
            matriz.index.equals(concentraciones.index)
            and matriz.columns.equals(concentraciones.columns)
        ):
            raise ValueError(
                ("La matriz completada cambió el orden de filas o columnas.")
            )
        observados = concentraciones.notna().to_numpy()
        np.testing.assert_allclose(
            matriz.to_numpy()[observados],
            concentraciones.to_numpy()[observados],
            rtol=1e-12,
            atol=1e-14,
        )
    if not candidatos.empty:
        estadisticas = pd.read_csv(
            salida / "imputacion/estadisticas_antes_despues.csv"
        )
        if not (
            estadisticas.observed + estadisticas.imputed + estadisticas.pending
            == len(ratios)
        ).all():
            raise ValueError(("Los conteos no suman el tamaño de la cohorte."))
        comunes = pd.read_csv(salida / "imputacion/casillas_comunes.csv")
        for _, tabla in comunes.groupby("policy"):
            conjuntos = [
                set(
                    zip(
                        grupo.loc[grupo.common_methods, "sample_id"],
                        grupo.loc[grupo.common_methods, "metabolite_id"],
                    )
                )
                for _, grupo in tabla.groupby("method")
            ]
            if not all((muestra == conjuntos[0] for muestra in conjuntos)):
                raise ValueError(
                    ("Las casillas comunes difieren entre métodos.")
                )
    print(
        f"OK {estudio}: base única y entrenamiento exclusivo; "
        "errores y observados verificados."
    )


if __name__ == "__main__":
    main()
