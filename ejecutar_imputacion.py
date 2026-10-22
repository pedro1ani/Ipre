"""Evalúa con máscaras y luego completa candidatos de un solo estudio.

Uso: python ejecutar_imputacion.py --dataset ishii --output resultados/ishii
MMC2 se procesa por separado cambiando únicamente --dataset y --output.
"""

import argparse
import importlib.metadata
import json
import platform
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

from apoyo import describir, guardar_csv, guardar_json, huella
from imputadores import METHODS, Imputador

RAIZ = Path(__file__).resolve().parent
POLITICAS = ("available", "supported_plus_baseline")


def validar_entrada(ratios, muestras):
    if not (ratios.index.is_unique and ratios.columns.is_unique):
        raise ValueError("IDs duplicados")
    if not ratios.index.isin(muestras.index).all():
        raise ValueError("Muestra sin metadatos")


def aceptadas(tabla, politica):
    finitos = np.isfinite(tabla.estimate_R)
    if politica == "available":
        return finitos
    return finitos & (
        (tabla.method == "SimpleImputer") | (tabla.support == "supported")
    )


def crear_mascaras(
    ratios,
    concentraciones,
    metadatos,
    series_analiticas,
    objetivos,
    fingerprint,
    dataset="ishii",
):
    """Fija ensayos compartidos; las verdades se guardan fuera de la
    consulta."""
    ensayos = []
    patrones = []
    grupos = metadatos.loc[ratios.index, "perturbation_group"]
    exclusiones = {}
    for muestra in ratios.index:
        exclusiones[muestra] = sorted(ratios.index[grupos == grupos[muestra]])

    for muestra in ratios.index:
        # Evaluar glucosa revelaría la respuesta mediante el denominador fijo.
        if dataset == "mmc2" and muestra == "MMC2_GLUCOSE":
            continue
        por_prueba = {"perturbation": exclusiones[muestra]}
        if series_analiticas is not None:
            familias = series_analiticas["series_family"]
            misma_serie = familias == familias[muestra]
            mismo_grupo = grupos == grupos[muestra]
            por_prueba["analytical_series"] = sorted(
                ratios.index[misma_serie | mismo_grupo]
            )
        for metabolito in objetivos:
            if pd.isna(ratios.loc[muestra, metabolito]):
                continue
            for prueba, excluidas in por_prueba.items():
                ensayos.append(
                    {
                        "trial_id": f"{prueba}:{muestra}:{metabolito}",
                        "experiment": prueba,
                        "sample_id": muestra,
                        "hidden": [metabolito],
                        "score": [metabolito],
                        "excluded": excluidas,
                        "pattern_id": "",
                    }
                )

    patrones_reales = {}
    incompletas = ratios.index[ratios[objetivos].isna().any(axis=1)]
    for muestra in incompletas:
        faltantes = tuple(ratios.columns[ratios.loc[muestra].isna()])
        patrones_reales.setdefault(faltantes, []).append(muestra)
    for numero, (faltantes, origenes) in enumerate(
        sorted(patrones_reales.items()), start=1
    ):
        patron = f"P{numero:02d}"
        ocultas = list(faltantes)
        evaluadas = [nombre for nombre in objetivos if nombre in ocultas]
        visibles = [nombre for nombre in ratios if nombre not in ocultas]
        # Exigir visibles medidas reproduce exactamente el patrón de ausencia.
        elegibles = ratios.index[
            ratios[evaluadas].notna().all(axis=1)
            & ratios[visibles].notna().all(axis=1)
        ]
        patrones.append(
            {
                "pattern_id": patron,
                "source_samples": origenes,
                "hidden": ocultas,
                "score": evaluadas,
                "eligible_samples": list(elegibles),
            }
        )
        for muestra in elegibles:
            if dataset == "mmc2" and muestra == "MMC2_GLUCOSE":
                continue
            ensayos.append(
                {
                    "trial_id": f"real_patterns:{muestra}:{patron}",
                    "experiment": "real_patterns",
                    "sample_id": muestra,
                    "hidden": ocultas,
                    "score": evaluadas,
                    "excluded": exclusiones[muestra],
                    "pattern_id": patron,
                }
            )

    ensayos.sort(key=lambda ensayo: ensayo["trial_id"])
    if len({ensayo["trial_id"] for ensayo in ensayos}) != len(ensayos):
        raise ValueError("Ensayos duplicados")
    verdades = []
    for ensayo in ensayos:
        muestra = ensayo["sample_id"]
        for metabolito in ensayo["score"]:
            valor = ratios.loc[muestra, metabolito]
            if pd.isna(valor):
                raise ValueError("La verdad de evaluación debe estar medida")
            verdades.append(
                {
                    "case_id": ensayo["trial_id"] + "|" + metabolito,
                    "trial_id": ensayo["trial_id"],
                    "sample_id": muestra,
                    "metabolite_id": metabolito,
                    "experiment": ensayo["experiment"],
                    "truth_R": float(valor),
                    "truth_mM": float(
                        concentraciones.loc[muestra, metabolito]
                    ),
                }
            )
    verdades.sort(key=lambda registro: registro["case_id"])
    coleccion = {
        "schema": 1,
        "fingerprint": fingerprint,
        "feature_order": list(ratios.columns),
        "sample_order": list(ratios.index),
        "trials": ensayos,
        "patterns": patrones,
    }
    return coleccion, pd.DataFrame(verdades)


def evaluar_predicciones(predicciones):
    """Compara errores propios y sobre las mismas celdas para los tres
    métodos."""
    filas, motivos, pertenencias = ([], [], [])
    for politica in POLITICAS:
        casos_por_metodo = {
            metodo: set(
                tabla_metodo.loc[aceptadas(tabla_metodo, politica), "case_id"]
            )
            for metodo, tabla_metodo in predicciones.groupby("method")
        }
        comunes = set.intersection(
            *(casos_por_metodo[metodo] for metodo in METHODS)
        )
        for metodo, tabla_metodo in predicciones.groupby("method"):
            for registro in tabla_metodo.itertuples():
                respondio = registro.case_id in casos_por_metodo[metodo]
                pertenencias.append(
                    dict(
                        policy=politica,
                        method=metodo,
                        case_id=registro.case_id,
                        responded=respondio,
                        common_methods=registro.case_id in comunes,
                    )
                )
                motivo = (
                    "responded"
                    if respondio
                    else f"{registro.support}: {registro.mechanism}"
                )
                if (
                    not respondio
                    and registro.method == "IterativeImputer"
                    and (not registro.converged)
                ):
                    motivo += "; not_converged"
                motivos.append(
                    dict(
                        policy=politica,
                        method=metodo,
                        experiment=registro.experiment,
                        metabolite_id=registro.metabolite_id,
                        reason=motivo,
                    )
                )
            for alcance, seleccionados in [
                ("own", casos_por_metodo[metodo]),
                ("common_methods", comunes),
            ]:
                for (prueba, metabolito), grupo in tabla_metodo.groupby(
                    ["experiment", "metabolite_id"]
                ):
                    respuestas = grupo[grupo.case_id.isin(seleccionados)]
                    filas.append(
                        dict(
                            policy=politica,
                            scope=alcance,
                            method=metodo,
                            experiment=prueba,
                            metabolite_id=metabolito,
                            attempts=len(grupo),
                            predictions=len(respuestas),
                            coverage_pct=100 * len(respuestas) / len(grupo),
                            MAE_mM=respuestas.absolute_error_mM.mean(),
                            MAE_R=respuestas.absolute_error_R.mean(),
                        )
                    )
    motivos = (
        pd.DataFrame(motivos)
        .value_counts(sort=False)
        .rename("cells")
        .reset_index()
    )
    return (pd.DataFrame(filas), motivos, pd.DataFrame(pertenencias))


def resumir_candidatos(concentraciones, candidatos, carpeta):
    """Conserva los observados; completa solo objetivos candidatos y calcula
    estadísticas."""
    estadisticas, pertenencias = ([], [])
    for politica in POLITICAS:
        casos_por_metodo = {
            metodo: set(
                zip(
                    tabla_metodo.loc[
                        aceptadas(tabla_metodo, politica), "sample_id"
                    ],
                    tabla_metodo.loc[
                        aceptadas(tabla_metodo, politica), "metabolite_id"
                    ],
                )
            )
            for metodo, tabla_metodo in candidatos.groupby("method")
        }
        comunes = set.intersection(
            *(casos_por_metodo[metodo] for metodo in METHODS)
        )
        for metodo in METHODS:
            tabla_metodo = candidatos[candidatos.method == metodo]
            for registro in tabla_metodo.itertuples():
                clave = (registro.sample_id, registro.metabolite_id)
                pertenencias.append(
                    dict(
                        policy=politica,
                        method=metodo,
                        sample_id=registro.sample_id,
                        metabolite_id=registro.metabolite_id,
                        selected=clave in casos_por_metodo[metodo],
                        common_methods=clave in comunes,
                    )
                )
            for alcance, seleccionados in [
                ("own", casos_por_metodo[metodo]),
                ("common_methods", comunes),
            ]:
                completada = concentraciones.copy()
                for registro in tabla_metodo.itertuples():
                    if (
                        registro.sample_id,
                        registro.metabolite_id,
                    ) in seleccionados:
                        if not pd.isna(
                            completada.loc[
                                registro.sample_id, registro.metabolite_id
                            ]
                        ):
                            raise ValueError("Observado modificado")
                        completada.loc[
                            registro.sample_id, registro.metabolite_id
                        ] = registro.estimate_mM
                if alcance == "own":
                    completada.to_csv(
                        carpeta / f"{politica}_{metodo}_mM.csv",
                        encoding="utf-8-sig",
                        float_format="%.17g",
                    )
                for etapa, matriz in [
                    ("before", concentraciones),
                    ("after", completada),
                ]:
                    for metabolito in matriz:
                        valores = matriz[metabolito].dropna()
                        observados = int(concentraciones[metabolito].count())
                        estadisticas.append(
                            dict(
                                policy=politica,
                                scope=alcance,
                                method=metodo,
                                stage=etapa,
                                metabolite_id=metabolito,
                                total_samples=len(matriz),
                                observed=observados,
                                imputed=len(valores) - observados,
                                pending=len(matriz) - len(valores),
                                **describir(valores),
                            )
                        )
    return (pd.DataFrame(estadisticas), pd.DataFrame(pertenencias))


def main():
    argumentos = argparse.ArgumentParser(description=__doc__)
    argumentos.add_argument(
        "--dataset", choices=["ishii", "mmc2"], default="ishii"
    )
    argumentos.add_argument("--output", type=Path)
    opciones = argumentos.parse_args()
    dataset = opciones.dataset
    salida = (opciones.output or RAIZ / "resultados" / dataset).resolve()
    if salida.exists():
        raise ValueError("Utilice una carpeta de salida nueva")
    # Preparación: una base por ejecución y una configuración común.
    archivo_parametros = RAIZ / "config_imputacion.json"
    archivo_seleccion = RAIZ / "config_normalizacion.json"
    parametros = json.loads(archivo_parametros.read_text())
    objetivos = json.loads(archivo_seleccion.read_text())["targets"]
    archivo_ratios = RAIZ / "normalizacion" / f"matriz_{dataset}_R.csv"
    archivo_concentraciones = (
        RAIZ / "normalizacion" / f"matriz_{dataset}_mM.csv"
    )
    archivo_muestras = RAIZ / "normalizacion/muestras_y_exclusiones.csv"
    archivo_referencias = RAIZ / "normalizacion/referencias_glucosa.csv"
    ratios = pd.read_csv(archivo_ratios, index_col=0)
    concentraciones = pd.read_csv(archivo_concentraciones, index_col=0)
    muestras = pd.read_csv(archivo_muestras, index_col=0)
    referencias = pd.read_csv(archivo_referencias, index_col=0).value_mM
    validar_entrada(ratios, muestras)
    if not (
        bool(objetivos)
        and len(set(objetivos)) == len(objetivos)
        and (set(objetivos) <= set(ratios.columns))
    ):
        raise ValueError("Objetivos inválidos")
    if not (
        ratios.index.equals(concentraciones.index)
        and ratios.columns.equals(concentraciones.columns)
    ):
        raise ValueError("Orden inconsistente")
    if not np.allclose(
        ratios.mul(referencias).to_numpy(),
        concentraciones.to_numpy(),
        equal_nan=True,
    ):
        raise ValueError("mM y ratios inconsistentes")
    entradas = [
        archivo_parametros,
        archivo_seleccion,
        RAIZ / "apoyo.py",
        archivo_ratios,
        archivo_concentraciones,
        archivo_muestras,
        archivo_referencias,
        RAIZ / "imputadores.py",
        Path(__file__),
    ]
    series_analiticas = None
    if dataset == "ishii":
        fuente = RAIZ / "sources/ecoli_auditada_v3.db"
        entradas.append(fuente)
        with sqlite3.connect(
            fuente.resolve().as_uri() + "?mode=ro", uri=True
        ) as conexion:
            series_analiticas = (
                pd.read_sql_query(
                    (
                        "SELECT sample_id_original AS sample_id, "
                        "source_series_id AS series_raw, source_c"
                        "olumn FROM ishii_metabolite_source_colum"
                        "ns WHERE source_sheet='metabolites'"
                    ),
                    conexion,
                )
                .set_index("sample_id")
                .loc[ratios.index]
            )
        series_analiticas["series_family"] = (
            series_analiticas.series_raw.astype(str).str.replace(
                "*", "", regex=False
            )
        )
    huellas = {
        str(archivo.relative_to(RAIZ)): huella(archivo) for archivo in entradas
    }
    mascaras, verdades = crear_mascaras(
        ratios,
        concentraciones,
        muestras,
        series_analiticas,
        objetivos,
        huellas,
        dataset,
    )
    if verdades.empty:
        raise ValueError("No hay valores conocidos elegibles para evaluación")
    salida.mkdir(parents=True)
    validacion, imputacion = (salida / "validacion", salida / "imputacion")
    validacion.mkdir()
    imputacion.mkdir()
    guardar_json(validacion / "mascaras.json", mascaras)
    guardar_csv(validacion / "verdades.csv", verdades)
    if series_analiticas is not None:
        guardar_csv(
            validacion / "series_analiticas.csv",
            series_analiticas.reset_index(),
        )
    # Reutilizar un ajuste evita entrenar de nuevo para cada valor ocultado.
    modelos, ajustes = ({}, [])

    def obtener_modelo(metodo, excluidas):
        clave = (metodo, tuple(sorted(excluidas)))
        if clave not in modelos:
            entrenamiento = ratios.drop(index=excluidas)
            if set(entrenamiento.index) & set(excluidas):
                raise ValueError("Filtración de muestras")
            modelo = Imputador(metodo, parametros).fit(entrenamiento)
            modelo.audit_fit_id = len(ajustes)
            ajustes.append(
                dict(
                    fit_id=len(ajustes),
                    excluded=sorted(excluidas),
                    **modelo.metadata(),
                )
            )
            modelos[clave] = modelo
        return modelos[clave]

    # 1. Evaluar mediciones conocidas: ocultar, estimar y comparar.
    predicciones = []
    verdades_por_caso = verdades.set_index("case_id")
    for ensayo in mascaras["trials"]:
        consulta = ratios.loc[[ensayo["sample_id"]]].copy()
        consulta.loc[:, ensayo["hidden"]] = np.nan
        for metodo in METHODS:
            modelo = obtener_modelo(metodo, ensayo["excluded"])
            _, diagnosticos = modelo.transform(consulta)
            for registro in diagnosticos[
                diagnosticos.metabolite_id.isin(ensayo["score"])
            ].to_dict("records"):
                case_id = ensayo["trial_id"] + "|" + registro["metabolite_id"]
                medicion = verdades_por_caso.loc[case_id]
                registro.update(
                    fit_id=modelo.audit_fit_id,
                    case_id=case_id,
                    trial_id=ensayo["trial_id"],
                    experiment=ensayo["experiment"],
                    truth_R=float(medicion.truth_R),
                    truth_mM=float(medicion.truth_mM),
                    estimate_mM=registro["estimate_R"]
                    * referencias[registro["metabolite_id"]],
                )
                registro["absolute_error_R"] = abs(
                    registro["estimate_R"] - registro["truth_R"]
                )
                registro["absolute_error_mM"] = abs(
                    registro["estimate_mM"] - registro["truth_mM"]
                )
                predicciones.append(registro)
    predicciones = pd.DataFrame(predicciones)
    if len(predicciones) != len(verdades) * len(METHODS):
        raise ValueError("Evaluación incompleta")
    guardar_csv(validacion / "predicciones.csv", predicciones)
    metricas, motivos, comunes = evaluar_predicciones(predicciones)
    guardar_csv(validacion / "metricas.csv", metricas)
    guardar_csv(validacion / "cobertura_y_motivos.csv", motivos)
    guardar_csv(validacion / "casos_comunes.csv", comunes)
    guardar_json(
        validacion / "estado.json",
        dict(status="COMPLETE", dataset=dataset, scored_cells=len(verdades)),
    )
    print(
        f"{dataset}: evaluación completa, {len(verdades)} celdas por método.",
        flush=True,
    )
    # 2. Completar únicamente los huecos reales de los objetivos.
    filas_candidatos = []
    for muestra in ratios.index[ratios[objetivos].isna().any(axis=1)]:
        grupo = muestras.loc[muestra, "perturbation_group"]
        excluidas = list(
            ratios.index[
                muestras.loc[ratios.index, "perturbation_group"] == grupo
            ]
        )
        consulta = ratios.loc[[muestra]].copy()
        for metodo in METHODS:
            modelo = obtener_modelo(metodo, excluidas)
            _, diagnosticos = modelo.transform(consulta)
            for registro in diagnosticos[
                diagnosticos.metabolite_id.isin(objetivos)
            ].to_dict("records"):
                registro.update(
                    fit_id=modelo.audit_fit_id,
                    excluded_group=grupo,
                    estimate_mM=registro["estimate_R"]
                    * referencias[registro["metabolite_id"]],
                )
                filas_candidatos.append(registro)
    candidatos = pd.DataFrame(
        filas_candidatos,
        columns=list(
            dict.fromkeys(list(predicciones.columns) + ["excluded_group"])
        ),
    )
    guardar_csv(imputacion / "candidatos_y_soporte.csv", candidatos)
    if not candidatos.empty:
        estadisticas, pertenencias = resumir_candidatos(
            concentraciones, candidatos, imputacion
        )
        guardar_csv(
            imputacion / "estadisticas_antes_despues.csv", estadisticas
        )
        guardar_csv(imputacion / "casillas_comunes.csv", pertenencias)
    guardar_json(salida / "ajustes.json", ajustes)
    if huellas != {
        str(archivo.relative_to(RAIZ)): huella(archivo) for archivo in entradas
    }:
        raise ValueError("Entradas modificadas")
    guardar_json(
        salida / "manifest.json",
        dict(
            dataset=dataset,
            status="SUCCESS",
            config=parametros,
            targets=objetivos,
            input_sha256=huellas,
            python=platform.python_version(),
            versions={
                nombre: importlib.metadata.version(nombre)
                for nombre in ["numpy", "pandas", "scipy", "scikit-learn"]
            },
            validation_attempts=len(predicciones),
            real_attempts=len(candidatos),
            unavailable_tests=(
                []
                if series_analiticas is not None
                else ["analytical_series: no metadata homologous to Ishii"]
            ),
            excluded_evaluation_samples=(
                []
                if dataset == "ishii"
                else ["MMC2_GLUCOSE: reference reveals truth"]
            ),
            output_sha256={
                str(archivo.relative_to(salida)): huella(archivo)
                for archivo in sorted(salida.rglob("*"))
                if archivo.is_file()
            },
        ),
    )
    print(
        f"{dataset}: {len(candidatos)} candidatos. " f"Salida: {salida}",
        flush=True,
    )


if __name__ == "__main__":
    main()
