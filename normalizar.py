"""Construye matrices en mM y ratios; no entrena ni imputa. Python estándar."""

import argparse
import csv
import json
import math
import sqlite3
from pathlib import Path

from apoyo import huella

RAIZ = Path(__file__).resolve().parent


def guardar_filas(path, filas, columnas):
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columnas, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(filas)


def main():
    argumentos = argparse.ArgumentParser(description=__doc__)
    argumentos.add_argument(
        "--output", type=Path, default=RAIZ / "normalizacion"
    )
    opciones = argumentos.parse_args()
    archivo_seleccion = RAIZ / "config_normalizacion.json"
    seleccion = json.loads(archivo_seleccion.read_text(encoding="utf-8"))
    entradas = [
        RAIZ / "data/concentraciones_mM.db",
        RAIZ / "sources/ecoli_auditada_v3.db",
        archivo_seleccion,
        Path(__file__).resolve(),
    ]
    huellas = {
        str(archivo.relative_to(RAIZ)): huella(archivo) for archivo in entradas
    }
    base = sqlite3.connect(
        entradas[0].resolve().as_uri() + "?mode=ro", uri=True
    )
    base.row_factory = sqlite3.Row
    fuente = sqlite3.connect(
        entradas[1].resolve().as_uri() + "?mode=ro", uri=True
    )
    fuente.row_factory = sqlite3.Row
    metabolitos = list(seleccion["ishii_analyte_ids"])
    objetivos = seleccion["targets"]
    registros = [
        dict(registro)
        for registro in base.execute(
            "SELECT * FROM metabolite_measurements_mM"
        )
    ]
    ishii = [
        registro
        for registro in registros
        if registro["dataset_id"] == "ISHII2007"
        and registro["compartment"] == "intracellular"
        and (registro["metabolite_id"] in metabolitos)
        and (
            registro["analyte_id"]
            == seleccion["ishii_analyte_ids"][registro["metabolite_id"]]
        )
    ]
    metadatos = {
        registro["sample_id"]: dict(registro)
        for registro in fuente.execute("SELECT * FROM ishii_samples")
    }
    disponibles = sorted({registro["sample_key"] for registro in ishii})
    cohorte = [
        muestra
        for muestra in disponibles
        if metadatos.get(muestra, {}).get("dilution_rate_h_inv") == 0.2
        or muestra in seleccion["source_only_samples"]
    ]
    if len(cohorte) != 30:
        raise ValueError("Se esperaban 30 muestras Ishii")
    if not ("KO05x" in cohorte and "KO05" in cohorte):
        raise ValueError("Mantener ambos cultivos pfkA")
    referencias = {}
    for metabolito in metabolitos:
        seleccionados = [
            registro
            for registro in registros
            if registro["dataset_id"] == "MMC2"
            and registro["condition_id"] == seleccion["reference_condition"]
            and (registro["metabolite_id"] == metabolito)
            and (registro["compartment"] == "intracellular")
        ]
        if len(seleccionados) != 1:
            raise ValueError("Referencia duplicada o ausente: " + metabolito)
        registro = seleccionados[0]
        if not (
            registro["value_mM"] is not None
            and math.isfinite(registro["value_mM"])
            and (registro["value_mM"] > 0)
        ):
            raise ValueError("Referencia no positiva: " + metabolito)
        if not (
            registro["absence_code"] == "reported"
            and (not registro["is_imputed"])
            and (not registro["composite"])
        ):
            raise ValueError("Referencia no admisible: " + metabolito)
        if not (
            registro["uncertainty_type"] == "SD"
            and registro["uncertainty_mM"] is not None
            and (registro["uncertainty_mM"] >= 0)
            and registro["source_uncertainty_cell"]
        ):
            raise ValueError("SD sin documentar: " + metabolito)
        referencias[metabolito] = registro
    condiciones = {
        registro["condition_id"]: dict(registro)
        for registro in fuente.execute(
            (
                "SELECT * FROM conditions WHERE dataset_i"
                "d='MMC2' AND condition_type='carbon_sour"
                "ce'"
            )
        )
    }
    registros_mmc2 = []
    filas_condiciones = []
    for condicion in seleccion["mmc2_conditions"]:
        panel = [
            registro
            for registro in registros
            if registro["dataset_id"] == "MMC2"
            and registro["condition_id"] == condicion
            and (registro["metabolite_id"] in metabolitos)
        ]
        if not (
            condicion in condiciones
            and len(panel) == 7
            and (len({registro["metabolite_id"] for registro in panel}) == 7)
        ):
            raise ValueError("Panel MMC2 incompleto/duplicado")
        if not all(
            (
                registro["compartment"] == "intracellular"
                and (not registro["composite"])
                and (not registro["is_imputed"])
                and (
                    registro["mapping_status"]
                    in ("exact_base_unique", "confirmed_intracellular_source")
                )
                for registro in panel
            )
        ):
            raise ValueError("Identidad MMC2 no admisible")
        registros_mmc2.extend(panel)
        filas_condiciones.append(
            dict(
                condiciones[condicion],
                admitted=True,
                observed=sum(
                    (registro["value_mM"] is not None for registro in panel)
                ),
                missing=sum(
                    (registro["value_mM"] is None for registro in panel)
                ),
                reference_row=condicion == seleccion["reference_condition"],
                reason=(
                    "Identidad resuelta, intracelular, estaci"
                    "onario; ausencias se conservan. No impli"
                    "ca equivalencia biológica con Ishii."
                ),
                evaluation_role="independent_MMC2_analysis_pending",
            )
        )
    seleccionados = [
        registro for registro in ishii if registro["sample_key"] in cohorte
    ] + registros_mmc2
    celdas = {}
    trazabilidad = []
    for registro in seleccionados:
        muestra = (
            registro["sample_key"]
            if registro["dataset_id"] == "ISHII2007"
            else registro["condition_id"]
        )
        metabolito = registro["metabolite_id"]
        referencia = referencias[metabolito]
        clave = (muestra, metabolito)
        if clave in celdas:
            raise ValueError("Celda duplicada: " + str(clave))
        valor = registro["value_mM"]
        ratio = None if valor is None else valor / referencia["value_mM"]
        if not (
            valor is None
            or (
                math.isfinite(valor)
                and valor >= 0
                and math.isclose(
                    ratio * referencia["value_mM"],
                    valor,
                    rel_tol=1e-12,
                    abs_tol=1e-14,
                )
            )
        ):
            raise ValueError("Reversibilidad/valor inválido")
        celdas[clave] = (valor, ratio)
        trazabilidad.append(
            dict(
                sample_id=muestra,
                dataset_id=registro["dataset_id"],
                metabolite_id=metabolito,
                role="target" if metabolito in objetivos else "predictor",
                value_mM=valor,
                ratio=ratio,
                observed=valor is not None,
                absence_code=registro["absence_code"],
                measurement_uid=registro["measurement_uid"],
                analyte_id=registro["analyte_id"],
                protocol=registro["protocol"],
                source_file=registro["source_file"],
                source_sheet=registro["source_sheet"],
                source_cell=registro["source_cell"],
                source_sha256=registro["source_file_sha256"],
                uncertainty_mM=registro["uncertainty_mM"],
                uncertainty_type=registro["uncertainty_type"],
                reference_uid=referencia["measurement_uid"],
                reference_mM=referencia["value_mM"],
                reference_SD_mM=referencia["uncertainty_mM"],
                is_reference_cell=registro["measurement_uid"]
                == referencia["measurement_uid"],
                ratio_uncertainty_status="not_propagated_shared_denominator",
                is_imputed=False,
            )
        )
    if len(celdas) != 266:
        raise ValueError("Se esperaban 266 celdas")
    if not all(
        (
            celdas[seleccion["reference_condition"], metabolito][1] == 1
            for metabolito in metabolitos
        )
    ):
        raise ValueError("Referencia normalizada distinta de 1")
    filas_muestras = []
    for muestra in disponibles:
        metadatos_muestra = metadatos.get(muestra, {})
        filas_muestras.append(
            dict(
                sample_id=muestra,
                dataset_id="ISHII2007",
                included=muestra in cohorte,
                perturbation_group=(
                    "pfkA"
                    if muestra in ("KO05", "KO05x")
                    else metadatos_muestra.get("knockout_gene") or "WT"
                ),
                dilution_h_inv=metadatos_muestra.get("dilution_rate_h_inv"),
                carbon_source="glucose",
                protocol_family="mixed_analyte_protocols",
                metadata_status=(
                    "source_only_culture_no_S2_row"
                    if muestra == "KO05x"
                    else "documented"
                ),
                reason=(
                    "retain_previous_cohort_pending_metadata"
                    if muestra == "KO05x"
                    else (
                        "dilution_0.2"
                        if muestra in cohorte
                        else "other_dilution_excluded"
                    )
                ),
            )
        )
    for registro in filas_condiciones:
        filas_muestras.append(
            dict(
                sample_id=registro["condition_id"],
                dataset_id="MMC2",
                included=True,
                perturbation_group=registro["condition_id"],
                dilution_h_inv=None,
                carbon_source=registro["carbon_source"],
                protocol_family="Gerosa",
                metadata_status="documented",
                reason="independent_MMC2_condition",
            )
        )
    resumen = {}
    salida = opciones.output.resolve()
    salida.mkdir(parents=True, exist_ok=False)
    for estudio, muestras in [
        ("ishii", cohorte),
        ("mmc2", seleccion["mmc2_conditions"]),
    ]:
        for unidad, posicion in [("mM", 0), ("R", 1), ("observado", None)]:
            matriz = [
                dict(
                    sample_id=muestra,
                    **{
                        metabolito: (
                            int(celdas[muestra, metabolito][0] is not None)
                            if posicion is None
                            else celdas[muestra, metabolito][posicion]
                        )
                        for metabolito in metabolitos
                    },
                )
                for muestra in muestras
            ]
            guardar_filas(
                salida / f"matriz_{estudio}_{unidad}.csv",
                matriz,
                ["sample_id"] + metabolitos,
            )
        observados = sum(
            (
                celdas[muestra, metabolito][0] is not None
                for muestra in muestras
                for metabolito in metabolitos
            )
        )
        objetivos_observados = sum(
            (
                celdas[muestra, metabolito][0] is not None
                for muestra in muestras
                for metabolito in objetivos
            )
        )
        resumen[estudio] = dict(
            rows=len(muestras),
            columns=7,
            observed=observados,
            missing=len(muestras) * 7 - observados,
            target_observed=objetivos_observados,
            target_missing=len(muestras) * 4 - objetivos_observados,
        )
    guardar_filas(
        salida / "trazabilidad_celdas.csv", trazabilidad, list(trazabilidad[0])
    )
    guardar_filas(
        salida / "muestras_y_exclusiones.csv",
        filas_muestras,
        list(filas_muestras[0]),
    )
    guardar_filas(
        salida / "condiciones_mmc2.csv",
        filas_condiciones,
        list(filas_condiciones[0]),
    )
    filas_referencias = [
        dict(
            metabolite_id=metabolito,
            reference_condition=seleccion["reference_condition"],
            value_mM=registro["value_mM"],
            SD_mM=registro["uncertainty_mM"],
            uncertainty_type=registro["uncertainty_type"],
            measurement_uid=registro["measurement_uid"],
            source_file=registro["source_file"],
            source_sheet=registro["source_sheet"],
            source_cell=registro["source_cell"],
            source_SD_cell=registro["source_uncertainty_cell"],
            source_sha256=registro["source_file_sha256"],
        )
        for metabolito, registro in referencias.items()
    ]
    guardar_filas(
        salida / "referencias_glucosa.csv",
        filas_referencias,
        list(filas_referencias[0]),
    )
    if huellas != {
        str(archivo.relative_to(RAIZ)): huella(archivo) for archivo in entradas
    }:
        raise ValueError("Entradas modificadas")
    manifest = dict(
        status="SUCCESS",
        stage="normalized_no_imputation",
        summary=resumen,
        input_sha256=huellas,
        checks=dict(
            positive_observed_references=True,
            reference_SD_documented=True,
            unique_cells=True,
            missing_preserved=True,
            reversible=True,
            glucose_ratios_one=True,
            inputs_unchanged=True,
        ),
        output_sha256={
            archivo.name: huella(archivo)
            for archivo in sorted(salida.iterdir())
            if archivo.is_file()
        },
    )
    (salida / "manifest_normalizacion.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(resumen, indent=2))
    base.close()
    fuente.close()


if __name__ == "__main__":
    main()
