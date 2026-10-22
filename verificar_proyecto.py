"""Verifica esta distribución en solo lectura; no normaliza ni imputa."""

import csv
import hashlib
import json
import math
import sqlite3
from pathlib import Path

RAIZ = Path(__file__).resolve().parent


def abrir_base(archivo):
    return sqlite3.connect(archivo.resolve().as_uri() + "?mode=ro", uri=True)


def main():
    inventario = json.loads(
        (RAIZ / "MANIFEST_PROYECTO.json").read_text(encoding="utf-8")
    )
    for nombre, esperado in inventario["files_sha256"].items():
        archivo = RAIZ / nombre
        if not archivo.is_file():
            raise ValueError("Falta: " + nombre)
        if hashlib.sha256(archivo.read_bytes()).hexdigest() != esperado:
            raise ValueError("Archivo modificado: " + nombre)
    original = json.loads((RAIZ / "manifest.json").read_text(encoding="utf-8"))
    for key, carpeta in [
        ("source_sha256", "sources"),
        ("output_sha256", "data"),
    ]:
        for nombre, esperado in original[key].items():
            if (
                hashlib.sha256(
                    (RAIZ / carpeta / nombre).read_bytes()
                ).hexdigest()
                != esperado
            ):
                raise ValueError(
                    "Diferencia respecto a conversión original: " + nombre
                )
    with abrir_base(RAIZ / "data/concentraciones_mM.db") as base:
        if base.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
            raise ValueError("Integridad SQLite")
        if base.execute("PRAGMA foreign_key_check").fetchall():
            raise ValueError("Claves foráneas")
        filas = base.execute(
            (
                "SELECT dataset_id, value_original, value"
                "_mM, uncertainty_original, uncertainty_m"
                "M, conversion_divisor, is_imputed FROM m"
                "etabolite_measurements_mM"
            )
        ).fetchall()
        if len(filas) != 20483:
            raise ValueError("Conteo inesperado")
        for estudio, valor, convertido, sd, sd_mM, divisor, imputado in filas:
            if divisor != (2.3 if estudio in ("MMC2", "MMC3") else 1.0):
                raise ValueError("Divisor incorrecto")
            if imputado != 0:
                raise ValueError("Imputación inesperada")
            for a, b in [(valor, convertido), (sd, sd_mM)]:
                if (a is None) != (b is None):
                    raise ValueError("Patrón de faltantes modificado")
                if a is not None:
                    if not math.isclose(
                        a / divisor, b, rel_tol=1e-12, abs_tol=1e-14
                    ):
                        raise ValueError("Conversión incorrecta")
        conteos = base.execute(
            (
                "SELECT dataset_id, compartment, count(*)"
                ", count(value_mM) FROM metabolite_measur"
                "ements_mM GROUP BY 1,2"
            )
        ).fetchall()
        claves_muestra = {
            x[0]
            for x in base.execute(
                (
                    "SELECT DISTINCT sample_key FROM metaboli"
                    "te_measurements_mM WHERE dataset_id='ISH"
                    "II2007'"
                )
            )
        }
    with abrir_base(RAIZ / "sources/ecoli_auditada_v3.db") as fuente:
        if fuente.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
            raise ValueError("Integridad fuente")
        consulta = fuente.execute(
            "SELECT * FROM ishii_samples ORDER BY sample_id"
        )
        filas_esperadas = [[x[0] for x in consulta.description]] + [
            [str(v) if v is not None else "" for v in row] for row in consulta
        ]
    with (RAIZ / "data/samples_ishii.csv").open(
        encoding="utf-8-sig", newline=""
    ) as file:
        actual = list(csv.reader(file))
    if actual != filas_esperadas:
        raise ValueError("Metadatos diferentes de la fuente")
    identificadores = {row[actual[0].index("sample_id")] for row in actual[1:]}
    sin_metadatos = claves_muestra - identificadores
    if not sin_metadatos <= {"KO05x"}:
        raise ValueError(
            "Muestras inesperadas sin metadatos: " + str(sin_metadatos)
        )
    if sin_metadatos:
        print(
            (
                "PENDIENTE DOCUMENTADO: KO05x tiene medic"
                "iones, pero no fila en ishii_samples; no"
                " se fusiona con KO05."
            )
        )
    print(
        "OK: hashes, SQLite, valores, incertidumbres, faltantes y metadatos."
    )
    print("Conjunto / compartimento / posiciones / valores observados")
    for row in conteos:
        print(" / ".join(map(str, row)))
    print(
        (
            "Estado: conversión verificada. Consultar"
            " normalizacion/manifest_normalizacion.js"
            "on para las matrices; evaluación e imput"
            "ación separadas en resultados/<dataset>/"
            "; consultar docs/VALIDACION_Y_RESULTADOS"
            ".md."
        )
    )


if __name__ == "__main__":
    main()
