"""Coteja las concentraciones materializadas con las dos bases SQLite, en modo de lectura."""
import sqlite3
from contextlib import closing

import numpy as np
import pandas as pd

from protocolo.fuentes import RAIZ


def verificar():
    archivos = RAIZ / "auditoria"
    medidas = pd.read_csv(archivos / "mediciones_auditadas.csv", index_col="measurement_uid",
                          low_memory=False, float_precision="round_trip").sort_index()
    origen = RAIZ / "sources" / "ecoli_auditada_v3.db"
    trabajo = RAIZ / "data" / "concentraciones_mM.db"
    with closing(sqlite3.connect(origen.as_uri() + "?mode=ro", uri=True)) as db:
        assert db.execute("PRAGMA integrity_check").fetchone() == ("ok",)
        originales = pd.read_sql_query(
            "SELECT measurement_uid, analyte_id, sample_key, value_original, unit_original, "
            "absence_code, source_sheet, source_cell, uncertainty_original, is_imputed "
            "FROM audit_measurements_readable WHERE entity_type='metabolite'",
            db, index_col="measurement_uid").sort_index()
    with closing(sqlite3.connect(trabajo.as_uri() + "?mode=ro", uri=True)) as db:
        assert db.execute("PRAGMA integrity_check").fetchone() == ("ok",)
        convertidas = pd.read_sql_query(
            "SELECT measurement_uid, value_mM, uncertainty_mM, conversion_divisor "
            "FROM metabolite_measurements_mM", db, index_col="measurement_uid").sort_index()

    assert medidas.index.is_unique and originales.index.is_unique and convertidas.index.is_unique
    pd.testing.assert_index_equal(medidas.index, originales.index)
    pd.testing.assert_index_equal(medidas.index, convertidas.index)
    for campo in ("analyte_id", "sample_key", "unit_original", "absence_code",
                  "source_sheet", "source_cell", "is_imputed"):
        pd.testing.assert_series_equal(medidas[campo], originales[campo], check_dtype=False)
    for campo in ("value_original", "uncertainty_original"):
        np.testing.assert_allclose(medidas[campo], originales[campo], equal_nan=True,
                                   rtol=1e-12, atol=1e-14)
    for campo, sql in (("valor_mM", "value_mM"), ("SD_mM", "uncertainty_mM"),
                       ("divisor_conversion", "conversion_divisor")):
        np.testing.assert_allclose(medidas[campo], convertidas[sql], equal_nan=True,
                                   rtol=1e-12, atol=1e-14)
    return len(medidas)


if __name__ == "__main__":
    print(f"{verificar()} posiciones cotejadas por identificador en ambas bases.")
