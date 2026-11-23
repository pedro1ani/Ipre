"""Lee captación, flujos y dilución por identificador de cultivo y conserva las celdas fuente."""
import sqlite3
from contextlib import closing

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter

from .preparar import RAIZ, SALIDA

REVISION = SALIDA / "revision"
UNIDAD_Q = "mmol sustrato/(g peso seco·h)"
UNIDAD_V = "mmol/(g peso seco·h)"
# Número de carbonos de la molécula suministrada; las fichas KEGG se identifican en la salida.
CARBONO = {"acetate": (2, "C00033"), "fructose": (6, "C00095"), "galactose": (6, "C00124"),
           "glucose": (6, "C00031"), "glycerol": (3, "C00116"), "gluconate": (6, "C00257"),
           "pyruvate": (3, "C00022"), "succinate": (4, "C00042")}


def leer_fila(hoja, fila, columnas, ids, tipo):
    """Distingue un cero publicado de una celda vacía, un ID ausente y un guion."""
    datos = []
    for muestra in ids:
        col = columnas.get(muestra)
        coord = f"{get_column_letter(col)}{fila}" if col else None
        bruto = hoja[coord].value if coord else None
        valor, sd = np.nan, np.nan
        if col is None:
            estado = "sin_columna_para_este_ID"
        elif bruto is None:
            estado = "vacia"
        elif str(bruto).strip() == "-":
            estado = "excluida_del_modelo" if tipo.startswith(("Ishii_flujo", "Ishii_coef")) else "guion_fuente"
        elif str(bruto).strip().lower() == "nan":
            estado = "no_reportada"
        else:
            texto = str(bruto).strip().split("±")
            valor = float(texto[0])
            sd = float(texto[1]) if len(texto) == 2 else np.nan
            assert np.isfinite(valor) and (np.isnan(sd) or sd >= 0)
            estado = "reportada"
        datos.append(dict(muestra=muestra, celda=coord, valor=valor, SD=sd, estado=estado))
    return pd.DataFrame(datos).set_index("muestra")


def auditar(ish, wt, mmc, meta, salida=REVISION):
    salida.mkdir(parents=True, exist_ok=True)
    ids = ish.append(wt)
    muestras = ids.append(mmc)
    fisiologia = pd.DataFrame(index=muestras).rename_axis("muestra")
    fisiologia["cohorte"] = ["Ishii30"] * len(ish) + ["WT4"] * len(wt) + ["MMC2"] * len(mmc)
    fuentes = RAIZ / "sources"
    libros = {n: load_workbook(fuentes / (n + ".xlsx"), data_only=True) for n in ("ishii", "mmc2")}
    w = libros["ishii"]
    assert w["generals"]["B11"].value == "%-substrate uptake"
    assert w["generals"]["B12"].value == "mmol/g-dry cell weight/h"
    assert "excluded from model" in w["generals"]["A64"].value
    assert w["final fluxes"]["A2"].value.strip() == "Reactions  (mmol/gdcw/h)"
    assert w["fluxes"]["A12"].value == w["final fluxes"]["A12"].value == "G6P -> 6PG"
    trazas = []

    def registrar(libro, hoja, fila, columnas, indice, tipo, unidad, clave=None):
        sh = libros[libro][hoja]
        t = leer_fila(sh, fila, columnas, indice, tipo)
        t["cohorte"] = fisiologia.loc[t.index, "cohorte"]
        trazas.append(t.assign(archivo=libro + ".xlsx", hoja=hoja, fila=fila,
                               etiqueta=sh.cell(fila, 1).value.strip(), tipo=tipo, unidad=unidad))
        if clave is not None:
            fisiologia.loc[indice, clave] = t.valor
            fisiologia.loc[indice, "celda_" + clave] = t.celda
        return t

    for hoja, filas, tipo in [("rates", (3, 4, 5, 7, 8, 9, 10, 11, 12), "Ishii_tasa"),
                              ("fluxes", range(3, 46), "Ishii_flujo_relativo"),
                              ("fluxes", range(47, 54), "Ishii_coef_intercambio"),
                              ("final fluxes", range(3, 46), "Ishii_flujo_absoluto")]:
        sh = w[hoja]
        encabezados = [c for c in sh[1] if c.value in ids]
        columnas = {c.value: c.column for c in encabezados}
        assert len(encabezados) == len(columnas)
        unidad = {"Ishii_tasa": UNIDAD_V, "Ishii_flujo_relativo": "% de captación de glucosa",
                  "Ishii_flujo_absoluto": UNIDAD_V, "Ishii_coef_intercambio": "adimensional"}[tipo]
        for fila in filas:
            clave = {("rates", 3): "q_sustrato", ("fluxes", 12): "v_rel",
                     ("final fluxes", 12): "v_abs"}.get((hoja, fila))
            registrar("ishii", hoja, fila, columnas, ids, tipo, unidad, clave)

    fis = libros["mmc2"]["Physiology"]
    columnas = {"MMC2_" + c.value.upper(): c.column for c in fis[1][1:9]}
    assert set(columnas) == set(mmc)
    assert "mmol" in fis["A4"].value and "Uptake rate carbon source" in fis["A4"].value
    for fila, unidad, clave in [(2, "1/h", "mu"), (4, UNIDAD_Q, "q_sustrato"),
                                (6, UNIDAD_V, None), (8, UNIDAD_V, None),
                                (10, UNIDAD_V, None), (11, "g/g", None)]:
        t = registrar("mmc2", "Physiology", fila, columnas, mmc, "MMC2_fisiologia", unidad, clave)
        if clave:
            fisiologia.loc[mmc, "SD_" + clave] = t.SD
    flujo = libros["mmc2"]["Metabolic fluxes"]
    columnas = {"MMC2_" + c.value.upper(): c.column for c in flujo[2][1:9]}
    assert set(columnas) == set(mmc) and flujo["A38"].value == "Growth_rate"
    for fila in range(3, 38):  # La fila 38 repite crecimiento; no es otro flujo químico.
        t = registrar("mmc2", "Metabolic fluxes", fila, columnas, mmc, "MMC2_flujo", UNIDAD_V)
        for muestra, col in columnas.items():
            col_sd = col + 9
            assert flujo.cell(2, col).value == flujo.cell(2, col_sd).value
            sd = flujo.cell(fila, col_sd)
            trazas[-1].at[muestra, "celda_SD"] = sd.coordinate
            trazas[-1].at[muestra, "SD"] = float(sd.value) if isinstance(sd.value, (int, float)) else np.nan

    # El dato primario es la etiqueta o el rango de muestras de Ishii; SQLite y el CSV se cotejan.
    with closing(sqlite3.connect(f"file:{fuentes / 'ecoli_auditada_v3.db'}?mode=ro", uri=True)) as db:
        d = pd.read_sql_query("SELECT sample_id, dilution_rate_h_inv FROM ishii_samples", db,
                              index_col="sample_id").dilution_rate_h_inv
        aux = pd.read_sql_query("SELECT condition_id, variable_original, value, sd "
                                "FROM physiology_measurements WHERE dataset_id='MMC2'", db)
    en_sql = ids.intersection(d.index)
    assert set(ids.difference(en_sql)) == {"KO05x"}  # El libro sí identifica el primer cultivo pfkA.
    np.testing.assert_allclose(d.loc[en_sql], meta.loc[en_sql, "dilution_rate_h_inv"])
    for s in ids:
        dilucion = meta.at[s, "dilution_rate_h_inv"]
        col = meta.at[s, "source_column"]
        assert w["metabolites"][col + "2"].value == s
        if s.startswith("GR"):
            fuente = "metabolites!" + col + "4"
            assert str(dilucion) + "h-1" in w["metabolites"][col + "4"].value
        else:
            fuente = "generals!C20" if s.startswith("KO") else "generals!C22"
            assert "0.2h-1" in w["generals"][fuente.split("!")[1]].value.replace(" ", "") and dilucion == 0.2
        fisiologia.at[s, "D"] = dilucion
        fisiologia.at[s, "fuente_D"] = fuente
    np.testing.assert_allclose(fisiologia.loc[mmc, "mu"], meta.loc[mmc, "growth_rate"])
    for nombre, fila in [("mu", 2), ("q_sustrato", 4)]:
        x = aux[aux.variable_original.str.strip().eq(fis.cell(fila, 1).value.strip())].set_index("condition_id")
        np.testing.assert_allclose(fisiologia.loc[mmc, nombre], x.loc[mmc, "value"])
        np.testing.assert_allclose(fisiologia.loc[mmc, "SD_" + nombre], x.loc[mmc, "sd"])

    fisiologia["sustrato"] = meta.loc[muestras, "carbono"].replace({"glucosa": "glucose"})
    fisiologia["carbonos"] = [CARBONO[s][0] for s in fisiologia.sustrato]
    fisiologia["KEGG_sustrato"] = [CARBONO[s][1] for s in fisiologia.sustrato]
    fisiologia["q_C"] = fisiologia.q_sustrato * fisiologia.carbonos
    fisiologia["SD_q_C"] = fisiologia.SD_q_sustrato * fisiologia.carbonos
    assert fisiologia.q_sustrato.dropna().gt(0).all()
    pub = fisiologia.v_abs.notna()
    np.testing.assert_allclose(fisiologia.loc[pub, "v_abs"],
                               (fisiologia.q_sustrato * fisiologia.v_rel / 100)[pub], atol=1e-10)
    bruto = pd.concat(trazas).reset_index()
    bruto.to_csv(salida / "fuentes_fisiologicas.csv.gz", index=False,
                 compression={"method": "gzip", "mtime": 0})
    inventario = bruto.groupby(["archivo", "hoja", "fila", "etiqueta", "tipo", "unidad", "cohorte"]).valor
    inventario.agg(posiciones="size", disponibles="count", minimo="min", maximo="max").to_csv(
        salida / "inventario.csv")
    fisiologia.to_csv(salida / "condiciones.csv")
    for libro in libros.values():
        libro.close()
    return fisiologia
