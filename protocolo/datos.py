"""Reconstrucción desde filas analíticas. No se imputan ni se fusionan señales aquí."""
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook

from .fisiologia import auditar
from .fuentes import EXTRAS, RAIZ, SEMILLA
from .proteinas import leer_proteinas

PRIV = Path(os.environ.get("IPRE_DETALLE", RAIZ / "salida_protocolo")).resolve()
PUB = RAIZ / "protocolo" / "resultados"
ANCLAS = ("fum", "pep", "s7p", "succ")


def huella(ruta):
    return hashlib.sha256(ruta.read_bytes()).hexdigest()


def guardar_json(ruta, dato):
    ruta.write_text(json.dumps(dato, indent=2, ensure_ascii=False,
                               sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def leer_matriz(nombre):
    return pd.read_csv(PRIV / f"{nombre}.csv", index_col=0, float_precision="round_trip")


def criterio_identidad(fila):
    """Una señal nominal puede servir dentro de Ishii sin definir una equivalencia externa."""
    if fila.fila_elegida_ishii == True:  # Solo correspondencias ya documentadas.
        return "compartida"
    if fila.composite or "+" in fila.original_label:
        return "mezcla_no_resuelta"
    if fila.source_row in {52, 132, 133, 524} or fila.nombre_repetido_en_base == True:
        return "identidad_duplicada_o_alternativa"
    if (fila.original_label.startswith(("Divalent ion", "?")) or
            fila.original_label in {"Glycerophosphorate", "Keramine"}):
        return "identidad_incierta"
    return "nominal_Ishii"


def preparar():
    PRIV.mkdir(parents=True, exist_ok=True)
    PUB.mkdir(parents=True, exist_ok=True)
    audit = RAIZ / "auditoria"
    meta = pd.read_csv(audit / "metadatos.csv", index_col=0).sort_index()
    ish = meta.index[meta.dataset_id.eq("ISHII2007") & meta.incluida]
    wt = meta.index[meta.dataset_id.eq("ISHII2007") & ~meta.incluida]
    aux = meta.index[meta.dataset_id.eq("MMC2") & meta.incluida & (meta.index != "MMC2_ALL")]
    assert (len(ish), len(wt), len(aux), meta.loc[ish, "grupo"].nunique()) == (30, 4, 8, 25)
    assert meta.loc[ish, "dilution_rate_h_inv"].eq(.2).all()
    h = pd.read_csv(audit / "homologacion_maestra.csv")
    h = h[h.analyte_id.str.startswith("ISHII:metabolites:", na=False) &
          h.compartment.eq("intracellular")].sort_values("source_row").copy()
    h["source_row"] = h.source_row.astype(int)
    assert len(h) == 579 and h.unidad_original.eq("mM").all()
    h["identidad"] = [criterio_identidad(r) for r in h.itertuples()]
    h["variable"] = [r.metabolito if r.identidad == "compartida" else
                     EXTRAS.get(r.source_row, (f"ish_{r.source_row}",))[0] for r in h.itertuples()]
    assert h.variable.is_unique
    comunes = h.loc[h.identidad.eq("compartida"), "variable"].tolist()
    medidas = pd.read_csv(audit / "mediciones_auditadas.csv", low_memory=False,
                          float_precision="round_trip")
    muestras = ish.append(wt)
    t = medidas[medidas.analyte_id.isin(h.analyte_id) & medidas.sample_key.isin(muestras)].copy()
    t["variable"] = t.analyte_id.map(h.set_index("analyte_id").variable)
    libros = {n: load_workbook(RAIZ / "sources" / f"{n}.xlsx", data_only=True)
              for n in ("ishii", "mmc2")}
    try:
        sh = libros["ishii"]["metabolites"]
        for f in h.itertuples():
            assert str(sh.cell(f.source_row, 1).value).strip() == f.original_label
        assert "not detectable" in libros["ishii"]["generals"]["A49"].value.lower()
        assert float(libros["mmc2"]["Conversion Factors"]["B6"].value) == 2.3
        g = medidas[medidas.dataset_id.eq("MMC2") & medidas.metabolito.isin(comunes) &
                    medidas.sample_key.isin(aux) & medidas.entity_type.eq("metabolite")].copy()
        g["variable"] = g.metabolito
        assert len(g) == 8 * 20 and len(t) == 34 * 579
        for tabla, libro, divisor in [(t, libros["ishii"], 1.), (g, libros["mmc2"], 2.3)]:
            assert tabla.is_imputed.eq(0).all() and tabla.divisor_conversion.eq(divisor).all()
            original = pd.to_numeric([libro[r.source_sheet][r.source_cell].value
                                      for r in tabla.itertuples()], errors="coerce")
            np.testing.assert_allclose(original / divisor, tabla.valor_mM, atol=1e-13, rtol=1e-13,
                                       equal_nan=True)
            assert tabla.valor_mM.dropna().ge(0).all() and np.isfinite(tabla.valor_mM.dropna()).all()
    finally:
        for libro in libros.values():
            libro.close()
    def matriz(tabla, campo, idx, cols):
        return tabla.pivot(index="sample_key", columns="variable", values=campo).loc[idx, cols]
    cruda = matriz(t, "valor_mM", muestras, h.variable)
    original = cruda.loc[ish]
    h["n_30"] = h.variable.map(original.count())
    h["n_grupos"] = h.variable.map(original.notna().groupby(meta.loc[ish, "grupo"]).any().sum())
    h["distintos"] = h.variable.map(original.nunique())
    nominal = h.identidad.isin(["nominal_Ishii", "compartida"])
    h["decision"] = np.select([~nominal, h.n_30.lt(10), h.n_grupos.lt(5), h.distintos.lt(2)],
                              [h.identidad, "menos_de_10", "menos_de_5_grupos", "sin_variacion"],
                              default="candidato_preliminar")
    cols = h.loc[h.decision.eq("candidato_preliminar"), "variable"]
    h[["source_row","variable","original_label","protocol","compartment","unidad_original",
       "identidad","n_30","n_grupos","distintos","decision"]].to_csv(PUB/"inventario.csv", index=False)
    cruda.to_csv(PRIV / "original_34_todas.csv")
    original[cols].to_csv(PRIV / "ishii.csv")
    cruda.loc[wt, cols].to_csv(PRIV / "wt.csv")
    matriz(g, "valor_mM", aux, comunes).to_csv(PRIV / "mmc2.csv")
    matriz(t, "absence_code", ish, cols).to_csv(PRIV / "ausencias.csv")
    meta.loc[muestras.append(aux)].to_csv(PRIV / "metadatos.csv")
    lin = ["sample_key","variable","analyte_id","source_sheet","source_cell","measurement_uid",
           "protocol","compartment","unit_original","value_original","divisor_conversion","valor_mM",
           "SD_mM","absence_code"]
    pd.concat([t, g])[lin].to_csv(PRIV/"trazabilidad.csv.gz", index=False,
                                    compression={"method":"gzip", "mtime":0})
    q = auditar(ish, wt, aux, meta, PRIV/"fisiologia")
    prot = leer_proteinas(ish, PRIV/"proteinas")
    pd.concat([q.q_C, prot.add_prefix("prot:")], axis=1).to_csv(PRIV/"auxiliares.csv")
    guardar_json(PRIV/"fuentes.json", {p.name:huella(p) for p in
                 [audit/"mediciones_auditadas.csv",audit/"homologacion_maestra.csv",audit/"metadatos.csv"] +
                 [RAIZ/"sources"/n for n in ("ishii.xlsx","mmc2.xlsx","mmc4.xlsx","ecoli_auditada_v3.db")]})
    print(f"{len(h)} señales, {len(cols)} objetivos preliminares, {original[cols].isna().sum().sum()} huecos.")


if __name__ == "__main__":
    preparar()
