"""Reconstruye el panel en mM y verifica cada posición contra su Excel original."""
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook

RAIZ = Path(__file__).resolve().parents[1]
SALIDA = RAIZ / "ampliacion" / "resultados"
AUDITORIA = RAIZ / "auditoria"

# Son señales de Ishii. Ninguna de estas claves designa una medición de Gerosa.
EXTRAS = {
    174: ("ish_g6p", "Glucose 6-phosphate"),
    173: ("ish_f6p", "Fructose 6-phosphate"),
    164: ("ish_ru5p", "Ribulose 5-phosphate"),
    53: ("ish_mal", "Malate"),
    124: ("ish_3pg", "3-Phosphoglycerate"),
    172: ("ish_g1p", "Glucose 1-phosphate"),
    165: ("ish_r5p", "Ribose 5-phosphate"),
}


def preparar():
    """Devuelve matrices separadas para los 30 cultivos, ocho auxiliares y cuatro WT."""
    SALIDA.mkdir(parents=True, exist_ok=True)
    historica = pd.read_csv(AUDITORIA / "panel_documental_mM.csv", index_col=0)
    meta = pd.read_csv(AUDITORIA / "metadatos.csv", index_col=0)
    inventario = pd.read_csv(AUDITORIA / "homologacion_maestra.csv")
    medidas = pd.read_csv(AUDITORIA / "mediciones_auditadas.csv", low_memory=False)
    ref = pd.read_csv(AUDITORIA / "referencias_glucosa.csv", index_col=0)
    comp = historica.columns
    ish = historica.index[meta.loc[historica.index, "dataset_id"].eq("ISHII2007")]
    mmc = historica.index.difference(ish, sort=False)
    wt = meta.index[meta.dataset_id.eq("ISHII2007") & ~meta.index.isin(ish)]
    assert (len(ish), len(mmc), len(wt)) == (30, 8, 4)
    assert meta.loc[ish, "dilution_rate_h_inv"].eq(0.2).all()
    assert meta.loc[wt, "grupo"].eq("WT").all() and "MMC2_ALL" not in mmc
    assert ref.loc[comp, "valida"].all()
    assert (ref.loc[comp, "valor_mM"] > 0).all() and (ref.loc[comp, "SD_mM"] >= 0).all()
    assert np.isfinite(ref.loc[comp, ["valor_mM", "SD_mM"]]).all().all()

    catalogo = inventario.query("dataset_id == 'ISHII2007' and compartment == 'intracellular'").copy()
    catalogo = catalogo.set_index("analyte_id")
    comunes = catalogo[catalogo.fila_elegida_ishii.eq(True)].set_index("metabolito")
    fuentes = {f"ISHII:metabolites:{int(comunes.at[m, 'source_row'])}": m for m in comp}
    for fila, (clave, nombre) in EXTRAS.items():
        analito = f"ISHII:metabolites:{fila}"
        assert catalogo.at[analito, "original_label"] == nombre
        assert catalogo.at[analito, "protocol"] == "Anion"
        assert not catalogo.at[analito, "composite"]
        fuentes[analito] = clave
    columnas = list(fuentes.values())
    fuentes.update({f"MMC2:{m}": m for m in comp})
    muestras = ish.append(wt).append(mmc)
    filas = medidas[medidas.analyte_id.isin(fuentes) & medidas.sample_key.isin(muestras)].copy()
    filas["variable"] = filas.analyte_id.map(fuentes)
    assert filas.compartment.eq("intracellular").all() and filas.is_imputed.eq(0).all()
    assert not filas.duplicated(["sample_key", "variable"]).any()
    assert len(filas) == 34 * len(columnas) + 8 * len(comp)
    libros = {d: load_workbook(RAIZ / "sources" / nombre, data_only=True)
              for d, nombre in [("ISHII2007", "ishii.xlsx"), ("MMC2", "mmc2.xlsx")]}
    volumen = libros["MMC2"]["Conversion Factors"]["B6"].value
    assert volumen == 2.3
    for f in filas.itertuples():
        hoja = libros[f.dataset_id][f.source_sheet]
        valor = hoja[f.source_cell].value
        divisor = volumen if f.dataset_id == "MMC2" else 1
        assert f.unit_original == ("umol gCDW^-1" if f.dataset_id == "MMC2" else "mM")
        if pd.notna(f.valor_mM):
            assert f.absence_code == "reported" and np.isfinite(f.valor_mM) and f.valor_mM >= 0
            np.testing.assert_allclose(f.valor_mM, valor / divisor, rtol=1e-12, atol=1e-14)
        else:
            assert f.absence_code != "reported" and (valor is None or str(valor).lower() == "nan")
        if pd.notna(f.SD_mM):
            np.testing.assert_allclose(f.SD_mM, hoja[f.source_uncertainty_cell].value / divisor)

    matriz = filas.pivot(index="sample_key", columns="variable", values="valor_mM")
    matriz = matriz.loc[muestras, columnas].rename_axis(index="sample_id", columns=None)
    pd.testing.assert_frame_equal(matriz.loc[historica.index, comp], historica, check_exact=False,
                                  rtol=1e-12, atol=1e-14)
    assert matriz.loc[mmc, list(columnas[len(comp):])].isna().all().all()
    gerosa = inventario[inventario.dataset_id.eq("MMC2")].set_index("metabolito")
    fichas = []
    for analito, clave in fuentes.items():
        if analito.startswith("MMC2:"):
            continue
        original = catalogo.loc[analito]
        pareja = clave if clave in comp else original.contraparte_propuesta
        g = gerosa.loc[pareja] if pd.notna(pareja) and pareja in gerosa.index else None
        fila_ishii = int(original.source_row)
        assert libros["ISHII2007"]["metabolites"][f"A{fila_ishii}"].value.strip() == original.original_label
        if g is not None:
            hoja = libros["MMC2"]["Metabolite concentrations"]
            assert hoja[f"B{int(g.source_row)}"].value.strip() == pareja
            assert hoja[f"C{int(g.source_row)}"].value.strip() == g.nombre_homologado
        fichas.append(dict(
            variable=clave, rol="objetivo_y_predictor" if clave in comp else "solo_predictor_Ishii",
            identidad_Ishii=original.original_label, fuente_Ishii=analito, protocolo_Ishii=original.protocol,
            compartimento=original.compartment, unidad_Ishii="mM",
            etiqueta_MMC2=pareja, identidad_MMC2=g.nombre_homologado if g is not None else None,
            fila_MMC2=int(g.source_row) if g is not None else None,
            protocolo_MMC2="LC-MS/MS" if g is not None else None,
            unidad_MMC2="µmol/gCDW" if g is not None else None,
            equivalencia=original.estado_homologacion, fundamento=original.justificacion,
            n_Ishii_30=matriz.loc[ish, clave].count(), n_WT_extra=matriz.loc[wt, clave].count(),
            n_MMC2_etiqueta_propuesta=(int(medidas.loc[
                medidas.analyte_id.eq(f"MMC2:{pareja}") & medidas.sample_key.isin(mmc), "valor_mM"
            ].count()) if g is not None else None),
            n_MMC2_utilizables=matriz.loc[mmc, clave].count(),
            referencia_mM=ref.at[clave, "valor_mM"] if clave in comp else None,
            SD_referencia_mM=ref.at[clave, "SD_mM"] if clave in comp else None,
        ))
    for libro in libros.values():
        libro.close()
    pd.DataFrame(fichas).to_csv(SALIDA / "diccionario.csv", index=False)
    filas[["sample_key", "variable", "analyte_id", "source_id", "source_sheet", "source_cell",
           "measurement_uid", "absence_code", "unit_original", "divisor_conversion",
           "valor_mM", "source_uncertainty_cell", "SD_mM"]].to_csv(SALIDA / "trazabilidad.csv", index=False)
    matriz.to_csv(SALIDA / "matriz_mM.csv")
    return matriz.loc[ish], matriz.loc[mmc], matriz.loc[wt], meta, ref.loc[comp, "valor_mM"]


if __name__ == "__main__":
    preparar()
