"""Abundancias por cultivo y selección de una proteína de control en el entrenamiento."""
import hashlib
import sqlite3
from contextlib import closing

import numpy as np
import pandas as pd
from openpyxl import load_workbook

from .fuentes import RAIZ, SEMILLA

# Pares propuestos en el plan; varias isoenzimas pueden catalizar la misma reacción.
PARES = [("eno", "pep", "ENO", "P0A6P9"), ("pykF", "pep", "PYK", "P0AD61"),
         ("pykF", "pyr", "PYK", "P0AD61"), ("talB", "s7p", "TALA", "P0A870"),
         ("fumA", "fum", "FUM", "P0AC33"), ("sucC", "succ", "SUCOAS", "P0A836")]


def leer_proteinas(ish, destino):
    detalle = destino / "detalle"
    detalle.mkdir(parents=True, exist_ok=True)
    fuentes = RAIZ / "sources"
    with closing(sqlite3.connect(f"file:{fuentes / 'ecoli_auditada_v3.db'}?mode=ro", uri=True)) as db:
        t = pd.read_sql_query("""
            SELECT r.sample_id AS muestra, r.protein_series_id AS serie, m.run_id AS corrida,
                   m.protein_id_original AS proteina, m.concentration_original AS abundancia,
                   m.cv_percent AS CV_tecnico_pct, m.csv_row_original AS fila, m.csv_conc_column,
                   m.csv_cv_column, m.concentration_unit AS unidad
            FROM ishii_protein_measurements m JOIN ishii_sample_runs r USING(run_id)
            ORDER BY m.csv_row_original, r.sample_id, r.protein_series_id""", db)
        anotaciones = pd.read_sql_query("""
            SELECT protein_id_original, functional_annotation_reaction_id, relationship_type,
                   final_status, gene_identity_validation FROM ishii_mmc4_mapping""", db)
    libro = load_workbook(fuentes / "ishii.xlsx", data_only=True)
    hoja = libro["proteins"]
    assert t.unidad.eq(libro["generals"]["B9"].value).all()
    for c, valor, campo in [("csv_conc_column","abundancia","Conc"), ("csv_cv_column","CV_tecnico_pct","CV")]:
        coord = t[c] + t.fila.astype(str)
        t["celda_"+campo] = coord
        desde_excel = pd.to_numeric([hoja[a].value for a in coord], errors="raise")
        np.testing.assert_allclose(desde_excel, t[valor], atol=1e-12, rtol=0, equal_nan=True)
        assert all(hoja[s+"4"].value == campo for s in t[c].unique())
    assert all(hoja[f"A{r}"].value == g for r,g in t[["fila","proteina"]].drop_duplicates().itertuples(index=False))
    for f in t[["csv_conc_column","muestra","serie"]].drop_duplicates().itertuples(index=False):
        assert hoja[f.csv_conc_column+"1"].value == f.muestra and hoja[f.csv_conc_column+"2"].value == f.serie
    # KO05x tiene cabecera, pero ni abundancias ni CV. RF04–06 no tienen columna.
    assert hoja["J1"].value == "KO05x" and all(hoja[f"{c}{r}"].value is None
                                               for c in ["J","K"] for r in t.fila.unique())
    libro.close()
    t = t[t.muestra.isin(ish)].copy()
    assert not t.duplicated(["corrida","proteina"]).any() and t.abundancia.dropna().ge(0).all()
    t.to_csv(detalle / "proteinas_por_corrida.csv.gz", index=False)
    datos = t.groupby(["muestra","proteina"]).abundancia
    resumen = datos.agg(["median","count","min","max","std"]).reset_index()
    resumen.to_csv(detalle / "proteinas_por_cultivo.csv.gz", index=False)
    proteinas = resumen.pivot(index="muestra",columns="proteina",values="median").reindex(ish)
    disponibilidad = proteinas.notna().sum().rename("cultivos_medidos").to_frame()
    disponibilidad["corridas_con_medicion"] = t.groupby("proteina").abundancia.count()
    disponibilidad["corridas_con_CV"] = t.groupby("proteina").CV_tecnico_pct.count()
    disponibilidad.to_csv(detalle / "inventario_proteico.csv")
    assert proteinas.shape == (30,67) and proteinas.notna().any().sum()==61
    assert proteinas.notna().any(axis=1).sum()==26 and t[t.muestra.eq("RF03")].corrida.nunique()==5
    mapa = load_workbook(fuentes / "mmc4.xlsx", read_only=True, data_only=True)
    reacciones = {fila[0]:(i,fila[2]) for i,fila in enumerate(mapa["reactions"].values,1)}
    red = []
    for gene, met, rxn, up in PARES:
        a = anotaciones[anotaciones.protein_id_original.eq(gene) &
                        anotaciones.functional_annotation_reaction_id.eq(rxn)]
        assert len(a)==1 and a.iloc[0].final_status=="functional_supported_S4"
        fila, formula = reacciones[rxn]
        assert met+"[c]" in formula.split() and proteinas[gene].notna().any()
        red.append(dict(proteina=gene,objetivo=met,reaccion=rxn,fila_mmc4=fila,formula=formula.strip(),
                        relacion=a.iloc[0].relationship_type,identidad=a.iloc[0].gene_identity_validation,
                        UniProt=up, fila_proteina=int(t.loc[t.proteina.eq(gene),"fila"].iloc[0]),
                        cultivos=int(proteinas[gene].count())))
    mapa.close()
    pd.DataFrame(red).to_csv(destino / "pares_de_red.csv", index=False)
    return proteinas


def elegir_sin_red(train, gene, objetivo):
    """Iguala el patrón de disponibilidad DENTRO del entrenamiento; el desempate ignora la red."""
    presente = train[gene].notna()
    candidatos = []
    for p in train.columns.drop(gene):
        validos = train[p].notna()
        if train.loc[presente & validos, p].nunique() < 2 or (presente & validos).sum() < 10:
            continue
        prioridad = hashlib.sha256(f"{SEMILLA}|{gene}|{objetivo}|{p}".encode()).hexdigest()
        candidatos.append((int((presente != validos).sum()), prioridad, p))
    if not candidatos:
        return None, None
    distancia, _, escogida = min(candidatos)
    return escogida, distancia
