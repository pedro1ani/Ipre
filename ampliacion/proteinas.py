"""Bloque 8: cotejo Excel–SQLite y elección de un control sin consultar el error."""
import hashlib
import sqlite3
from contextlib import closing

import numpy as np
import pandas as pd
from openpyxl import load_workbook

from .modelos import SEMILLA
from .preparar import RAIZ

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


def sensibilidad_rf03(destino):
    """Repite el contraste individual usando una sola corrida de RF03 a la vez."""
    from .bloques_7_8 import evaluar

    detalle = destino / "detalle"
    corridas = pd.read_csv(detalle / "proteinas_por_corrida.csv.gz", float_precision="round_trip")
    medianas = corridas.groupby(["muestra","proteina"]).abundancia.median().unstack()
    rf = corridas[corridas.muestra.eq("RF03")].pivot(index="serie", columns="proteina", values="abundancia")
    reglas = pd.read_csv(detalle / "proteinas_elegidas.csv")
    fuente = pd.read_csv(detalle / "predicciones.csv.gz", float_precision="round_trip")
    pred = fuente[fuente.familia.isin(reglas.par) & fuente.tipo.eq("individual")]
    ajustes = pd.read_csv(detalle / "ajustes.csv.gz")
    ajustes = ajustes[ajustes.familia.isin(reglas.par)].set_index(["familia","bloque","brazo","metodo"]).sort_index()
    assert ajustes.index.is_unique
    ish = pd.read_csv(detalle / "matriz_mM.csv",index_col=0,float_precision="round_trip")
    columnas_id = ["familia","metodo","bloque","grupo","celda"]
    originales = pred.pivot(index=columnas_id, columns="alternativa", values="pred_mM")
    verdaderas = pred.drop_duplicates(columnas_id).set_index(columnas_id).verdad_mM
    simulados = []
    for regla in reglas.itertuples():
        casos = pred[pred.familia.eq(regla.par) & pred.bloque.eq(regla.bloque)]
        q = casos[casos.metodo.eq("Mediana") & casos.alternativa.eq("control")].drop(columns=["pred_mM","motivo"])
        for metodo in ["KNN","Iterative"]:
            ajuste = ajustes.loc[regla.par,regla.bloque,"control",metodo]
            x = ish.loc[ajuste.ids_train.split("|")]
            cols = ajuste.variables.split("|")
            for serie, fila in rf.iterrows():
                dispon = ajuste.ids_con_dato.split("|")
                # Una comparación entre resúmenes de RF03 solo tiene sentido con idéntico soporte.
                if fila[[regla.red,regla.sin_red]].isna().any():
                    raise ValueError("Falta una abundancia de RF03: no mezclar subconjuntos distintos.")
                for brazo, proteina in [("red",regla.red),("sin_red",regla.sin_red)]:
                    cov = medianas[proteina].copy()
                    cov.at["RF03"] = fila[proteina]
                    p, _ = evaluar(x,dispon,q,ish,cols,cov,metodo,int(ajuste.k))
                    simulados.append(p[columnas_id+["pred_mM"]].assign(metodo=metodo,brazo=brazo,serie=serie))
    alt = pd.concat(simulados,ignore_index=True).pivot(
        index=columnas_id+["serie"],columns="brazo",values="pred_mM")
    pares = alt.join(originales.add_suffix("_mediana"),on=columnas_id).join(verdaderas,on=columnas_id).dropna()
    registros = []
    for (par,metodo,serie), v in pares.groupby(level=["familia","metodo","serie"]):
        for nom in ("sin_red","red"):
            def mae(c):
                return (v[c] - v.verdad_mM).abs().mean()
            registros.append(dict(par=par,metodo=metodo,serie_RF03=serie,brazo=nom,pares=len(v),
                                  MAE_control_mM=mae("control_mediana"),MAE_mediana_RF03_mM=mae(nom+"_mediana"),
                                  MAE_una_corrida_mM=mae(nom)))
    pd.DataFrame(registros).to_csv(detalle/"sensibilidad_RF03.csv",index=False)
