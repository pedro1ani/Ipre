"""Audita, reproduce y evalúa fuentes, tamaño y patrones. No imputa huecos reales."""
import hashlib
import json
import platform
from importlib.metadata import version

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from .mascaras_conjuntas import TIPOS, coleccion, consultas, elegir_conjunto
from .medir_bloques import medir
from .modelos import CONFIG, SEMILLA, elegibilidad, predecir
from .preparar import AUDITORIA, RAIZ, SALIDA, preparar
from .validar import ARCHIVO_MASCARAS, bloques, elegir

DESTINO = SALIDA / "bloques_1_4"
FIJOS = {"Mediana": 0, "KNN": 1, "Iterative": 0}  # 5 uniformes; 50 rondas, tol 0,001.
REPETICIONES = 8


def huella(archivo):
    return hashlib.sha256(archivo.read_bytes()).hexdigest()


def inventariar(ish, meta, ref, detalle):
    cat = pd.read_csv(AUDITORIA / "homologacion_maestra.csv")
    cat = cat.query("dataset_id == 'ISHII2007' and compartment == 'intracellular'").copy()
    med = pd.read_csv(AUDITORIA / "mediciones_auditadas.csv", low_memory=False)
    med = med[med.analyte_id.isin(cat.analyte_id) & med.sample_key.isin(ish.index)].copy()
    assert med.unit_original.eq("mM").all() and med.is_imputed.eq(0).all()
    assert not med.duplicated(["sample_key", "analyte_id"]).any()
    med["grupo_medido"] = med.sample_key.map(meta.grupo).where(med.valor_mM.notna())
    obs = med.groupby("analyte_id").agg(n_30=("valor_mM", "count"),
        n_grupos=("grupo_medido", "nunique"), distintos=("valor_mM", "nunique"))
    cat = cat.join(obs, on="analyte_id", validate="one_to_one")
    fichas = pd.read_csv(detalle / "diccionario.csv").set_index("fuente_Ishii")
    cat["variable_panel"] = cat.analyte_id.map(fichas.variable)
    cat["decision"] = np.select(
        [cat.variable_panel.isin(ref.index), cat.variable_panel.notna(), cat.composite.eq(1),
         cat.nombre_repetido_en_base.eq(True), cat.n_30.lt(10), cat.n_grupos.lt(5), cat.distintos.lt(2)],
        ["objetivo_compartido", "predictor_solo_Ishii", "senal_compuesta",
         "etiqueta_repetida", "pocas_observaciones", "pocos_grupos", "sin_variacion"],
        default="candidata_para_revision_de_identidad")
    assert len(cat) == 579 and cat.n_30.ge(10).sum() == 149
    cat[["analyte_id", "source_row", "original_label", "protocol", "compartment", "unidad_original",
         "estado_homologacion", "n_30", "n_grupos", "distintos", "variable_panel", "decision"]
        ].to_csv(detalle / "inventario_Ishii.csv", index=False)
    return cat.groupby(["decision", "protocol"]).size().rename("filas").reset_index()


def evaluar(x, ish, q, cols, metodo, k, escalas, grupos, ref, exp, alt, esc, registro,
            dividir=False, rep=0, bloque=None):
    """El IQR común describe Ishii entrenado; MMC2 nunca aporta una verdad de prueba."""
    if q.empty:
        return None
    assert set(grupos.loc[x.index]).isdisjoint(grupos.loc[q.muestra])
    if k is None:
        p = q.assign(pred_mM=np.nan, motivo="sin_soporte_de_seleccion_interna")
        info = {"rondas": None, "aviso_convergencia": False}
    else:
        p, info = predecir(x, ish, q, cols, metodo, k, referencia=ref if dividir else None)
    b = int(q.bloque.iloc[0]) if bloque is None else bloque
    registro.append(dict(experimento=exp, alternativa=alt, escenario=esc, metodo=metodo,
                          bloque=b, rep=rep, k=k, n_train=len(x), n_columnas=len(cols),
                          se_ajusto=k is not None and len(cols)>0, **info))
    return p.assign(experimento=exp, alternativa=alt, escenario=esc, metodo=metodo,
                    rep=rep, evento=p.evento if "evento" in p else p.celda,
                    grupo=p.muestra.map(grupos), IQR_train=p.metabolito.map(escalas),
                    referencia_mM=p.metabolito.map(ref))


def curva(train, ish, mmc, q, grupos, ref, b, registro, muestreos):
    salida = []
    disponibles = np.array(sorted(grupos.loc[train.index].unique()))
    tamanos = [("50pct", len(disponibles)//2), ("75pct", len(disponibles)*3//4),
               ("100pct", len(disponibles))]
    for rep in range(REPETICIONES):
        rng = np.random.default_rng(np.random.SeedSequence([SEMILLA, b, rep]))
        orden = rng.permutation(disponibles)
        for nombre, n in tamanos:
            sub = train[grupos.loc[train.index].isin(orden[:n])]
            eleg = elegibilidad(sub[ref.index], grupos)
            if nombre == "50pct":
                columnas = eleg.index[eleg.elegible]
                iqr = sub.quantile(.75) - sub.quantile(.25)
            assert eleg.loc[columnas, "elegible"].all()
            muestreos.append(dict(bloque=b, rep=rep, tamano=nombre, n_grupos=n, n_cultivos=len(sub),
                                  n_elegibles=eleg.elegible.sum(), panel_fijo=len(columnas),
                                  grupos="+".join(orden[:n]), columnas="+".join(columnas)))
            for metodo, k in FIJOS.items():
                salida.append(evaluar(sub, ish, q, columnas, metodo, k, iqr, grupos, ref,
                                     "curva", nombre, "A", registro, rep=rep))
            if nombre == "50pct":
                # El mismo número de filas añadidas: nunca partir un grupo para llegar a ocho.
                adicionales = []
                for g in orden[n:]:
                    ids = train.index[grupos.loc[train.index].eq(g)].tolist()
                    if len(adicionales) + len(ids) <= len(mmc):
                        adicionales.extend(ids)
                if len(adicionales) != len(mmc):
                    raise ValueError("No hay ocho cultivos adicionales en grupos completos.")
                assert set(adicionales).isdisjoint(sub.index)
                muestreos[-1]["ocho_Ishii"] = "+".join(adicionales)
                comparadores = [("base","A",sub), ("mas_8_Ishii","A",pd.concat([sub,train.loc[adicionales]])),
                                ("mas_8_MMC2","B",pd.concat([sub,mmc]))]
                for alt, esc, x in comparadores:
                    for metodo, k in FIJOS.items():
                        salida.append(evaluar(x, ish, q, columnas, metodo, k, iqr, grupos, ref,
                                              "ocho_filas", alt, esc, registro, rep=rep))
    return salida


def matriz_por_serie(ish, mmc, meta, ref, externas, registro, series):
    familia = meta.loc[ish.index, "serie_familia"]
    salida = []
    for serie in sorted(familia.dropna().unique()):
        prueba = familia.index[familia.eq(serie)]
        grupos_fuera = meta.loc[prueba, "grupo"]
        entrena = familia.notna() & familia.ne(serie) & ~meta.loc[ish.index, "grupo"].isin(grupos_fuera)
        train = ish.loc[entrena]
        q = externas[externas.muestra.isin(prueba)].copy()
        q["bloque"] = int(serie)
        e = elegibilidad(train[ref.index], meta.grupo)
        cols = e.index[e.elegible]
        assert set(familia.loc[train.index]).isdisjoint([serie])
        series.append(dict(serie=int(serie), n_prueba=len(prueba), n_train=len(train),
                           entrenamiento=train.index.tolist(), prueba=prueba.tolist()))
        iqr = train.quantile(.75) - train.quantile(.25)
        for esc, x in [("A", train), ("B", pd.concat([train, mmc]))]:
            for metodo, k in FIJOS.items():
                salida.append(evaluar(x, ish, q, cols, metodo, k, iqr, meta.grupo, ref,
                                     "serie", "dejar_serie_y_grupo_fuera", esc, registro))
    return salida


def comprobar_referencia(pred):
    ids = ["escenario", "metodo", "celda"]
    anterior = pd.read_csv(RAIZ / "evaluacion/validacion_celdas.csv").set_index(ids)
    reciente = pd.read_csv(SALIDA / "predicciones.csv.gz")
    reciente = reciente[reciente.diseno.eq("seleccion_mM")].set_index(ids)
    pares = []
    for fuente, archivo in [("referencia_historica", anterior), ("referencia_actual", reciente)]:
        nuevas = pred[pred.experimento.eq(fuente)].pivot(index=ids, columns="alternativa", values="pred_mM")
        assert nuevas.index.sort_values().equals(archivo.index.sort_values())
        nuevas["archivo"] = archivo.pred_mM
        for (esc, metodo), g in nuevas.groupby(level=ids[:2]):
            assert g.isna().eq(g.mM.isna(), axis=0).all().all()
            dif1, dif2 = (g.mM - g.archivo).abs().max(), (g.mM - g.ratio).abs().max()
            assert dif1 < 1e-10 and dif2 < 1e-10 and g.mM.count() == 416
            pares.append(dict(fuente=fuente, escenario=esc, metodo=metodo, intentos=len(g),
                              pares=g.mM.count(), max_delta_archivo_mM=dif1, max_delta_escala_mM=dif2))
    return pd.DataFrame(pares)


def ejecutar(destino=DESTINO):
    detalle = destino / "detalle"
    detalle.mkdir(parents=True, exist_ok=True)
    ish, mmc, wt, meta, ref = preparar(detalle)
    grupos = meta.grupo
    assert grupos.loc[ish.index].nunique() == 25
    inventariar(ish, meta, ref, detalle).to_csv(destino / "inventario_resumen.csv", index=False)
    archivo = json.loads(ARCHIVO_MASCARAS.read_text())
    h = pd.read_csv(RAIZ / "evaluacion/parametros_seleccionados.csv", dtype={"bloque_externo": str})
    historico = h.set_index(["bloque_externo", "escenario", "metodo"]).configuracion
    mascaras = coleccion(ish, grupos, ref.index, detalle / "mascaras_conjuntas.json")
    salidas, registro, muestreos, busqueda, ajuste_conjunto = [], [], [], [], []
    soporte, patrones, z, series = [], [], [], []
    for b, train, q in bloques(archivo["externas"], ish, grupos, ref.index, 5):
        print(f"Bloque {b+1}/5: referencia, fuentes, tamaño y patrones", flush=True)
        internos = archivo["internas"][str(b)]
        elegidos = elegir(train, internos, mmc, grupos, ref, False, b, busqueda)
        division = mascaras["bloques"][b]
        cols = division["externo"]["columnas"]
        iqr = train.quantile(.75) - train.quantile(.25)
        for paso, sub in [("externa", train)] + [(f"interna_{i}", s) for i, s, _ in
                                                  bloques(internos, train, grupos, ref.index, 3)]:
            for umbral in (8, 10, 12):
                for panel, nombres in [("comun", ref.index), ("ampliado", ish.columns)]:
                    e = elegibilidad(sub[nombres], grupos, umbral)
                    soporte.append(e.rename_axis("variable").reset_index().assign(
                        bloque=b, etapa=paso, panel=panel, umbral=umbral))
        for esc, x in [("A", train), ("B", pd.concat([train, mmc]))]:
            z1 = StandardScaler().fit_transform(x[cols])
            z2 = StandardScaler().fit_transform(x[cols].div(ref.loc[cols]))
            z.append(dict(bloque=b, escenario=esc, max_delta_Z=np.nanmax(np.abs(z1-z2))))
            for metodo in CONFIG:
                previo = 0 if metodo == "Mediana" else int(historico.loc[str(b), esc, metodo])
                nuevo = elegidos[esc, metodo, "mM"]
                for familia, k in [("referencia_historica", previo), ("referencia_actual", nuevo)]:
                    for alt in ("mM", "ratio"):
                        salidas.append(evaluar(x, ish, q, cols, metodo, k, iqr, grupos, ref,
                                               familia, alt, esc, registro, dividir=alt=="ratio"))
                for alt, k in [("opciones_A", elegidos["A", metodo, "mM"]), ("reajuste", nuevo)]:
                    salidas.append(evaluar(x, ish, q, cols, metodo, k, iqr, grupos, ref,
                                           "fuentes", alt, esc, registro))
        salidas.extend(curva(train, ish, mmc, q, grupos, ref, b, registro, muestreos))
        for etapa, d in [("externa", division["externo"])] + [
                (f"interna_{i}", s) for i, s in enumerate(division["internos"])]:
            for tipo in TIPOS:
                ps = d[tipo]
                patrones.append(dict(bloque=b, etapa=etapa, tipo=tipo, descubiertos=len(ps),
                    aplicables=sum(bool(p["aplicar"]) for p in ps),
                    ocasiones=sum(len(p["aplicar"]) for p in ps),
                    intentos=sum(len(p["ocultar"])*len(p["aplicar"]) for p in ps),
                    cultivos_prueba=len(d["prueba"]),
                    cultivos_con_patron=len({s for p in ps for s in p["aplicar"]})))
        for tipo in TIPOS:
            conj = consultas(division["externo"], tipo, ish, b)
            for esc, x in [("A", train), ("B", pd.concat([train, mmc]))]:
                seleccion = elegir_conjunto(division["internos"], tipo, ish, mmc, esc, ajuste_conjunto, b)
                if conj.empty:
                    continue
                for metodo in CONFIG:
                    k = elegidos[esc, metodo, "mM"]
                    solo = conj.drop(columns="ocultar")
                    salidas.append(evaluar(x, ish, solo, cols, metodo, k, iqr, grupos, ref,
                                           "patron_"+tipo, "individual", esc, registro))
                    fijo = evaluar(x, ish, conj, cols, metodo, k, iqr, grupos, ref,
                                   "patron_"+tipo, "conjunto", esc, registro)
                    salidas += [fijo, fijo.assign(experimento="ajuste_"+tipo, alternativa="fijo")]
                    salidas.append(evaluar(x, ish, conj, cols, metodo, seleccion[metodo], iqr, grupos, ref,
                                           "ajuste_"+tipo, "reajuste", esc, registro))
    salidas.extend(matriz_por_serie(ish, mmc, meta, ref, pd.DataFrame(archivo["externas"]), registro, series))
    pred = pd.concat(salidas, ignore_index=True)
    glob, met, por_grupo = medir(pred)
    comprobar_referencia(pred).to_csv(destino / "reproduccion.csv", index=False)
    for nombre, tabla in [("global", glob), ("por_metabolito", met), ("patrones_por_bloque", pd.DataFrame(patrones)),
                           ("equivalencia_Z", pd.DataFrame(z))]:
        tabla.to_csv(destino / f"{nombre}.csv", index=False)
    sp = pd.concat(soporte, ignore_index=True)
    agg = sp.groupby(["panel", "etapa", "bloque", "umbral"]).elegible.sum().rename("variables").reset_index()
    agg.to_csv(destino / "elegibilidad.csv", index=False)
    registro = pd.DataFrame(registro)
    registro.groupby(["experimento", "alternativa", "escenario", "metodo"]).agg(
        intentos_de_ajuste=("bloque","size"), ajustes=("se_ajusto","sum"),
        avisos_convergencia=("aviso_convergencia","sum")
    ).reset_index().to_csv(destino / "convergencia.csv", index=False)
    pred.groupby(["experimento","alternativa","escenario","metodo","motivo"]).size().rename(
        "eventos").to_csv(destino / "motivos.csv")
    for nombre, t in [("predicciones", pred), ("por_grupo_y_metabolito", por_grupo), ("ajustes", registro),
                      ("busqueda_individual", pd.DataFrame(busqueda)), ("busqueda_conjunta", pd.DataFrame(ajuste_conjunto)),
                      ("submuestras", pd.DataFrame(muestreos)), ("elegibilidad_variable", sp)]:
        t.to_csv(detalle / f"{nombre}.csv.gz", index=False, compression={"method":"gzip", "mtime":0})
    (detalle / "series.json").write_text(json.dumps(series, indent=2), encoding="utf-8")
    from .figuras_bloques import contrastes_y_figura
    contrastes_y_figura(pred, glob, met, por_grupo, destino)
    entradas = [RAIZ/"sources/ishii.xlsx", RAIZ/"sources/mmc2.xlsx", ARCHIVO_MASCARAS,
                *[AUDITORIA/f for f in ("panel_documental_mM.csv","metadatos.csv","referencias_glucosa.csv",
                                       "homologacion_maestra.csv","mediciones_auditadas.csv")],
                RAIZ/"evaluacion/parametros_seleccionados.csv", RAIZ/"evaluacion/validacion_celdas.csv",
                SALIDA/"predicciones.csv.gz"]
    manifest = dict(base_auditada="f3d3b669b2d6b121649c9b23f97fec79decc89f5", python=platform.python_version(),
                    paquetes={p:version(p) for p in ("numpy","pandas","scipy","scikit-learn","openpyxl","threadpoolctl")},
                    semilla=SEMILLA, hilos_BLAS=1, repeticiones_curva=REPETICIONES, opciones=CONFIG,
                    opciones_fijas_curva_y_serie=FIJOS, n_Ishii=len(ish), grupos=grupos.loc[ish.index].nunique(),
                    WT_sensibilidad=len(wt), auxiliares=len(mmc), matriz_documental=[len(ish)+len(wt)+len(mmc),ish.shape[1]],
                    objetivos=len(ref), sha256_entradas={str(p.relative_to(RAIZ)):huella(p) for p in entradas},
                    sha256_mascaras_nuevas=huella(detalle/"mascaras_conjuntas.json"),
                    sha256_matriz_mM=huella(detalle/"matriz_mM.csv"))
    (destino/"ejecucion.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        ejecutar()
