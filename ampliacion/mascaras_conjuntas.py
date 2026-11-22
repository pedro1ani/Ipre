"""Patrones exactos observados en Ishii. No se inventan subpatrones para ganar casos."""
import hashlib
import json

import pandas as pd

from .modelos import CONFIG, elegibilidad, predecir
from .validar import ARCHIVO_MASCARAS, bloques

TIPOS = ("documental", "elegible")
CLAVES = ["muestra", "metabolito", "celda", "ocultar", "evento", "verdad_mM", "bloque"]


def catalogar(train, prueba, panel):
    """El entrenamiento propone; la prueba solo indica si existe la verdad completa."""
    panel = list(panel)
    patrones = {}
    for cultivo, fila in train[panel].isna().iterrows():
        faltan = tuple(fila.index[fila])
        if len(faltan) >= 2:
            patrones.setdefault(faltan, []).append(cultivo)
    return [
        {"ocultar": list(faltan), "origen": origen,
         "aplicar": prueba.index[prueba[list(faltan)].notna().all(axis=1)].tolist()}
        for faltan, origen in patrones.items()
    ]


def particion(train, prueba, grupos, objetivos):
    assert set(grupos.loc[train.index]).isdisjoint(grupos.loc[prueba.index])
    e = elegibilidad(train[objetivos], grupos)
    panel = e.index[e.elegible].tolist()
    return {"entrenamiento": train.index.tolist(), "prueba": prueba.index.tolist(), "columnas": panel,
            "documental": catalogar(train, prueba, objetivos),
            "elegible": catalogar(train, prueba, panel)}


def coleccion(ish, grupos, objetivos, ruta):
    archivo = json.loads(ARCHIVO_MASCARAS.read_text())
    tramos = []
    for b, train, q in bloques(archivo["externas"], ish, grupos, objetivos, 5):
        externo = particion(train, ish.loc[q.muestra.unique()], grupos, objetivos)
        internos = []
        for _, sub, qi in bloques(archivo["internas"][str(b)], train, grupos, objetivos, 3):
            internos.append(particion(sub, train.loc[qi.muestra.unique()], grupos, objetivos))
        tramos.append({"externo": externo, "internos": internos})
    d = {"regla": "v1; patrón exacto de al menos dos ausencias; una vez por cultivo de prueba",
         "sha256_matriz": hashlib.sha256(ish.to_csv().encode()).hexdigest(),
         "sha256_grupos": hashlib.sha256(grupos.loc[ish.index].to_csv().encode()).hexdigest(),
         "sha256_individual": hashlib.sha256(ARCHIVO_MASCARAS.read_bytes()).hexdigest(),
         "objetivos": list(objetivos), "bloques": tramos}
    contenido = json.dumps(d, ensure_ascii=False, indent=2) + "\n"
    if ruta.exists() and ruta.read_text(encoding="utf-8") != contenido:
        raise ValueError(f"Máscaras incompatibles: conserve la corrida anterior en otra carpeta: {ruta}")
    if not ruta.exists():
        ruta.write_text(contenido, encoding="utf-8")
    return json.loads(ruta.read_text(encoding="utf-8"))


def consultas(division, tipo, ish, bloque):
    filas = []
    for patron in division[tipo]:
        ocultar = "+".join(patron["ocultar"])
        for s in patron["aplicar"]:
            for m in patron["ocultar"]:
                celda = f"{s}|{m}"
                filas.append((s, m, celda, ocultar, f"{celda}|{ocultar}", ish.at[s, m], bloque))
    q = pd.DataFrame(filas, columns=CLAVES)
    assert q.evento.is_unique and q.verdad_mM.notna().all()
    return q


def elegir_conjunto(internos, tipo, ish, mmc, fuente, registro, bloque):
    """Selecciona las opciones dentro de los grupos de entrenamiento, nunca el método."""
    partes = []
    for i, division in enumerate(internos):
        q = consultas(division, tipo, ish, i)
        if q.empty:
            continue
        x = ish.loc[division["entrenamiento"]]
        if fuente == "B":
            x = pd.concat([x, mmc])
        for metodo in ("KNN", "Iterative"):
            for k in range(len(CONFIG[metodo])):
                p, _ = predecir(x, ish, q, division["columnas"], metodo, k)
                partes.append(p.assign(k=k, metodo=metodo, interno=i))
    elegidos = {"Mediana": 0}
    tabla = pd.concat(partes, ignore_index=True) if partes else pd.DataFrame()
    for metodo in ("KNN", "Iterative"):
        votos = {}
        if not tabla.empty:
            t = tabla[tabla.metodo.eq(metodo)]
            comun = t.pivot(index="evento", columns="k", values="pred_mM").notna().all(axis=1)
            t = t[t.evento.isin(comun.index[comun])].copy()
            for k, v in t.groupby("k"):
                error = (v.pred_mM - v.verdad_mM).abs()
                # Primero promediar los patrones de cada celda: sigue siendo una sola medición.
                mae = error.groupby([v.metabolito, v.celda]).mean().groupby(level=0).mean().mean()
                registro.append(dict(bloque=bloque, tipo=tipo, escenario=fuente, metodo=metodo,
                                     k=k, intentos=len(comun), eventos=len(v), celdas=v.celda.nunique(),
                                     internos=v.interno.nunique(), MAE_macro_mM=mae))
                if v.interno.nunique() >= 2:
                    votos[k] = mae
        elegidos[metodo] = min(votos, key=votos.get) if votos else None
        if not votos:
            registro.append(dict(bloque=bloque, tipo=tipo, escenario=fuente, metodo=metodo,
                                 k=None, nota="sin_dos_bloques_internos_con_pares"))
    return elegidos
