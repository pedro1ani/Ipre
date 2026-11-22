"""Diseño común: las particiones, máscaras y opciones ya elegidas en la validación interna."""
import hashlib
import json

import pandas as pd

from .mascaras_conjuntas import TIPOS, coleccion, consultas
from .modelos import elegibilidad
from .preparar import SALIDA, preparar
from .validar import ARCHIVO_MASCARAS, bloques

DESTINO = SALIDA / "bloques_5_6"
PANELES = {"comun": "seleccion_mM", "ampliado": "principal"}
TIPOS_PRUEBA = ("individual", *TIPOS)


def huella(archivo):
    return hashlib.sha256(archivo.read_bytes()).hexdigest()


def diseno(destino=DESTINO):
    detalle = destino / "detalle"
    detalle.mkdir(parents=True, exist_ok=True)
    ish, mmc, _, meta, ref = preparar(detalle)
    anterior = json.loads((SALIDA / "bloques_1_4/ejecucion.json").read_text())
    ruta = detalle / "mascaras_conjuntas.json"
    if huella(ARCHIVO_MASCARAS) != anterior["sha256_entradas"]["evaluacion/mascaras.json"]:
        raise ValueError("La colección individual ya no coincide con el bloque 4.")
    mascaras = coleccion(ish, meta.grupo, ref.index, ruta)
    if huella(ruta) != anterior["sha256_mascaras_nuevas"]:
        raise ValueError("Cambió la matriz, su orden, los grupos o los patrones del bloque 4.")
    archivo = json.loads(ARCHIVO_MASCARAS.read_text())
    opciones = pd.read_csv(SALIDA / "ajustes.csv").set_index(
        ["diseno", "bloque", "escenario", "metodo"]).k
    tramos = []
    for b, train, q in bloques(archivo["externas"], ish, meta.grupo, ref.index, 5):
        div = mascaras["bloques"][b]
        externa = (train, preguntas(q, div["externo"], ish, b))
        internas = []
        for i, sub, qi in bloques(archivo["internas"][str(b)], train, meta.grupo, ref.index, 3):
            internas.append((sub, preguntas(qi, div["internos"][i], train, i)))
        tramos.append((externa, internas))
    return ish, mmc, meta.grupo, ref, opciones, tramos


def preguntas(individual, division, ish, bloque):
    partes = [individual.assign(ocultar=individual.metabolito, evento=individual.celda,
                                tipo="individual")]
    for tipo in TIPOS:
        q = consultas(division, tipo, ish, bloque)
        if not q.empty:
            partes.append(q.assign(tipo=tipo))
    return pd.concat(partes, ignore_index=True)


def panel_elegible(train, grupos, ref, panel):
    nombres = ref.index if panel == "comun" else train.columns
    e = elegibilidad(train[nombres], grupos)
    return e.index[e.elegible].tolist()
