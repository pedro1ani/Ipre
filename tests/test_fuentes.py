"""Controles de regresión de los insumos y del resumen científico publicado."""
import json
import unittest

import pandas as pd

from auditoria.verificar import verificar
from protocolo.datos import PUB


class Fuentes(unittest.TestCase):
    def test_concentraciones_y_ausencias_en_ambas_bases(self):
        self.assertEqual(verificar(), 20483)

    def test_balance_del_inventario_y_las_casillas(self):
        publicacion = json.loads((PUB / "ejecucion.json").read_text())
        inventario = pd.read_csv(PUB / "inventario.csv")
        seleccion = inventario[inventario.decision.eq("candidato_preliminar")]
        cultivos = publicacion["cultivos"]
        observadas = int(seleccion.n_30.sum())
        faltantes = cultivos * len(seleccion) - observadas
        self.assertEqual((cultivos, len(inventario), len(seleccion), observadas, faltantes),
                         (30, 579, 135, 3402, 648))
        self.assertEqual((observadas, faltantes), (publicacion["medidos"], publicacion["huecos"]))
        estados = publicacion["por_estado"]
        self.assertEqual(estados, {"validada": 221, "respaldo_general": 418, "pendiente": 9})
        self.assertEqual(sum(estados.values()), faltantes)
