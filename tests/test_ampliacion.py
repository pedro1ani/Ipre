"""Controles científicos del análisis. Ejecutar desde la raíz con unittest discover."""
import json
import unittest

import numpy as np
import pandas as pd

from ampliacion.modelos import CONFIG, predecir
from ampliacion.preparar import SALIDA, preparar
from ampliacion.validar import ARCHIVO_MASCARAS, bloques


class ValidacionMetabolica(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ish, cls.mmc, cls.wt, cls.meta, cls.ref = preparar()
        cls.mascaras = json.loads(ARCHIVO_MASCARAS.read_text())

    def test_no_hay_grupos_compartidos_ni_verdades_externas_al_entrenamiento_interno(self):
        vistos = set()
        for b, train, q in bloques(self.mascaras["externas"], self.ish, self.meta.grupo, self.ref.index, 5):
            self.assertFalse(vistos & set(q.celda))
            vistos.update(q.celda)
            for _, sub, qi in bloques(self.mascaras["internas"][str(b)], train, self.meta.grupo, self.ref.index, 3):
                self.assertFalse(set(q.muestra) & (set(qi.muestra) | set(sub.index)))
                self.assertFalse(set(self.meta.loc[qi.muestra, "grupo"]) & set(self.meta.loc[sub.index, "grupo"]))
        self.assertEqual(len(vistos), 457)

    def test_ningun_metodo_puede_leer_el_valor_oculto(self):
        train = pd.DataFrame({"x": [1, 2, 3, 4, 5, 6], "y": [3, 5, 7, 9, 11, 13],
                              "z": [1, 4, 9, 16, 25, 36]}, dtype=float)
        cultivo = pd.DataFrame({"x": [1.], "y": [6.], "z": [8.]}, index=["prueba"])
        q = pd.DataFrame({"muestra": ["prueba"], "metabolito": ["x"]})
        modificado = cultivo.assign(x=1000000.)
        for metodo in CONFIG:
            a, _ = predecir(train, cultivo, q, train.columns, metodo, 0)
            b, _ = predecir(train, modificado, q, train.columns, metodo, 0)
            self.assertTrue(a.pred_mM.notna().all())
            np.testing.assert_array_equal(a.pred_mM, b.pred_mM)

    def test_escalas_y_historico_sobre_las_mismas_celdas(self):
        r = pd.read_csv(SALIDA / "equivalencia_predicciones.csv")
        self.assertTrue(r.n_pares.eq(416).all())
        self.assertLess(r[["max_delta_escala_mM", "max_delta_archivo_mM"]].to_numpy().max(), 1e-8)
        self.assertFalse(r[["cambios_cobertura_escala", "cambios_cobertura_archivo"]].to_numpy().any())
        self.assertLess(pd.read_csv(SALIDA / "equivalencia_Z.csv").max_delta_Z.max(), 1e-12)

    def test_datos_originales_y_sensibilidad_WT(self):
        solo_ish = self.ish.columns.difference(self.ref.index)
        self.assertTrue(self.mmc[solo_ish].isna().all().all())
        reg = pd.read_csv(SALIDA / "ajustes.csv")
        n = reg[reg.diseno.isin(["principal", "WT_extra"])].pivot(
            index=["bloque", "escenario", "metodo"], columns="diseno", values="n_train")
        for b, _, q in bloques(self.mascaras["externas"], self.ish, self.meta.grupo, self.ref.index, 5):
            esperado = 0 if "WT" in set(self.meta.loc[q.muestra, "grupo"]) else len(self.wt)
            self.assertTrue((n.loc[b, "WT_extra"] - n.loc[b, "principal"]).eq(esperado).all())
        original = self.ish[self.ref.index]
        for esc in ["A", "B"]:
            for metodo in CONFIG:
                aplicadas = pd.read_csv(SALIDA / f"matriz_{esc}_{metodo}_mM.csv", index_col=0)
                pd.testing.assert_frame_equal(aplicadas.where(original.notna()), original,
                                              check_names=False, rtol=1e-12, atol=1e-14)


if __name__ == "__main__":
    unittest.main()
