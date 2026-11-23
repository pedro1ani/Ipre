import json
import unittest

import numpy as np
import pandas as pd

from ampliacion.fisiologia import REVISION
from ampliacion.modelos import predecir
from ampliacion.preparar import AUDITORIA, SALIDA
from ampliacion.validar import ARCHIVO_MASCARAS


class RevisionFisiologica(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # En un checkout limpio el detalle por cultivo se regenera desde los datos del proyecto.
        if not all((REVISION / n).is_file() for n in
                   ("condiciones.csv", "ajustes.csv", "predicciones.csv.gz", "global.csv", "pesos_Gerosa.csv")):
            from ampliacion.revisar_escalas import revisar
            revisar()
        cls.tasas = pd.read_csv(REVISION / "condiciones.csv", index_col="muestra")
        cls.ajustes = pd.read_csv(REVISION / "ajustes.csv")
        cls.pred = pd.read_csv(REVISION / "predicciones.csv.gz")
        cls.matriz = pd.read_csv(SALIDA / "matriz_mM.csv", index_col=0)
        cls.meta = pd.read_csv(AUDITORIA / "metadatos.csv", index_col=0)

    def test_cruce_de_hojas_y_unidades(self):
        t = self.tasas
        self.assertEqual(self.meta.at["RF05", "source_column"], "AK")
        self.assertEqual(t.at["RF05", "celda_q_sustrato"], "AI3")
        self.assertAlmostEqual(t.at["RF05", "q_sustrato"], 2.86)
        self.assertTrue(pd.isna(t.at["RF04", "q_sustrato"]))
        self.assertEqual(t.at["RF04", "v_rel"], 14)
        self.assertTrue(pd.isna(t.at["KO15", "v_rel"]))  # La reacción se excluyó del modelo.
        self.assertEqual(t.at["KO17", "v_rel"], 0.0)  # Este cero sí figura en la hoja.
        self.assertAlmostEqual(t.at["MMC2_GLUCOSE", "q_C"], 6 * 9.65)
        self.assertAlmostEqual(t.at["MMC2_ACETATE", "q_C"], 2 * 13.58)
        self.assertTrue(t[t.cohorte.eq("MMC2")].D.isna().all())

    def test_particion_y_datos_disponibles(self):
        mascaras = pd.DataFrame(json.loads(ARCHIVO_MASCARAS.read_text())["externas"])
        for fila in self.ajustes.itertuples():
            train = fila.ids_entrenamiento.split("|")
            prueba = mascaras[mascaras.bloque.eq(fila.bloque)].muestra.unique()
            self.assertTrue(set(self.meta.loc[train, "grupo"]).isdisjoint(self.meta.loc[prueba, "grupo"]))
            if fila.estrategia == "indice_q":
                self.assertAlmostEqual(fila.referencia_q_Ishii_train, self.tasas.loc[train, "q_sustrato"].median())
            if fila.familia == "flujo":
                self.assertTrue(self.tasas.loc[train, ["v_abs", "v_rel"]].notna().all().all())
            if fila.estrategia == "cov_D" and "WT" in set(self.meta.loc[prueba, "grupo"]):
                self.assertFalse(fila.covariable_activa)
        for _, filas in self.pred.groupby(["familia", "estrategia", "escenario", "metodo"]):
            self.assertEqual(set(filas.celda), set(mascaras.celda))
        cap = self.pred[self.pred.familia.eq("captacion")]
        sin_q = cap.muestra.map(self.tasas.q_C).isna()
        self.assertTrue(cap.loc[sin_q, "pred_mM"].isna().all())
        self.assertTrue(cap.loc[sin_q, "motivo"].eq("sin_fisiologia_requerida").all())

    def test_retorno_del_indice_a_mM(self):
        # Para Mediana, la predicción física debe ser q_test × mediana(C_train/q_train).
        f = self.pred.query("estrategia=='indice_q' and metodo=='Mediana' and metabolito=='atp' "
                            "and pred_mM==pred_mM").iloc[0]
        ajuste = self.ajustes.query("estrategia=='indice_q' and metodo=='Mediana' and bloque==@f.bloque")
        ids = ajuste.iloc[0].ids_entrenamiento.split("|")
        valor = (self.matriz.loc[ids, "atp"] / self.tasas.loc[ids, "q_sustrato"]).median()
        self.assertAlmostEqual(f.pred_mM, valor * self.tasas.at[f.muestra, "q_sustrato"], places=12)

    def test_la_fisiologia_no_sustituye_el_soporte_de_metabolitos(self):
        train = pd.DataFrame({"x": [1., 2, 3, 4], "y": [3., 4, 5, 6], "z": [2., 6, 3, 4]})
        cultivo = pd.DataFrame({"x": [50.], "y": [5.], "z": [np.nan]}, index=["test"])
        q = pd.DataFrame({"muestra": ["test"], "metabolito": ["x"]})
        cov = pd.Series([1., 2., 3., 4., 5.], index=[0, 1, 2, 3, "test"])
        for metodo in ["KNN", "Iterative"]:
            p, _ = predecir(train, cultivo, q, list(train), metodo, 0, covariable=cov)
            self.assertTrue(p.pred_mM.isna().all())
            self.assertEqual(p.at[0, "motivo"], "menos_de_dos_predictores")

    def test_mascara_global_y_peso_Gerosa(self):
        n = self.pred.pivot(index="celda", columns=["familia", "escenario", "metodo", "estrategia"],
                            values="pred_mM").notna().all(axis=1).sum()
        globales = pd.read_csv(REVISION / "global.csv")
        self.assertTrue(globales.query("alcance=='comunes_todas'").n_pares.eq(n).all())
        pesos = pd.read_csv(REVISION / "pesos_Gerosa.csv").dropna(subset=["MAE_mM"])
        np.testing.assert_allclose(pesos.MAE_R, pesos.MAE_mM / pesos.referencia_mM, atol=1e-12)
        grupos = pesos.groupby(["familia", "escenario", "metodo"])
        np.testing.assert_allclose(grupos.aporte_pct_MAE_R.sum(), 100)
        np.testing.assert_allclose(grupos.aporte_pct_MAE_mM.sum(), 100)


if __name__ == "__main__":
    unittest.main()
