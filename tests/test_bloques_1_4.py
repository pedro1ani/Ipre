import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from ampliacion.mascaras_conjuntas import TIPOS, coleccion
from ampliacion.medir_bloques import medir
from ampliacion.modelos import CONFIG, predecir
from ampliacion.preparar import preparar


class PruebasDeBloques(unittest.TestCase):
    def test_todo_el_patron_queda_invisible_a_los_tres_metodos(self):
        rng = np.random.default_rng(831)
        train = pd.DataFrame(rng.uniform(.1, 3, size=(14, 5)), columns=list("abcde"))
        cultivo = pd.DataFrame([[.8, .9, 1, 2, 1.3]], index=["KO"], columns=train.columns)
        q = pd.DataFrame({"muestra":["KO","KO"], "metabolito":["a","b"], "ocultar":["a+b","a+b"]})
        otro = cultivo.assign(a=6000, b=55000)
        for metodo in CONFIG:
            p, _ = predecir(train, cultivo, q, train.columns, metodo, 0)
            cambio, _ = predecir(train, otro, q, train.columns, metodo, 0)
            self.assertTrue(p.n_predictores.eq(3).all() and p.pred_mM.notna().all())
            np.testing.assert_array_equal(p.pred_mM, cambio.pred_mM)
        self.assertEqual(cultivo.at["KO","a"], .8)

    def test_catalogo_exacto_sin_grupos_compartidos_y_verdades_medidas(self):
        with tempfile.TemporaryDirectory() as tmp:
            raiz = Path(tmp)
            ish, _, _, meta, ref = preparar(raiz)
            ruta = raiz/"mascaras.json"
            c = coleccion(ish, meta.grupo, ref.index, ruta)
            for b in c["bloques"]:
                for d in [b["externo"], *b["internos"]]:
                    train = ish.loc[d["entrenamiento"]]
                    prueba = ish.loc[d["prueba"]]
                    self.assertFalse(set(meta.loc[train.index,"grupo"]) & set(meta.loc[prueba.index,"grupo"]))
                    if d is not b["externo"]:
                        self.assertTrue(set(train.index)|set(prueba.index) <= set(b["externo"]["entrenamiento"]))
                    for tipo in TIPOS:
                        panel = list(ref.index) if tipo == "documental" else d["columnas"]
                        esperado = {tuple(row.index[row.isna()]) for _,row in train[panel].iterrows()}
                        esperado = {m for m in esperado if len(m)>=2}
                        self.assertEqual({tuple(p["ocultar"]) for p in d[tipo]}, esperado)
                        for p in d[tipo]:
                            ms = p["ocultar"]
                            self.assertTrue(set(p["origen"]) <= set(train.index))
                            self.assertTrue(prueba.loc[p["aplicar"],ms].notna().all().all())
                            self.assertEqual(set(p["aplicar"]), set(prueba.index[prueba[ms].notna().all(axis=1)]))
            self.assertEqual(coleccion(ish,meta.grupo,ref.index,ruta), c)
            ruta.write_text(json.dumps({"otro":"estudio"}))
            with self.assertRaisesRegex(ValueError, "incompatibles"):
                coleccion(ish,meta.grupo,ref.index,ruta)

    def test_celda_repetida_y_abstencion_no_mejoran_el_promedio(self):
        rows = []
        for metodo in ["M","K"]:
            for celda, evento, estimado in [("c1","c1|a",2.),("c1","c1|b",4.),
                                             ("c2","c2|a",8.),("c3","c3|a",0.)]:
                if metodo == "K" and celda == "c3":
                    estimado = np.nan
                rows.append(dict(celda=celda,evento=evento,pred_mM=estimado,verdad_mM=0.,metodo=metodo,
                                 experimento="x", alternativa="fijo", escenario="A", rep=0, bloque=0,
                                 IQR_train=2., referencia_mM=1., grupo=celda, metabolito="m"))
        globales, _, _ = medir(pd.DataFrame(rows))
        comun = globales[globales.alcance.eq("comunes")]
        self.assertTrue(comun.celdas_predichas.eq(2).all() and comun.eventos_predichos.eq(3).all())
        self.assertTrue(comun.MAE_mM_macro.eq(5.5).all())
        self.assertTrue(comun.cobertura_eventos.eq(.75).all())
