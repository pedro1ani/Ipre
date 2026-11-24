"""Controles de selección, independencia y concordancia con los archivos originales."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from protocolo.fuentes import RAIZ
from protocolo.proteinas import elegir_sin_red, leer_proteinas


class Reglas(unittest.TestCase):
    def test_control_ciego_desempata_sin_usar_abundancia_u_orden_de_columnas(self):
        base = pd.DataFrame({"red": range(16), "a": range(16), "b": range(16), "otra": range(16)},dtype=float)
        base.loc[0,["red","a","b"]]=np.nan
        antes=elegir_sin_red(base,"red","objetivo")
        self.assertIn(antes[0],["a","b"])
        self.assertEqual(antes[1],0)
        cambio=base.loc[:,list(reversed(base.columns))].copy()*901+7
        self.assertEqual(antes,elegir_sin_red(cambio,"red","objetivo"))
        # Una columna constante no permite estimar un efecto aunque tenga cobertura perfecta.
        cambio.loc[cambio["a"].notna(),"a"]=1
        self.assertEqual(elegir_sin_red(cambio,"red","objetivo")[0],"b")

    def test_corridas_y_unidades_coinciden_con_excel(self):
        meta=pd.read_csv(RAIZ/"auditoria/metadatos.csv",index_col=0)
        ids=meta.index[meta.incluida & meta.dataset_id.eq("ISHII2007")]
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            m=leer_proteinas(ids,p)
            t=pd.read_csv(p/"detalle/proteinas_por_corrida.csv.gz")
            self.assertEqual(set(m.index[m.isna().all(axis=1)]),{"RF04","RF05","RF06","KO05x"})
            r=t.query("muestra=='RF03' and proteina=='sucC'")
            self.assertEqual(set(r.serie),set(range(2,7)))
            self.assertAlmostEqual(m.at["RF03","sucC"],r.abundancia.median())
            self.assertEqual(m.index.tolist(),ids.tolist())
            self.assertTrue(r.CV_tecnico_pct.notna().all())
