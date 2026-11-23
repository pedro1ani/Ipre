"""Controles de selección, independencia y concordancia con los archivos originales."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from ampliacion.bloques_7_8 import DESTINO
from ampliacion.posterior_proteinas import entradas
from ampliacion.preparar import RAIZ
from ampliacion.proteinas import elegir_sin_red, leer_proteinas


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

    def test_muestreo_no_ve_ni_objetivo_ni_otros_componentes_ocultos(self):
        x=pd.DataFrame(np.arange(1.,25.).reshape(6,4),columns=list("abcd"),index=list("123456"))
        consulta=pd.DataFrame([[1000.,2000.,3.,4.]],index=["prueba"],columns=x.columns)
        q=pd.DataFrame({"muestra":["prueba","prueba"],"ocultar":["a+b","a+b"],"metabolito":["a","b"]})
        cov=pd.Series([1.,2.,3.,4.,5.,6.,19.],index=[*x.index,"prueba"])
        a=entradas(x,consulta,q,list(x.columns),0,cov)
        consulta.loc["prueba",["a","b"]]=[-999.,999999.]
        b=entradas(x,consulta,q,list(x.columns),0,cov)
        for i in (0,1,2,4,5,6,7):
            np.testing.assert_array_equal(a[i],b[i])
        self.assertEqual(a[1].shape,(1,5))
        self.assertTrue(np.isnan(a[1][0,:2]).all())
        self.assertAlmostEqual(a[6][-1],cov.loc[x.index].mean())
        np.testing.assert_array_equal(a[4],[0,0])
        consulta.loc["prueba","c"]=np.nan
        with self.assertRaises(AssertionError):
            entradas(x,consulta,q,list(x.columns),0,cov)

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


class Resultados(unittest.TestCase):
    @unittest.skipUnless((DESTINO/"detalle/predicciones.csv.gz").is_file(),"Primero ejecute el análisis 7–8.")
    def test_igualdad_de_brazos_grupos_retirados_y_mediana(self):
        j=pd.read_csv(DESTINO/"detalle/ajustes.csv.gz")
        p=pd.read_csv(DESTINO/"detalle/predicciones.csv.gz")
        e=pd.read_csv(DESTINO/"detalle/proteinas_elegidas.csv")
        series=pd.read_csv(RAIZ/"auditoria/metadatos.csv",index_col=0).grupo
        for (f,b,m),t in j[j.familia.isin(e.par)].groupby(["familia","bloque","metodo"]):
            self.assertEqual(set(t.brazo),{"control","sin_red","red"})
            for k in ["ids_train","variables","ids_con_dato","k"]:
                self.assertEqual(t[k].nunique(),1)
            verdad=p[(p.familia==f)&(p.bloque==b)].grupo
            self.assertTrue(set(series.loc[t.ids_train.iloc[0].split("|")]).isdisjoint(verdad))
        med=p[p.metodo.eq("Mediana") & p.familia.isin(e.par)].pivot(
            index=["familia","tipo","evento"],columns="alternativa",values="pred_mM")
        for brazo in ["sin_red","red"]:
            np.testing.assert_array_equal(med.control,med[brazo])

    @unittest.skipUnless((DESTINO/"tasas_resultados.csv").is_file(),"Primero ejecute el análisis 7–8.")
    def test_reproduce_control_fisiologico_y_dilucion_constante(self):
        g=pd.read_csv(DESTINO/"tasas_resultados.csv")
        h=pd.read_csv(DESTINO.parent/"revision/resumen_comparaciones.csv")
        for (f,esc,met),t in g.query("tipo=='individual'").groupby(["familia","escenario","metodo"]):
            r=t.iloc[0]
            if f=="D_principal":
                self.assertAlmostEqual(r.control,r.D,places=12)
                continue
            old={"captacion":"captacion","flujo":"flujo","D_con_WT":"dilucion_WT"}[f]
            sub=h.query("alcance=='comunes_familia' and familia==@old and escenario==@esc and metodo==@met")
            nombres={"mM_Z":"control","cov_q_C":"q_C","cov_v_abs":"v_abs","cov_v_rel":"v_rel",
                     "solo_origen":"solo_origen","cov_D":"D"}
            for v in sub[sub.estrategia.isin(nombres)].itertuples():
                self.assertEqual(r.celdas_predichas,v.n_pares)
                self.assertAlmostEqual(r[nombres[v.estrategia]],v.MAE_mM_macro,delta=6e-7)

    @unittest.skipUnless((DESTINO/"detalle/intervalos.csv.gz").is_file(),"Primero ejecute los intervalos del bloque 8.")
    def test_intervalos_desde_sorteos_y_disponibilidad(self):
        detalle=DESTINO/"detalle"
        it=pd.read_csv(detalle/"intervalos.csv.gz").query("sorteos==100 and nivel==95").copy()
        p=pd.read_csv(detalle/"predicciones.csv.gz")
        pred=p[p.metodo.eq("Iterative") & p.alternativa.eq("red") & p.pred_mM.notna()]
        ids=["familia","tipo","evento","bloque"]
        self.assertTrue(it[it.alternativa.eq("red")].set_index(ids).sort_index().index.equals(
            pred.set_index(ids).sort_index().index))
        for lote,t in it.groupby("lote"):
            orig=np.concatenate([np.load(f,allow_pickle=False) for f in sorted((detalle/"sorteos").glob(lote+"_*.npy"))])
            self.assertEqual(orig.shape[0],100)
            self.assertTrue(np.isfinite(orig).all() and orig.min()>=0)
            sub=orig[:,t.columna_sorteo.astype(int)]
            esperado=np.quantile(sub,[.025,.975],axis=0)
            np.testing.assert_allclose(esperado, t[["inferior_mM","superior_mM"]].to_numpy().T,atol=1e-12,rtol=1e-12)
        resumen=pd.read_csv(DESTINO/"intervalos_proteinas.csv").query("nivel==95")
        for r in resumen.itertuples():
            t=it[(it.familia==r.familia)&(it.alternativa==r.alternativa)&it.tipo.eq("individual")]
            self.assertEqual(len(t),r.celdas)
            self.assertEqual(t.verdad_mM.between(t.inferior_mM,t.superior_mM).sum(),r.aciertos_evento)


if __name__=="__main__":
    unittest.main()
