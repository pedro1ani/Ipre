import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from ampliacion.incertidumbre import intervalos, preparar_ajuste, simular_lote
from ampliacion.selector import aplicar, decidir


def ejemplos(errores):
    filas = []
    for i in range(12):
        f = i % 3
        for metodo, valores in errores.items():
            filas.append(dict(interno=f, grupo=f"g{i%6}", metabolito="a", celda=f"c{i}",
                              evento=f"e{i}", metodo=metodo, verdad_mM=2, pred_mM=2+valores[f]))
    return pd.DataFrame(filas)


class PruebasSelector(unittest.TestCase):
    def test_ganar_el_promedio_no_basta_ni_se_cuentan_patrones_como_muestras(self):
        p = ejemplos({"KNN":[1,1,3], "Iterative":[2,2,2], "Mediana":[3,3,3]})
        regla = decidir(p, ["a"]).iloc[0]
        self.assertEqual((regla.elegido, regla.victorias, regla.celdas), ("KNN", 2, 12))
        pocas = p[p.celda.isin([f"c{i}" for i in range(6)])]
        repetidas = pd.concat([pocas, pocas.assign(evento=pocas.evento+"|otro")])
        self.assertEqual(decidir(repetidas, ["a"]).iloc[0].razon, "soporte_interno_insuficiente")
        p = ejemplos({"KNN":[0,2,2], "Iterative":[5,1,1], "Mediana":[6,6,6]})
        self.assertEqual(decidir(p, ["a"]).iloc[0].razon, "ventaja_no_repetida")
        p = ejemplos({"KNN":[1,1,1], "Iterative":[1,1,1], "Mediana":[3,3,3]})
        self.assertEqual(decidir(p, ["a"]).iloc[0].razon, "empate")

    def test_eleccion_mediana_respaldo_y_abstencion_son_distintos(self):
        f = pd.DataFrame(dict(evento=list("abcd")*3, metabolito=list("abcd")*3,
                              metodo=np.repeat(["Mediana","KNN","Iterative"],4),
                              pred_mM=[1,2,3,np.nan, 7,np.nan,8,np.nan, 9,9,9,np.nan],
                              motivo="con_soporte"))
        r = pd.DataFrame(dict(metabolito=list("abcd"), elegido=["Mediana","KNN",None,None]))
        out = aplicar(f,r)
        self.assertEqual(out.politica.tolist(), ["eleccion_interna",
            "mediana_general_por_falta_de_soporte_externo","mediana_general_sin_eleccion","abstencion"])
        np.testing.assert_allclose(out.pred_mM, [1,2,3,np.nan], equal_nan=True)


class PruebasIncertidumbre(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(14)
        self.x = pd.DataFrame(rng.uniform(.15,2,size=(18,5)), columns=list("abcde"))
        self.x.loc[[1,4],"a"] = np.nan
        self.x.loc[[0,10],"c"] = np.nan
        self.ish = pd.DataFrame([[1.,2.,3.,4.,5.], [2.,3.,1.,3.,5.]],
                                columns=self.x.columns, index=["test","sin_soporte"])
        self.q = pd.DataFrame(dict(muestra=["test","test","sin_soporte"], metabolito=["a","b","e"],
                                   ocultar=["a+b","a+b","b+c+d+e"], verdad_mM=[1.,2.,5.]))
        self.grupos = pd.Series({"test":"gX", "sin_soporte":"gY"})

    def test_sorteos_sin_verdad_y_con_semillas_reproducibles(self):
        f, datos = preparar_ajuste(self.x,self.ish,self.q,list(self.x),self.grupos,0)
        otra = self.ish.copy()
        otra.loc["test", ["a","b"]] = [1000.,5000.]
        otra.loc["sin_soporte", ["b","e"]] = [3000.,9000.]
        f2, d2 = preparar_ajuste(self.x,otra,self.q.assign(verdad_mM=9999),list(self.x),self.grupos,0)
        self.assertTrue(f.SD_condicional_mM.iloc[:2].gt(0).all())
        self.assertTrue(np.isnan(f.determinista_mM.iat[2]))
        np.testing.assert_allclose(f.determinista_mM, f2.determinista_mM, equal_nan=True)
        with tempfile.TemporaryDirectory() as tmp:
            a,b = Path(tmp)/"a.npy",Path(tmp)/"b.npy"
            simular_lote((a,1100,1103,datos))
            simular_lote((b,1100,1103,d2))
            dib = np.load(a)
            np.testing.assert_array_equal(dib,np.load(b))
            self.assertEqual(dib.shape,(3,2))
            self.assertTrue((dib>=0).all() and (dib.std(axis=0)>0).all())
            t = intervalos(f,dib,3,"prueba")
            self.assertTrue(t[t.metabolito.eq("e")].cubierto.isna().all())
            u = t[t.metabolito.eq("a")].sort_values("nivel")
            self.assertTrue(u.ancho_mM.is_monotonic_increasing)
            self.assertAlmostEqual(u.error_mM.iat[0], dib[:,0].mean()-1.)


class PruebasDeResultados(unittest.TestCase):
    """Al reproducir los bloques, contrasta las salidas contra los controles ya archivados."""

    def test_controles_puntuales_conservados(self):
        from ampliacion.diseno_5_6 import DESTINO, PANELES
        fich = DESTINO / "detalle/selector_predicciones.csv.gz"
        if not fich.exists():
            self.skipTest("Requiere ejecutar el bloque 5.")
        pred = pd.read_csv(fich)
        casos = ["escenario", "bloque", "panel", "tipo", "evento"]
        knn = pred[pred.metodo.eq("KNN")].set_index(casos)
        elegidos = pred.query("metodo=='Selector' and eleccion=='KNN' and politica=='eleccion_interna'").set_index(casos)
        np.testing.assert_array_equal(elegidos.n_donantes, knn.loc[elegidos.index].n_donantes)
        control = pd.read_csv(DESTINO.parent / "predicciones.csv.gz")
        indices = ["escenario","metodo","bloque","celda"]
        for panel, estudio in PANELES.items():
            q = pred.query("panel==@panel and tipo=='individual' and metodo!='Selector'").set_index(indices)
            antes = control[control.diseno.eq(estudio)].set_index(indices)
            self.assertTrue(q.index.is_unique)
            self.assertTrue(q.index.sort_values().equals(antes.index.sort_values()))
            np.testing.assert_allclose(q.pred_mM, antes.loc[q.index].pred_mM,
                                       equal_nan=True, atol=1e-12, rtol=0)

    def test_intervalos_y_cobertura_desde_sorteos(self):
        from ampliacion.diseno_5_6 import DESTINO
        detalle = DESTINO / "detalle"
        if not (DESTINO / "muestreo.json").exists():
            self.skipTest("Requiere ejecutar el bloque 6.")
        import json
        d = pd.read_csv(detalle / "intervalos.csv.gz")
        d = d.query("serie=='prefijo' and sorteos==200").copy()
        intervalo = d.inferior_mM.notna()
        verdad = d.verdad_mM.between(d.inferior_mM, d.superior_mM)
        np.testing.assert_array_equal(d.loc[intervalo,"cubierto"].astype(bool), verdad[intervalo])
        self.assertTrue(d.loc[~intervalo, ["media_mM","superior_mM","cubierto"]].isna().all().all())
        self.assertTrue(d.loc[intervalo,"inferior_mM"].ge(0).all())
        ids = ["panel", "escenario", "tipo", "bloque", "evento"]
        bajos = d.pivot(index=ids,columns="nivel",values="inferior_mM").dropna()
        altos = d.pivot(index=ids,columns="nivel",values="superior_mM").dropna()
        self.assertTrue(((bajos[95]<=bajos[80]) & (bajos[80]<=bajos[50])).all())
        self.assertTrue(((altos[95]>=altos[80]) & (altos[80]>=altos[50])).all())
        determinista = d.query("nivel==95 and tipo=='individual'")[ids+["determinista_mM"]]
        sel = pd.read_csv(detalle / "selector_predicciones.csv.gz").query(
            "metodo=='Iterative' and tipo=='individual'")
        par = determinista.merge(sel[ids+["pred_mM"]],on=ids,validate="one_to_one")
        np.testing.assert_allclose(par.determinista_mM,par.pred_mM,equal_nan=True,atol=1e-12,rtol=0)
        informe = pd.read_csv(DESTINO / "intervalos_global.csv").query(
            "serie=='prefijo' and sorteos==200 and nivel==95 and alcance=='comunes'")
        for s in informe.itertuples():
            casos = d[(d.panel==s.panel) & (d.tipo==s.tipo) & (d.escenario==s.escenario) &
                      d.nivel.eq(95) & d.media_mM.notna()].copy()
            self.assertEqual(casos.celda.nunique(),s.celdas_con_intervalo)
            casos["acierto"] = casos.verdad_mM.between(casos.inferior_mM,casos.superior_mM)
            self.assertAlmostEqual(casos.groupby("celda").acierto.mean().mean(),s.cobertura)
        # Reconstruir un intervalo desde sus 200 sorteos originales, sin utilizar límites guardados.
        h = json.loads((DESTINO/"muestreo.json").read_text())["huella_corrida"]
        cache = detalle/"sorteos"/h
        a = np.concatenate([np.load(f,allow_pickle=False) for f in sorted(cache.glob("comun_A_0_*.npy"))])
        real = d.query("panel=='comun' and escenario=='A' and bloque==0 and nivel==95").dropna(subset=["media_mM"])
        self.assertEqual(a.shape,(200,len(real)))
        np.testing.assert_allclose(np.quantile(a,[.025,.975],axis=0).T,real[["inferior_mM","superior_mM"]])
