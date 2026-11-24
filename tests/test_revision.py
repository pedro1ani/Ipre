import unittest

import numpy as np
import pandas as pd

from protocolo.comparacion import bootstrap
from protocolo.datos import leer_matriz
from protocolo.distribuciones import describir, matrices
from protocolo.estrategias import DET, RES, FIJOS, aplicar_estrategias, elegir_config


class Decisiones(unittest.TestCase):
    def test_bootstrap_agrupa_todas_las_mediciones_del_mismo_grupo(self):
        # Dos grupos extremos producen réplicas formadas a veces por un único grupo.
        grupos = ["a"] * 20 + ["b"] * 20
        d = [0.] * 20 + [4.] * 20
        self.assertEqual(bootstrap(d, grupos), (0., 4.))
        self.assertTrue(np.isnan(bootstrap(d, ["a"] * 40)[0]))

    def test_verdad_oculta_no_interviene_en_eleccion(self):
        q = pd.DataFrame({"evento":["x","v"],"objetivo":["a","a"],"perfil":["b","b"]})
        datos = [(ev,"A_amplio",m, 1. if m == "Mediana" or ev == "x" else np.nan)
                 for ev in q.evento for m in FIJOS]
        p = pd.DataFrame(datos,columns=["evento","opcion","metodo","pred"])
        original = q.assign(clave="A_amplio/Mediana",estado="respaldo_general",pred=1.)
        ev = {"g":(1, 20, 8, 3, [], [], []), "s":{}, "m":{}, "p":{}}
        base = aplicar_estrategias(q, p, original, ev, "Global")
        cambiado = aplicar_estrategias(q.assign(verdad=[1000.,0.]),p,original.assign(verdad=[3.,4.]),ev,"Global")
        pd.testing.assert_frame_equal(base,cambiado)
        z=base[base.estrategia.eq("Respaldo_KNN")].set_index("evento")
        self.assertEqual(z.at["x","metodo"],"KNN3")
        self.assertEqual(z.at["v","metodo"],"Mediana")
        self.assertEqual(z.at["v","estado"],"respaldo")
        self.assertFalse(base.estrategia.str.contains("Ideal").any())

    def test_muestra_insuficiente_y_denominadores_casi_cero(self):
        from protocolo.estrategias import CONFIG
        r=pd.DataFrame({m:[1.,2.] for m in CONFIG}).assign(grupo=["a","b"])
        self.assertEqual(elegir_config(r),"Global")
        x=pd.DataFrame({"cero":[0.]*10,"muy_poco":[1e-18]*10,"breve":[1.,2.,3.]+[np.nan]*7})
        e=describir(x)
        self.assertTrue(e.loc[["cero","muy_poco"],"CV_pct"].isna().all())
        self.assertTrue(np.isnan(e.at["breve","asimetria"]))


@unittest.skipUnless((DET/"validacion.csv.gz").exists(),"Requiere ejecutar protocolo.revision")
class Resultados(unittest.TestCase):
    def test_exceso_de_mae_es_aditivo(self):
        z=pd.read_csv(RES/"descomposicion_regla.csv").query("prueba=='perturbacion/individual' and dimension=='estado'")
        m=pd.read_csv(RES/"comparacion_estrategias.csv").query("prueba=='perturbacion' and tipo=='individual' and alcance=='comunes'").set_index("estrategia")
        self.assertEqual(z.n.sum(), m.at["KNN3","predichas"])
        self.assertAlmostEqual(z.aporte_al_MAE_mM.sum(),m.at["Original","MAE_micro_mM"]-m.at["KNN3","MAE_micro_mM"],places=8)
        self.assertLess(m.at["Ideal_no_aplicable","MAE_micro_mM"],m.at["KNN5","MAE_micro_mM"])

    def test_matrices_y_comunes(self):
        ish=leer_matriz("ishii")
        h=pd.read_csv(RES/"estimaciones_por_casilla.csv.gz",float_precision="round_trip")
        self.assertFalse(h.duplicated(["cultivo","objetivo","estrategia"]).any())
        for m,x in matrices(ish,h):
            archivo=pd.read_csv(DET/"matrices"/(m+".csv"),index_col=0,float_precision="round_trip")
            np.testing.assert_array_equal(x.to_numpy(),archivo.to_numpy())
        comunes = h.pivot(index=["cultivo","objetivo"],columns="estrategia",values="pred").notna().all(axis=1)
        for m,x in matrices(ish,h,solo_comunes=True):
            self.assertEqual(int((ish.isna() & x.notna()).sum().sum()),int(comunes.sum()))

    def test_distribucion_posterior_y_cobertura(self):
        from protocolo.datos import PRIV
        r=pd.read_csv(DET/"incertidumbre_externa.csv.gz").iloc[0]
        op=r.clave.split("/")[0]
        with np.load(PRIV/"posteriores"/r.fold/(r.objetivo+"__"+op+".npz")) as f:
            v=f["valores"][:,np.flatnonzero(f["eventos"]==r.evento)[0]]
        self.assertAlmostEqual(v.std(ddof=1),r.sd_posterior_mM)
        np.testing.assert_allclose(np.quantile(v,[.025,.975]),[r.lo_mM,r.hi_mM])
        c=pd.read_csv(RES/"calibracion_posterior.csv")
        t=c.query("prueba=='perturbacion' and tipo=='individual' and nominal==0.95").iloc[0]
        self.assertEqual((t.dentro,t.n),(704,809))

