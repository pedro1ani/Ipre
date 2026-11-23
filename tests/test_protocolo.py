"""Controles del protocolo que pueden ejecutarse con o sin resultados individuales."""
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from protocolo.aplicar import comprobar
from protocolo.datos import ANCLAS, PRIV, PUB, leer_matriz
from protocolo.diseno import consultas
from protocolo.modelos import estimar, matriz_consultas
from protocolo.regla import BASE, aprender, decidir
from protocolo.validar import particiones


class Diseno(unittest.TestCase):
    def test_una_mascara_oculta_todos_sus_componentes_y_no_depende_de_la_verdad(self):
        cols=["y",*ANCLAS]
        tr=pd.DataFrame([[np.nan,np.nan,3,4,5],[1,2,3,4,5]],columns=cols)
        te=pd.DataFrame([[6.,7.,8.,9.,10.]],index=["cultivo"],columns=cols)
        q=consultas(tr,te)
        q=q[(q.objetivo=="y")&(q.tipo=="conjunta_anclas")]
        self.assertEqual(len(q),1)
        x=matriz_consultas(te,q,cols)
        self.assertTrue(np.isnan(x[0,:2]).all())
        te.loc["cultivo",["y","fum"]]=[1e9,-1e9]
        np.testing.assert_equal(x,matriz_consultas(te,q,cols))
        te.at["cultivo","fum"]=np.nan
        self.assertFalse(consultas(tr,te).query("objetivo=='y'").tipo.eq("conjunta_anclas").any())

    def test_vecinos_requieren_solapamiento_y_no_una_sustitucion_automatica(self):
        tr=np.array([[1,1,2,3],[2,2,2,3],[4,3,2,3],[999,9,np.nan,np.nan]],float)
        q=np.array([[np.nan,3,2,np.nan],[np.nan,3,np.nan,np.nan]])
        modelos=estimar(tr,q,4,("Mediana","KNN3"))
        pred,mot,_,don,*_=modelos["KNN3"]
        self.assertAlmostEqual(pred[0],7/3)
        self.assertEqual(don[0],3)
        self.assertTrue(np.isnan(pred[1]))
        self.assertEqual(mot[1],"menos_de_2_acompanantes")
        self.assertEqual(modelos["Mediana"][0][1],3)

    def test_mmc2_o_auxiliar_no_puede_ganar_solo_por_cobertura(self):
        filas=[]
        for k in range(14):
            for op,met,v in [("A_amplio","Mediana",4.),("A_amplio","KNN3",1.),
                             ("A_amplio","KNN5",2.),("A_amplio","Iterativo",3.),
                             ("A_comun","KNN3",.2),("B_comun","KNN3",.5)]:
                filas.append(dict(objetivo="y",perfil="p",cultivo=str(k),grupo=str(k),fold=k%3,
                                  opcion=op,metodo=met,pred=v,verdad=0))
        r=aprender(pd.DataFrame(filas))
        self.assertEqual(r["y|p"]["opcion"],"A_amplio/KNN3")
        chico=aprender(pd.DataFrame(filas).query("fold==0"))
        self.assertEqual(chico["y|p"]["estado"],"respaldo_general")

    def test_si_falta_auxiliar_se_exige_base_validada(self):
        f=SimpleNamespace(objetivo="y",perfil="p")
        alts={BASE:SimpleNamespace(pred=1.)}
        r={"opcion":"A_amplio_qC/KNN3","estado":"validada","base":BASE}
        r["estado_base"]="respaldo_general"
        self.assertEqual(decidir(f,alts,{"y|p":r})[1],"pendiente")
        r["estado_base"]="validada"
        self.assertEqual(decidir(f,alts,{"y|p":r})[0],BASE)

    def test_divisiones_nunca_comparten_grupos_ni_usan_auxiliares_como_verdad(self):
        idx=pd.Index([f"x{i}" for i in range(30)])
        meta=pd.DataFrame({"grupo":[str(i//2) for i in range(30)],
                           "serie_familia":np.tile(range(5),6)},index=idx)
        splits,ext=particiones(pd.DataFrame(index=idx),meta)
        self.assertEqual(len(ext),10)
        for tr,te in splits.values():
            self.assertFalse(set(meta.loc[tr,"grupo"]) & set(meta.loc[te,"grupo"]))
        externos=[*sum([list(te) for n,(_,te) in ext.items() if n.startswith("E")],[])]
        self.assertCountEqual(externos,idx)


@unittest.skipUnless((PRIV/"resultado_30_completo.csv").exists(),"Requiere ejecutar el protocolo.")
class Salidas(unittest.TestCase):
    def test_preserva_cada_medicion_y_cada_columna_excluida(self):
        comprobar()
        orig=leer_matriz("original_30_completa")
        matriz=leer_matriz("resultado_30_completo")
        np.testing.assert_array_equal(orig.to_numpy()[orig.notna()],matriz.to_numpy()[orig.notna()])
        inv=pd.read_csv(PUB/"inventario.csv")
        excl=inv.loc[inv.decision.ne("candidato_preliminar"),"variable"]
        pd.testing.assert_frame_equal(orig[excl],matriz[excl])
        self.assertTrue({"dhap","nadh","nadph"}.issubset(excl))
        z=pd.read_csv(PRIV/"auditoria_huecos.csv")
        num=matriz.notna().sum().sum()-orig.notna().sum().sum()
        self.assertEqual(z.pred.count(),num)
        self.assertEqual(len(z),len(z.drop_duplicates(["cultivo","objetivo"])))
        k=z[z.metodo.str.startswith("KNN",na=False)]
        self.assertTrue((k.ids_donantes_KNN.str.split("|").str.len()==k.donantes).all())

    def test_controles_auxiliares_tienen_las_mismas_filas(self):
        a=pd.read_csv(PRIV/"particiones/E1/ajustes.csv")
        for f in a.itertuples():
            if f.agregar:
                nom="ctrl_"+f.opcion.replace("_ciega_","_red_")
                par=a[a.objetivo.eq(f.objetivo)&a.opcion.eq(nom)].iloc[0]
                self.assertEqual(par.ids,f.ids)
                self.assertEqual(par.columnas,f.columnas)

    @unittest.skipUnless((PRIV/"intervalos.csv.gz").exists(),"Requiere ejecutar los intervalos.")
    def test_intervalos_solo_iterativos_y_dos_denominadores(self):
        e=pd.read_csv(PRIV/"externas.csv.gz")
        i=pd.read_csv(PRIV/"intervalos.csv.gz")
        z=e.merge(i,on=["fold","objetivo","evento"],how="left",validate="one_to_one")
        np.testing.assert_array_equal(z.lo_mM.notna(),z.clave.str.endswith("/Iterativo",na=False))
        ok=z.dropna(subset=["lo_mM"])
        self.assertTrue((ok.lo_mM>=0).all() and (ok.hi_mM>=ok.lo_mM).all())
        p=pd.read_csv(PUB/"intervalos.csv").query("prueba=='perturbacion' and tipo=='individual'").iloc[0]
        casos=z[z.fold.str.startswith("E")&z.tipo.eq("individual")]
        n=casos.lo_mM.notna().sum()
        self.assertEqual(n,p.con_intervalo)
        self.assertEqual(p.candidatas,len(casos))
        self.assertAlmostEqual((casos.verdad.between(casos.lo_mM,casos.hi_mM)).sum()/n,
                               p.cobertura_entre_intervalos)


if __name__=="__main__":
    unittest.main()
