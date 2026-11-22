"""Evalúa una tasa o una proteína por vez sobre mediciones ocultas; nunca completa huecos reales."""
import argparse
import json
import platform
from importlib.metadata import version

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from .diseno_5_6 import diseno, huella, panel_elegible
from .modelos import CONFIG, SEMILLA, predecir
from .preparar import RAIZ, SALIDA
from .proteinas import PARES, elegir_sin_red, leer_proteinas
from .validar import ARCHIVO_MASCARAS

DESTINO = SALIDA / "bloques_7_8"


def opciones_del_control(internas, ids, ish, grupos, ref):
    """La opción se ajusta al control restringido. Las etiquetas de prueba externa no aparecen aquí."""
    partes = []
    for i, (sub, preguntas) in enumerate(internas):
        x = sub[sub.index.isin(ids)]
        q = preguntas[preguntas.tipo.eq("individual") & preguntas.muestra.isin(ids)]
        if q.empty:
            continue
        assert set(grupos.loc[x.index]).isdisjoint(grupos.loc[q.muestra])
        columnas = panel_elegible(x, grupos, ref, "comun")
        for metodo in ("KNN","Iterative"):
            for k in range(len(CONFIG[metodo])):
                p, _ = predecir(x, ish, q, columnas, metodo, k)
                partes.append(p[["celda","metabolito"]].assign(
                    interno=i, metodo=metodo, k=k, error=(p.pred_mM - p.verdad_mM).abs()))
    tabla = pd.concat(partes, ignore_index=True) if partes else pd.DataFrame()
    opciones, scores = {"Mediana": 0}, []
    for metodo in ("KNN","Iterative"):
        opciones[metodo] = None
        if tabla.empty:
            continue
        ancho = tabla[tabla.metodo.eq(metodo)].pivot(
            index=["interno","metabolito","celda"],columns="k",values="error").dropna()
        n_fold = ancho.index.get_level_values("interno").nunique()
        if n_fold < 2 or ancho.empty:
            continue
        error = ancho.groupby(level="metabolito").mean().mean()
        opciones[metodo] = int(error.idxmin())
        for k,v in error.items():
            scores.append(dict(metodo=metodo,k=k,MAE_macro_interno_mM=v,pares=len(ancho),internos=n_fold))
    return opciones, scores


def evaluar(train, disponibles, q, ish, cols, cov, metodo, k):
    """Conserva todos los intentos, incluidos los cultivos que no tienen la covariable."""
    claves = ["tipo","evento"]
    sub = q[q.muestra.isin(disponibles)]
    info = dict(rondas=0, aviso_convergencia=False)
    if k is None:
        p = sub.assign(pred_mM=np.nan, motivo="sin_dos_bloques_internos")
    elif sub.empty:
        p = sub.assign(pred_mM=np.nan,motivo="sin_covariable")
    else:
        p, info = predecir(train, ish, sub, cols, metodo, k, covariable=cov)
    p = q.merge(p[claves+["pred_mM","motivo"]],on=claves,how="left",validate="one_to_one")
    p["motivo"] = p.motivo.fillna("sin_covariable_emparejada")
    return p, info


def tasas_documentadas(fis, destino):
    r = [
        ("captacion","Ishii","rates!3","mmol glucosa/(gPS·h)","q_C",
         "Ensayo enzimático de glucosa extracelular; tasa específica del cultivo, no concentración interna."),
        ("captacion","MMC2","Physiology!4","mmol sustrato/(gPS·h)","q_C",
         "Depleción por HPLC, biomasa y crecimiento; media por condición y SD publicadas."),
        ("flujo_rel","Ishii","fluxes!12","% de captación","v_rel",
         "Estimación G6P → 6PG a partir de marcaje de aminoácidos proteogénicos y restricciones del modelo."),
        ("flujo_abs","Ishii","final fluxes!12","mmol/(gPS·h)","v_abs",
         "Estimación, no medición directa: flujo relativo × captación de glucosa / 100."),
        ("dilucion","Ishii","generals!C20,C22; metabolites!fila 4","h⁻¹","D",
         "Consigna del quimiostato; constante en Ishii30. Los cuatro WT adicionales solo entran en entrenamiento."),
        ("crecimiento","MMC2","Physiology!2","h⁻¹","mu",
         "Pendiente de log OD en lote; no se homologó con la dilución de un quimiostato.")]
    filas = []
    for clave, origen, fuente, unidad, col, trazado in r:
        cohortes = ["MMC2"] if origen=="MMC2" else ["Ishii30","WT4"]
        for coh in cohortes:
            v = fis.loc[fis.cohorte.eq(coh),col]
            filas.append(dict(variable=clave,fuente=origen,hoja_fila=fuente,unidad_fuente=unidad,
                              cohorte=coh,cultivos=len(v),con_dato=v.count(),minimo_modelo=v.min(),maximo_modelo=v.max(),
                              unidad_modelo="mmol C/(gPS·h)" if col=="q_C" else unidad,procedencia=trazado,
                              disponibilidad_al_imputar={
                                  "q_C":"Si se midió consumo extracelular y biomasa del cultivo o condición.",
                                  "v_rel":"Solo si hay marcaje independiente y se resolvió el modelo de flujo.",
                                  "v_abs":"Requiere marcaje, ajuste del modelo y captación del cultivo.",
                                  "D":"Si existe la consigna de dilución; no aplica al lote.",
                                  "mu":"Solo como metadato de MMC2, no como sustituto de D."}[col]))
    pd.DataFrame(filas).to_csv(destino/"tasas_documentadas.csv",index=False)


def ejecutar(procesos=1, destino=DESTINO):
    from .posterior_proteinas import simular

    ish, mmc, grupos, ref, opciones, tramos = diseno(destino)
    detalle = destino / "detalle"
    proteinas = leer_proteinas(ish.index, destino)
    # Estas salidas de auditoría se reconstruyen del Excel si no están incluidas en el checkout.
    from .fisiologia import REVISION, auditar
    ficha = json.loads((REVISION/"ejecucion.json").read_text())
    for ruta in ["sources/ishii.xlsx","sources/mmc2.xlsx","sources/ecoli_auditada_v3.db"]:
        if huella(RAIZ/ruta) != ficha["sha256"][ruta]:
            raise ValueError(f"Fuente fisiológica distinta: {ruta}")
    matriz = pd.read_csv(detalle/"matriz_mM.csv",index_col=0)
    meta = pd.read_csv(RAIZ/"auditoria/metadatos.csv",index_col=0)
    wt = matriz.loc[matriz.index.difference(ish.index.union(mmc.index),sort=False)]
    fis = auditar(ish.index, wt.index, mmc.index, meta, salida=detalle/"fisiologia")
    assert len(wt)==4 and fis.loc[ish.index,"D"].nunique()==1
    tasas_documentadas(fis, destino)
    predicciones, ajustes, elecciones, historial, simulaciones = [], [], [], [], []
    for b, ((train, q), internas) in enumerate(tramos):
        print(f"Bloque externo {b+1}/5: tasas y proteómica", flush=True)
        # La escala de error siempre corresponde a Ishii medido en el entrenamiento común de cada contraste.
        def comparar(familia, esc, sub, x, disponibles, preguntas, brazos, ks, objetivo=None):
            cols = panel_elegible(sub, grupos, ref, "comun")
            iqr = sub.quantile(.75)-sub.quantile(.25)
            assert set(grupos.loc[x.index]).isdisjoint(grupos.loc[preguntas.muestra])
            for brazo, cov in brazos.items():
                for metodo,k in ks.items():
                    p, info = evaluar(x,disponibles,preguntas,ish,cols,cov,metodo,k)
                    p = p.assign(familia=familia,escenario=esc,alternativa=brazo,metodo=metodo,rep=0,
                                 grupo=p.muestra.map(grupos), IQR_train=p.metabolito.map(iqr),
                                 referencia_mM=p.metabolito.map(ref))
                    p["experimento"] = familia+"__"+esc+"__"+p.tipo
                    predicciones.append(p)
                    ajustes.append(dict(familia=familia,escenario=esc,brazo=brazo,metodo=metodo,bloque=b,
                        k=k,ids_train="|".join(x.index),ids_con_dato="|".join(disponibles),
                        variables="|".join(cols),covariable=cov.name if cov is not None else None,**info))
                    if objetivo and metodo=="Iterative" and p.pred_mM.notna().any():
                        simulaciones.append((familia,brazo,b,x,cols,k,cov,p[p.pred_mM.notna()].copy()))
        for familia in ["captacion","flujo","D_principal","D_con_WT"]:
            sub = train
            disponibles = ish.index
            brazos = {"control": None}
            if familia in ("captacion","flujo"):
                nombres = ["q_C"] if familia=="captacion" else ["v_rel","v_abs"]
                medidos = fis[nombres].notna().all(axis=1)
                sub = train.loc[medidos.loc[train.index]]
                disponibles = ish.index[medidos.loc[ish.index]]
                brazos.update({n:fis[n] for n in nombres})
            else:
                if familia=="D_con_WT" and "WT" in set(grupos.loc[train.index]):
                    sub = pd.concat([train,wt])
                brazos["D"] = fis.D
            for esc in (["A","B"] if familia=="captacion" else ["A"]):
                x = sub if esc=="A" else pd.concat([sub,mmc])
                variantes = dict(brazos)
                if esc=="B":
                    variantes["solo_origen"] = fis.cohorte.eq("MMC2").astype(float).rename("origen")
                ks = {m:int(opciones.loc["seleccion_mM",b,esc,m]) for m in CONFIG}
                comparar(familia,esc,sub,x,disponibles,q,variantes,ks)
        cache = {}
        for gene, objetivo, _, _ in PARES:
            familia = gene+"_"+objetivo
            ciega, distancia = elegir_sin_red(proteinas.loc[train.index],gene,objetivo)
            if ciega is None:
                raise ValueError(f"Revisar el control emparejado: {familia}, bloque {b}")
            medidas = proteinas[[gene,ciega]].notna().all(axis=1)
            disponibles = ish.index[medidas]
            sub = train.loc[train.index.isin(disponibles)]
            assert len(sub)>=10 and grupos.loc[sub.index].nunique()>=5
            assert (proteinas.loc[sub.index,[gene,ciega]].nunique()>1).all()
            key = tuple(sub.index)
            if key not in cache:
                cache[key] = opciones_del_control(internas,sub.index,ish,grupos,ref)
            ks, busqueda = cache[key]
            historial.extend(dict(par=familia,bloque=b,**r) for r in busqueda)
            preguntas = q[q.metabolito.eq(objetivo)]
            elecciones.append(dict(par=familia,objetivo=objetivo,red=gene,sin_red=ciega,bloque=b,
                diferencia_patron_train=distancia,n_train=len(sub),grupos_train=grupos.loc[sub.index].nunique(),
                muestras_train="|".join(sub.index),n_prueba_con_dato=preguntas[preguntas.muestra.isin(disponibles)].muestra.nunique(),
                KNN=ks["KNN"],Iterative=ks["Iterative"]))
            comparar(familia,"A",sub,sub,disponibles,preguntas,{"control":None,
                     "sin_red":proteinas[ciega],"red":proteinas[gene]},ks,objetivo)
    pred = pd.concat(predicciones,ignore_index=True)
    pred.to_csv(detalle/"predicciones.csv.gz",index=False)
    pd.DataFrame(ajustes).to_csv(detalle/"ajustes.csv.gz",index=False)
    pd.DataFrame(elecciones).to_csv(detalle/"proteinas_elegidas.csv",index=False)
    pd.DataFrame(historial).to_csv(detalle/"seleccion_interna.csv.gz",index=False)
    print("Iniciando intervalos con los cultivos y opciones de cada control",flush=True)
    intervalos, n_fit = simular(simulaciones, ish, destino, procesos)
    from .resumir_7_8 import resumir
    resumir(pred, intervalos, destino)
    archivos = [RAIZ/"sources"/n for n in ("ishii.xlsx","mmc2.xlsx","mmc4.xlsx","ecoli_auditada_v3.db",
                                           "Ishii_methods.pdf","Gerosa_methods.pdf")]
    archivos += [ARCHIVO_MASCARAS,SALIDA/"ajustes.csv"]
    registro = dict(base="0da3f454275a433e8b8737af3d94dc5b849e0713",python=platform.python_version(),
        paquetes={p:version(p) for p in ["numpy","pandas","scipy","scikit-learn","openpyxl","threadpoolctl"]},
        semilla=SEMILLA,semillas_posteriores=[1100,1199],ajustes_posteriores=n_fit,grupos_externos=25,
        fuentes_sha256={str(p.relative_to(RAIZ)):huella(p) for p in archivos},
        sha256_mascaras_conjuntas=huella(detalle/"mascaras_conjuntas.json"))
    (destino/"ejecucion.json").write_text(json.dumps(registro,indent=2,ensure_ascii=False)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--procesos",type=int,default=1,help="paralelismo de los sorteos; no cambia los resultados")
    args=parser.parse_args()
    if args.procesos<1:
        parser.error("Se requiere al menos un proceso.")
    with threadpool_limits(limits=1):
        ejecutar(args.procesos)
