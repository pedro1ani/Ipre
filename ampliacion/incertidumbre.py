"""Bloque 6: intervalos predictivos contra mediciones ocultas, siempre en mM."""
import argparse
import hashlib
import json
import os
import warnings
from concurrent.futures import ProcessPoolExecutor
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from .diseno_5_6 import DESTINO, PANELES, diseno, huella, panel_elegible
from .modelos import CONFIG, SEMILLA

SORTEOS = (30, 100, 200)
NIVELES = (50, 80, 95)
SEMILLA_INICIAL = 1100


def perfiles(ish, q, columnas):
    """Un perfil por cultivo y patrón: dos objetivos del mismo patrón comparten el sorteo."""
    perfil = q[["muestra", "ocultar"]].drop_duplicates().reset_index(drop=True)
    valores = ish.loc[perfil.muestra, columnas].to_numpy(dtype=float, copy=True)
    col = pd.Index(columnas)
    for i, patron in enumerate(perfil.ocultar):
        j = col.get_indexer(patron.split("+"))
        valores[i, j[j >= 0]] = np.nan
    indice = pd.MultiIndex.from_frame(perfil).get_indexer(
        pd.MultiIndex.from_frame(q[["muestra", "ocultar"]]))
    return valores, indice, col.get_indexer(q.metabolito)


def imputador(k, cero, semilla, posterior):
    return IterativeImputer(estimator=BayesianRidge(), initial_strategy="median",
                            imputation_order="ascending", skip_complete=False,
                            min_value=cero, random_state=semilla, sample_posterior=posterior,
                            **CONFIG["Iterative"][k])


def preparar_ajuste(x, ish, q, columnas, grupos, k):
    scaler = StandardScaler().fit(x[columnas].to_numpy())
    z = scaler.transform(x[columnas].to_numpy())
    valores, filas, cols = perfiles(ish, q, columnas)
    consultas = scaler.transform(valores)
    cero = -scaler.mean_ / scaler.scale_
    visibles = np.isfinite(valores[filas]).sum(axis=1)
    valido = (cols >= 0) & (visibles >= 2)
    ficha = q.copy().assign(grupo=q.muestra.map(grupos), n_predictores=visibles,
                            n_objetivo_train=q.metabolito.map(x.count()),
                            motivo=np.select([cols < 0, visibles < 2],
                                             ["objetivo_no_elegible", "menos_de_dos_predictores"],
                                             default="con_soporte"))
    with warnings.catch_warnings(record=True) as avisos:
        warnings.simplefilter("always", ConvergenceWarning)
        fijo = imputador(k, cero, SEMILLA, False).fit(z)
        completo = fijo.transform(consultas)
    ficha["rondas_deterministas"] = fijo.n_iter_
    ficha["aviso_determinista"] = any(issubclass(a.category, ConvergenceWarning) for a in avisos)
    ficha["rondas_posteriores"] = CONFIG["Iterative"][k]["max_iter"]
    ficha["determinista_mM"] = np.nan
    ficha.loc[valido, "determinista_mM"] = scaler.inverse_transform(completo)[filas[valido], cols[valido]]
    ficha["SD_condicional_mM"] = np.nan
    # Esta SD corresponde a UNA regresión con los otros metabolitos fijados en
    # el perfil determinista completado. Es previa al truncamiento en cero.
    ultimas = {paso.feat_idx: paso for paso in fijo.imputation_sequence_}
    for j in np.unique(cols[valido]):
        sel = valido & (cols == j)
        paso = ultimas[j]
        _, sd = paso.estimator.predict(completo[filas[sel]][:, paso.neighbor_feat_idx], return_std=True)
        ficha.loc[sel, "SD_condicional_mM"] = sd * scaler.scale_[j]
    datos = (z, consultas, cero, k, filas[valido], cols[valido], scaler.mean_, scaler.scale_)
    return ficha, datos


def simular_lote(tarea):
    """Cada semilla vuelve a ajustar el entrenamiento y completa los perfiles externos una vez."""
    archivo, inicio, fin, datos = tarea
    if archivo.exists():
        return
    z, consultas, cero, k, filas, cols, media, escala = datos
    dibujos = []
    with threadpool_limits(limits=1):
        for semilla in range(inicio, fin):
            modelo = imputador(k, cero, semilla, True).fit(z)
            est = modelo.transform(consultas)
            assert modelo.n_iter_ == CONFIG["Iterative"][k]["max_iter"]
            dibujos.append(est[filas, cols] * escala[cols] + media[cols])
    a = np.asarray(dibujos)
    if not np.isfinite(a).all() or (a < -1e-12).any():
        raise ValueError("El muestreo produjo un resultado no finito o negativo.")
    temporal = archivo.with_suffix(".tmp.npy")
    np.save(temporal, np.maximum(a, 0.0), allow_pickle=False)
    temporal.replace(archivo)


def intervalos(ficha, dibujos, n, etiqueta):
    """Los casos sin soporte permanecen en la tabla y en los denominadores de disponibilidad."""
    valido = ficha.motivo.eq("con_soporte")
    if dibujos.shape != (n, int(valido.sum())):
        raise ValueError("El número de sorteos o el orden de los casos es incompatible.")
    t = ficha.assign(sorteos=n, serie=etiqueta, media_mM=np.nan, SD_multiple_mM=np.nan)
    t.loc[valido, "media_mM"] = dibujos.mean(axis=0)
    t.loc[valido, "SD_multiple_mM"] = dibujos.std(axis=0, ddof=1)
    t["error_mM"] = t.media_mM - t.verdad_mM
    partes = []
    for nivel in NIVELES:
        p = (1 - nivel / 100) / 2
        lim = np.quantile(dibujos, [p, 1-p], axis=0)
        v = t.assign(nivel=nivel, inferior_mM=np.nan, superior_mM=np.nan)
        v.loc[valido, ["inferior_mM", "superior_mM"]] = lim.T
        v["ancho_mM"] = v.superior_mM - v.inferior_mM
        v["cubierto"] = v.verdad_mM.between(v.inferior_mM, v.superior_mM).where(valido)
        v["bajo"] = (v.verdad_mM < v.inferior_mM).where(valido)
        v["sobre"] = (v.verdad_mM > v.superior_mM).where(valido)
        partes.append(v)
    return pd.concat(partes, ignore_index=True)


def ejecutar(procesos=1, destino=DESTINO):
    ish, mmc, grupos, ref, opciones, tramos = diseno(destino)
    detalle = destino / "detalle"
    # Los lotes se pueden reanudar; otra matriz, código, entorno u opciones nunca reutiliza sorteos.
    archivos = [Path(__file__), Path(__file__).with_name("diseno_5_6.py"),
                Path(__file__).with_name("modelos.py"), detalle / "mascaras_conjuntas.json",
                detalle / "matriz_mM.csv", DESTINO.parent / "ajustes.csv"]
    paquetes = {p: version(p) for p in ["numpy", "scipy", "pandas", "scikit-learn", "threadpoolctl"]}
    id_corrida = hashlib.sha256(json.dumps([list(map(huella, archivos)), paquetes,
                                           SEMILLA_INICIAL, SORTEOS]).encode()).hexdigest()
    cache = detalle / "sorteos" / id_corrida
    cache.mkdir(parents=True, exist_ok=True)
    ajustes, tareas = {}, []
    for b, ((train, q), _) in enumerate(tramos):
        for panel, estudio in PANELES.items():
            cols = panel_elegible(train, grupos, ref, panel)
            for esc, x in [("A", train), ("B", pd.concat([train, mmc]))]:
                k = int(opciones.loc[estudio, b, esc, "Iterative"])
                ficha, datos = preparar_ajuste(x, ish, q, cols, grupos, k)
                clave = f"{panel}_{esc}_{b}"
                rutas = [cache / f"{clave}_{i:03}.npy" for i in range(0, 200, 10)]
                ajustes[clave] = ficha.assign(panel=panel, escenario=esc, k=k), rutas
                tareas.extend((ruta, SEMILLA_INICIAL+i*10, SEMILLA_INICIAL+(i+1)*10, datos)
                              for i, ruta in enumerate(rutas))
    tareas.sort(key=lambda tarea: tarea[1])
    print(f"20 entrenamientos × 200 semillas; {procesos} procesos; lotes ya listos: "
          f"{sum(t[0].exists() for t in tareas)}/{len(tareas)}", flush=True)
    with ProcessPoolExecutor(max_workers=procesos) as pool:
        for i, _ in enumerate(pool.map(simular_lote, tareas), 1):
            if i % 25 == 0:
                print(f"{i}/{len(tareas)} lotes terminados", flush=True)
    partes = []
    for ficha, rutas in ajustes.values():
        sorteos = np.concatenate([np.load(p, allow_pickle=False) for p in rutas])
        for n in SORTEOS:
            partes.append(intervalos(ficha, sorteos[:n], n, "prefijo"))
        partes.append(intervalos(ficha, sorteos[100:], 100, "segunda_mitad"))
    salida = pd.concat(partes, ignore_index=True)
    salida.to_csv(detalle / "intervalos.csv.gz", index=False)
    from .calibracion import resumir
    resumir(salida, destino)
    (destino / "muestreo.json").write_text(json.dumps(dict(
        semilla_primera=SEMILLA_INICIAL, semilla_ultima=SEMILLA_INICIAL+199,
        sorteos=list(SORTEOS), niveles=list(NIVELES), reajustes_independientes=4000,
        rondas="max_iter de cada opción interna; sin parada temprana posterior",
        paquetes=paquetes, huella_corrida=id_corrida), indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--procesos", type=int, default=1, choices=range(1, (os.cpu_count() or 1)+1),
                        metavar="N", help="procesos en paralelo; no cambia los sorteos")
    with threadpool_limits(limits=1):
        ejecutar(parser.parse_args().procesos)
