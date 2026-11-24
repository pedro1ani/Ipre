"""Las consultas se construyen una vez por partición y sin usar sus valores de prueba."""
import hashlib
import json

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

from .datos import ANCLAS, SEMILLA


def bloques(ids, grupos, k):
    return [(ids[a], ids[b]) for a, b in GroupKFold(k, shuffle=True, random_state=SEMILLA)
            .split(ids, groups=grupos.loc[ids])]


def perfil(visibilidad, objetivo):
    faltan = [a for a in ANCLAS if a != objetivo and not visibilidad[a]]
    return "+".join(faltan) or "anclas_presentes"


def consultas(train, test, real=False):
    """Proyectar patrones en cuatro anclas prefijadas permite estudiar también patrones escasos."""
    observada = test.notna()
    consultas = []
    full = set(tuple(train.columns[~fila]) for fila in train.notna().to_numpy())
    for objetivo in test:
        anclas = [a for a in ANCLAS if a != objetivo]
        patrones = {tuple(a for a in anclas if pd.isna(f[a]))
                    for _, f in train.loc[train[objetivo].isna(), anclas].iterrows()}
        clases = [(o,"completa") for o in sorted(full) if objetivo in o and len(o) > 1]
        clases += [(tuple(sorted((objetivo, *p))),"conjunta_anclas") for p in sorted(patrones) if p]
        for cultivo, f in observada.iterrows():
            def anotar(tipo, ocultas):
                llave = (cultivo, objetivo, tipo, "|".join(ocultas))
                evento = hashlib.sha256(json.dumps(llave).encode()).hexdigest()[:20]
                visible = {a: bool(f[a] and a not in ocultas) for a in ANCLAS}
                consultas.append((evento,cultivo,objetivo,tipo,llave[-1],perfil(visible,objetivo)))
            if real:
                if not f[objetivo]:
                    anotar("real",(objetivo,))
                continue
            if not f[objetivo]:
                continue
            anotar("individual",(objetivo,))
            for ocultas, tipo in clases:
                if all(f[a] for a in ocultas):
                    # La combinación publicada debe quedar exactamente ausente en el panel definido.
                    universo = anclas if tipo == "conjunta_anclas" else test.columns
                    if all(f[a] or a in ocultas for a in universo):
                        anotar(tipo,ocultas)
    return pd.DataFrame(consultas, columns=["evento","cultivo","objetivo","tipo","ocultas","perfil"])


def escoger_companeros(train, eleg, y, maximo=3):
    """Asociación ordinal con penalización por escasos pares; exclusivamente en entrenamiento."""
    candidatos = [c for c in train if c != y and c in eleg]
    if not candidatos:
        return []
    z = train.loc[train[y].notna(), candidatos + [y]]
    pares = z[candidatos].count()
    # Ranking conjunto por columna; recalcular rangos en cada par (Spearman con datos ausentes).
    corr = z.corr(method="spearman", min_periods=5)[y].drop(y).abs()
    orden = (corr * pares / max(1,len(z))).dropna().sort_index().sort_values(
        ascending=False, kind="stable")
    return orden.head(maximo).index.tolist()


def elegibilidad(train, grupos, minimo=10):
    medidas = train.notna().astype(int)
    n = medidas.sum()
    ng = medidas.groupby(grupos.loc[train.index]).max().sum()
    basicas = (n >= minimo) & (ng >= 5) & (train.nunique() >= 2)
    pares = (medidas.T @ medidas).to_numpy(copy=True)
    np.fill_diagonal(pares, 0)  # Un metabolito no puede respaldarse a sí mismo.
    respaldo = (pares[:, basicas] >= 5).sum(axis=1)
    return pd.DataFrame({"n": n, "grupos": ng, "predictores_con_5_pares": respaldo,
                         "elegible": basicas & (respaldo >= 2)})
