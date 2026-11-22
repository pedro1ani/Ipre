"""Cobertura, ancho y sensibilidad a semillas; los contextos repetidos pesan una sola celda."""
import numpy as np
import pandas as pd

CLAVES = ["panel", "tipo", "escenario", "serie", "sorteos", "nivel"]
INDICADORES = ["cobertura", "ancho_mM", "sesgo_mM", "MAE_media_mM", "bajo", "sobre",
               "SD_condicional_mM", "SD_multiple_mM", "n_predictores", "n_objetivo_train"]


def estabilidad_montecarlo(tabla):
    id_caso = ["panel", "tipo", "escenario", "nivel", "bloque", "celda", "evento"]
    parametros = [("30_a_100", "prefijo", 30, "prefijo", 100),
                  ("100_a_200", "prefijo", 100, "prefijo", 200),
                  ("100_vs_otros_100", "prefijo", 100, "segunda_mitad", 100)]
    resultados = []
    for nombre, s1, n1, s2, n2 in parametros:
        def casos(serie, n):
            return tabla[tabla.serie.eq(serie) & tabla.sorteos.eq(n) & tabla.media_mM.notna()].set_index(id_caso)
        a, b = casos(s1, n1), casos(s2, n2)
        if not a.index.equals(b.index):
            b = b.loc[a.index]
        cambio = pd.DataFrame({"cambio_cobertura": b.cubierto.astype(float) - a.cubierto.astype(float),
                               "cambio_de_estado": b.cubierto.ne(a.cubierto),
                               "desplazamiento_inferior_mM": (b.inferior_mM-a.inferior_mM).abs(),
                               "desplazamiento_superior_mM": (b.superior_mM-a.superior_mM).abs()})
        c = cambio.groupby(level=id_caso[:4]+["celda"]).mean()
        g = c.groupby(level=id_caso[:4])
        res = g.mean().assign(celdas=g.size()).reset_index().assign(contraste=nombre)
        resultados.append(res)
    return pd.concat(resultados, ignore_index=True)


def resumir(tabla, destino):
    detalle = destino / "detalle"
    id_caso = ["tipo", "serie", "sorteos", "nivel", "bloque", "evento"]
    amplia = tabla.pivot(index=id_caso, columns=["panel", "escenario"], values="media_mM")
    if len(tabla) != amplia.size or amplia.shape[1] != 4:
        raise ValueError("Faltan intentos en un panel o escenario.")
    comunes = amplia.notna().all(axis=1)
    emparejado = comunes.reindex(pd.MultiIndex.from_frame(tabla[id_caso])).to_numpy()
    base = tabla.assign(cobertura=tabla.cubierto.astype(float), bajo=tabla.bajo.astype(float),
                        sobre=tabla.sobre.astype(float), sesgo_mM=tabla.error_mM,
                        MAE_media_mM=tabla.error_mM.abs(), intervalo=tabla.media_mM.notna())
    partes = []
    for alcance, sel in [("disponibles", base.intervalo), ("comunes", emparejado)]:
        t = base.copy()
        t.loc[~sel, INDICADORES] = np.nan
        t["evento_cubierto"] = t.cobertura.fillna(0)
        t["intervalo"] = sel
        e = t.groupby(CLAVES + ["grupo", "metabolito", "celda"], sort=False)
        c = e[INDICADORES].mean().assign(intentos=e.size(), intervalos=e.intervalo.sum(),
                                         aciertos_evento=e.evento_cubierto.sum()).reset_index()
        partes.append(c.assign(alcance=alcance, grupo_valido=c.grupo.where(c.intervalos.gt(0))))
    celdas = pd.concat(partes, ignore_index=True)
    claves = CLAVES + ["alcance"]
    def agregar(sub):
        g = celdas.groupby(claves+sub)
        s = g[INDICADORES].mean().assign(intentos=g.intentos.sum(), eventos_con_intervalo=g.intervalos.sum(),
            aciertos_evento=g.aciertos_evento.sum(), celdas_intentadas=g.celda.size(),
            celdas_con_intervalo=g.cobertura.count(), grupos_evaluados=g.grupo_valido.nunique(),
            mediana_ancho_mM=g.ancho_mM.median()).reset_index()
        s["disponibilidad_eventos"] = s.eventos_con_intervalo/s.intentos
        return s
    met = agregar(["metabolito"])
    globales = agregar([])
    macro = met.groupby(claves).cobertura.mean().rename("cobertura_macro")
    globales.merge(macro, on=claves, validate="one_to_one").to_csv(destino / "intervalos_global.csv", index=False)
    met.to_csv(detalle / "intervalos_por_metabolito.csv.gz", index=False)
    agregar(["grupo"]).to_csv(detalle / "intervalos_por_grupo.csv.gz", index=False)
    met[met.serie.eq("prefijo") & met.sorteos.eq(200) & met.alcance.eq("comunes")].to_csv(
        destino / "intervalos_por_metabolito_200.csv", index=False)
    estabilidad_montecarlo(tabla).to_csv(destino / "estabilidad_sorteos.csv", index=False)
    ajustes = tabla[["bloque", "panel", "escenario", "k", "rondas_deterministas",
                    "aviso_determinista", "rondas_posteriores"]].drop_duplicates()
    assert len(ajustes) == 20
    ajustes.to_csv(destino / "convergencia.csv", index=False)
