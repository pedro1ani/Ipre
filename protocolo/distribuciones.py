"""Describe lo medido, lo estimado y lo que se recupera cuando la verdad sí se conoce."""
import numpy as np
import pandas as pd

from .comparacion import publicar
from .datos import PRIV, leer_matriz
from .estrategias import DET, RES, PRINCIPALES, guardar

ESTADISTICOS = ["media", "mediana", "SD", "varianza", "CV_pct", "Q1", "Q3", "IQR", "asimetria"]


def describir(x):
    n = x.count()
    media, sd = x.mean(), x.std(ddof=1)
    q1, q3 = x.quantile(.25), x.quantile(.75)
    tol = np.maximum(1e-12, x.abs().max() * 1e-6)
    return pd.DataFrame(dict(n=n, faltantes=len(x)-n, faltantes_pct=(1-n/len(x))*100,
        media=media, mediana=x.median(), SD=sd, varianza=x.var(ddof=1),
        CV_pct=100*sd/media.where(media.abs()>tol), minimo=x.min(), Q1=q1, Q3=q3,
        IQR=q3-q1, maximo=x.max(), asimetria=x.skew().where((n>=8)&(sd>tol)),
        extremos_Tukey=((x.lt(q1-1.5*(q3-q1))) | (x.gt(q3+1.5*(q3-q1)))).sum())).rename_axis("objetivo")


def matrices(ish, huecos=None, solo_comunes=False):
    if huecos is None:
        huecos = pd.read_csv(DET/"huecos.csv.gz", float_precision="round_trip")
    p = huecos[huecos.estrategia.isin(PRINCIPALES)].pivot(
        index=["cultivo","objetivo"], columns="estrategia", values="pred")
    if solo_comunes:
        p = p.dropna()
    for m in PRINCIPALES:
        vacios = p[m].unstack().reindex_like(ish)
        assert not (vacios.notna() & ish.notna()).any().any()
        x = ish.fillna(vacios)
        np.testing.assert_array_equal(x.to_numpy()[ish.notna()], ish.to_numpy()[ish.notna()])
        yield m, x


def relativo(despues, antes, escala):
    # Un denominador pequeño en relación con su escala convierte el porcentaje en ruido.
    umbral = np.maximum(1e-12, escala.abs()*1e-6)
    return 100*(despues-antes)/antes.abs().where(antes.abs()>umbral)


def correlaciones(x, tipo):
    n = x.notna().astype(int)
    pares = n.T @ n
    c = x.corr(method=tipo, min_periods=8)
    a, b = np.triu_indices(len(c), k=1)
    ix = pd.MultiIndex.from_arrays([c.columns[a],c.columns[b]],names=["a","b"])
    return pd.DataFrame({"r":c.to_numpy()[a,b], "pares":pares.to_numpy()[a,b]}, index=ix)


def origen_extremos(ish):
    original = describir(ish)
    flag = ish.lt(original.Q1-1.5*original.IQR) | ish.gt(original.Q3+1.5*original.IQR)
    ix = flag.stack()
    ix = ix[ix].index.set_names(["sample_key","variable"])
    r = pd.read_csv(PRIV/"trazabilidad.csv.gz",low_memory=False).set_index(["sample_key","variable"])
    t = r.loc[ix,["measurement_uid","source_sheet","source_cell","protocol","valor_mM","unit_original"]]
    np.testing.assert_allclose(t.valor_mM,[ish.at[s,y] for s,y in ix],rtol=1e-10)
    guardar(t.reset_index(),"extremos_trazados.csv.gz")
    return flag


def intervalo():
    i = pd.read_csv(PRIV/"intervalos.csv.gz")
    z = pd.read_csv(DET/"validacion.csv.gz")
    z = z[z.estrategia.eq("Original")].drop(columns="objetivo").merge(
        i,on=["fold","evento"],validate="one_to_one")
    z["prueba"] = np.where(z.fold.str.startswith("E"),"perturbacion","serie")
    original = pd.read_csv(PRIV/"externas.csv.gz",usecols=["fold","evento","visibles","clave"])
    z = z.drop(columns="predictores").merge(original,on=["fold","evento"],validate="one_to_one")
    z = z.rename(columns={"visibles":"predictores"})
    z["ancho"] = z.hi_mM-z.lo_mM
    z["error_abs"] = (z.pred-z.verdad).abs()
    z["ancho_IQR"] = z.ancho/z.IQR.where(z.IQR.gt(1e-12))
    z["error_IQR"] = z.error_abs/z.IQR.where(z.IQR.gt(1e-12))
    calibracion, asociaciones = [], []
    for (cv,tipo), g in z.groupby(["prueba","tipo"]):
        for nivel in (50,80,90,95):
            suf = "mM" if nivel==95 else str(nivel)
            dentro = g.verdad.between(g["lo_"+suf],g["hi_"+suf])
            calibracion.append(dict(prueba=cv,tipo=tipo,nominal=nivel/100,n=len(g),
                grupos=g.grupo.nunique(), dentro=int(dentro.sum()),observada=dentro.mean(),
                ancho_mediano_mM=(g["hi_"+suf]-g["lo_"+suf]).median()))
        for y,h in [("todos",g), *g.groupby("objetivo")]:
            def rho(a,b):
                return h[a].corr(h[b],method="spearman") if len(h)>=8 and min(h[a].nunique(),h[b].nunique())>1 else np.nan
            asociaciones.append(dict(prueba=cv,tipo=tipo,objetivo=y,n=len(h),grupos=h.grupo.nunique(),
                SD_posterior_mediana_mM=h.sd_posterior_mM.median(), ancho_mediano_mM=h.ancho.median(),
                rho_ancho_error=rho("ancho","error_abs"),rho_ancho_error_IQR=rho("ancho_IQR","error_IQR"),
                cobertura95=h.verdad.between(h.lo_mM,h.hi_mM).mean()))
    publicar(pd.DataFrame(calibracion),"calibracion_posterior")
    publicar(pd.DataFrame(asociaciones),"ancho_y_error")
    guardar(z,"incertidumbre_externa.csv.gz")


def ejecutar():
    ish = leer_matriz("ishii")
    h = pd.read_csv(DET/"huecos.csv.gz",float_precision="round_trip")
    orig = describir(ish)
    publicar(orig.reset_index(),"estadistica_original")
    atipico = origen_extremos(ish)
    publicas, todas, cambios_corr, resumen_corr = [], [], [], []
    for solo in (False, True):
        alcance = "huecos_comunes" if solo else "disponibles"
        for met,x in matrices(ish,h,solo):
            d = describir(x)
            d["observados"] = orig.n
            d["imputados"] = d.n - orig.n
            for k in ESTADISTICOS:
                d["delta_"+k] = d[k]-orig[k]
                escala = ish.abs().max()
                if k == "varianza":
                    escala = escala**2
                elif k in ("CV_pct","asimetria"):
                    escala = pd.Series(100 if k=="CV_pct" else 1,index=ish.columns)
                d["cambio_pct_"+k] = relativo(d[k],orig[k],escala)
            d["cociente_SD"] = 1+d.cambio_pct_SD/100
            hueco = x.where(ish.isna())
            d["imputadas_en_mediana"] = hueco.sub(orig.mediana).abs().le(1e-12).sum()
            d["imputadas_extremas"] = (hueco.lt(orig.Q1-1.5*orig.IQR)|hueco.gt(orig.Q3+1.5*orig.IQR)).sum()
            d["estrategia"],d["alcance"] = met,alcance
            todas.append(d.reset_index())
            tocadas = d[d.imputados.gt(0)]
            publicas.append(dict(alcance=alcance,estrategia=met,observados=int(orig.n.sum()),
                 imputados=int(d.imputados.sum()),pendientes=int(d.faltantes.sum()),metabolitos_con_imputacion=len(tocadas),
                 mediana_cociente_SD=tocadas.cociente_SD.median(),SD_disminuye=int((tocadas.cociente_SD<1-1e-12).sum()),
                 mediana_delta_media_mM=tocadas.delta_media.median(),mediana_delta_IQR_mM=tocadas.delta_IQR.median(),
                 en_mediana=int(d.imputadas_en_mediana.sum()),extremas=int(d.imputadas_extremas.sum())))
            for tipo in ("pearson","spearman"):
                antes = correlaciones(ish,tipo).add_suffix("_antes")
                c = antes.join(correlaciones(x,tipo).add_suffix("_despues"))
                c["delta_r"] = c.r_despues - c.r_antes
                c["cambio_relativo_pct"] = 100*c.delta_r / c.r_antes.abs().where(c.r_antes.abs()>=.1)
                c["pares_agregados"] = c.pares_despues-c.pares_antes
                c["estrategia"],c["alcance"],c["correlacion"] = met,alcance,tipo
                cambios_corr.append(c.reset_index())
                afectadas = c[c.pares_agregados.gt(0)].dropna(subset=["delta_r"])
                resumen_corr.append(dict(alcance=alcance,estrategia=met,correlacion=tipo,
                     pares_modificados=len(afectadas),mediana_delta_r=afectadas.delta_r.median(),
                     mediana_abs_delta_r=afectadas.delta_r.abs().median(),P90_abs_delta_r=afectadas.delta_r.abs().quantile(.9)))
    pd.concat(todas).to_csv(RES/"estadistica_por_estrategia.csv.gz",index=False,
                            compression={"method":"gzip","mtime":0},float_format="%.12g")
    guardar(pd.concat(cambios_corr),"correlaciones_por_par.csv.gz")
    publicar(pd.DataFrame(publicas),"efecto_distribucion")
    publicar(pd.DataFrame(resumen_corr),"efecto_correlaciones")
    z = pd.read_csv(DET/"validacion.csv.gz",float_precision="round_trip")
    z = z[z.fold.str.startswith("E") & z.tipo.eq("individual") & z.estrategia.isin(PRINCIPALES)]
    p = z.pivot(index=["cultivo","objetivo"],columns="estrategia",values="pred").dropna()
    verdad = z.drop_duplicates(["cultivo","objetivo"]).set_index(["cultivo","objetivo"]).verdad.reindex(p.index)
    observado = verdad.unstack().reindex_like(ish)
    # Un método puede producir una columna constante. Comparamos correlaciones solo si todos
    # tienen coeficiente definido en exactamente los mismos cultivos.
    metodos_corr = ("pearson", "spearman")
    corr_medida = {tipo: correlaciones(observado, tipo) for tipo in metodos_corr}
    corr_predicha = {(m, tipo): correlaciones(p[m].unstack().reindex_like(ish), tipo)
                    for m in PRINCIPALES for tipo in metodos_corr}
    pares_comunes = {tipo: pd.concat([corr_medida[tipo].r, *[corr_predicha[m,tipo].r for m in PRINCIPALES]],
                                     axis=1).notna().all(axis=1) for tipo in metodos_corr}
    rec, ident = [], []
    for met in PRINCIPALES:
        pred = p[met].unstack().reindex_like(ish)
        so, sp = observado.std(), pred.std()
        ratio = sp/so.where(so.gt(1e-12))
        fila = dict(estrategia=met, celdas=len(p),n_objetivos=int(ratio.count()),
                    mediana_cociente_SD=ratio.median(),mediana_abs_delta_SD=(sp-so).abs().median())
        for tipo in ("pearson","spearman"):
            original = corr_medida[tipo]
            estimada = corr_predicha[met,tipo]
            assert original.pares.equals(estimada.pares)
            delta = (estimada.r-original.r).where(pares_comunes[tipo])
            fila["pares_"+tipo] = delta.count()
            fila["mediana_abs_delta_"+tipo] = delta.abs().median()
        rec.append(fila)
        a = p[met].sub(verdad).abs()
        flag = pd.Series([atipico.at[s,y] for s,y in p.index],index=p.index)
        for nombre, ix in [("todas", flag.notna()),("extremas",flag),("no_extremas",~flag)]:
            ident.append(dict(estrategia=met,verdades=nombre,n=int(ix.sum()),MAE_mM=a[ix].mean()))
    publicar(pd.DataFrame(rec),"recuperacion_experimental")
    publicar(pd.DataFrame(ident),"influencia_extremos")
    info = pd.read_csv(RES/"evidencia_por_metabolito.csv").set_index("objetivo").join(orig)
    asoc = []
    for k in ["faltantes_pct","mediana","SD","CV_pct","media"]:
        t = info[[k,"MAE_KNN3_mM"]].dropna()
        asoc.append(dict(variable=k,metabolitos=len(t),rho_Spearman=t.iloc[:,0].corr(t.iloc[:,1],method="spearman")))
    publicar(pd.DataFrame(asoc),"dificultad_por_metabolito")
    intervalo()


if __name__ == "__main__":
    ejecutar()

