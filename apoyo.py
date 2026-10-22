"""Operaciones compartidas de archivo y estadística descriptiva."""

import hashlib
import json
import math


def huella(archivo):
    """Detecta cambios de contenido entre ejecuciones."""
    return hashlib.sha256(archivo.read_bytes()).hexdigest()


def guardar_json(archivo, contenido):
    texto = json.dumps(
        contenido, indent=2, ensure_ascii=False, allow_nan=False
    )
    archivo.write_text(texto + "\n", encoding="utf-8")


def guardar_csv(archivo, tabla):
    tabla.to_csv(
        archivo, index=False, encoding="utf-8-sig", float_format="%.17g"
    )


def describir(valores):
    """Omite ausencias y señala cuándo el CV carece de interpretación."""
    valores = valores.dropna()
    media = valores.mean()
    desviacion = valores.std(ddof=1)
    if len(valores) < 2:
        estado_cv = "n_less_than_2"
    elif media <= 0:
        estado_cv = "mean_nonpositive"
    elif abs(media) <= 1e-12:
        estado_cv = "mean_near_zero"
    else:
        estado_cv = "interpretable"
    cv = math.nan
    if estado_cv == "interpretable":
        cv = 100 * desviacion / media
    return {
        "mean_mM": media,
        "median_mM": valores.median(),
        "SD_mM": desviacion,
        "IQR_mM": valores.quantile(0.75) - valores.quantile(0.25),
        "CV_pct": cv,
        "CV_status": estado_cv,
    }
