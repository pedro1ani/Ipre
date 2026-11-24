"""Ejecutar desde la raíz: python -m protocolo.revision --procesos 4."""
import argparse
import json
from pathlib import Path

from . import comparacion, distribuciones, estrategias, figuras
from .aplicar import aplicar, comprobar
from .datos import PRIV, PUB, guardar_json, huella, preparar
from .informe import informe
from .intervalos import intervalos, intervalos_iterativo_fijo
from .validar import validar, versiones


def registrar():
    base = json.loads((PRIV / "congelado.json").read_text())
    tablas = {p.name: huella(p) for p in sorted(estrategias.RES.glob("*.csv*"))}
    codigo = {p.name: huella(p) for p in sorted(Path(__file__).parent.glob("*.py"))}
    guardar_json(estrategias.RES / "reproducibilidad.json", dict(
        versiones=versiones(), semilla=20261109, bootstrap_por_grupo=2000,
        fuentes=base["entradas"], mascaras=base["mascaras"], codigo=codigo, tablas=tablas,
        tablas_base={p.name: huella(p) for p in sorted(PUB.glob("*.csv"))},
        figuras={p.name: huella(p) for p in sorted(figuras.FIG.glob("*.png"))}))


def resultados_vigentes():
    """Evita mezclar resultados parciales o de otras fuentes al abrir un cuaderno."""
    marca = estrategias.RES / "reproducibilidad.json"
    if not marca.is_file():
        return False
    from .validar import huellas_entradas
    try:
        datos = json.loads(marca.read_text(encoding="utf-8"))
        if datos["versiones"] != versiones() or datos["fuentes"] != huellas_entradas():
            return False
        for carpeta, clave in [(PUB.parent, "codigo"), (estrategias.RES, "tablas"),
                               (PUB, "tablas_base"), (figuras.FIG, "figuras")]:
            if any(huella(carpeta / nombre) != h for nombre, h in datos[clave].items()):
                return False
        return True
    except (OSError, KeyError, ValueError):
        return False


def ejecutar(procesos=4):
    preparar()
    validar(procesos)
    aplicar()
    comprobar()
    estrategias.ejecutar(procesos)
    comparacion.ejecutar()
    intervalos(procesos)
    intervalos_iterativo_fijo(procesos)
    informe()
    distribuciones.ejecutar()
    figuras.ejecutar()
    registrar()


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--procesos", type=int, default=4)
    ejecutar(cli.parse_args().procesos)

