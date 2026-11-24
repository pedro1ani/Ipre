"""Figuras visibles en GitHub. Ejecutar: python -m protocolo.graficos_resumen"""
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator

TABLAS = Path(__file__).resolve().parent / "resultados"
FIGURAS = TABLAS / "figuras"
TINTA, VERDE, AZUL, OCRE = "#163b43", "#197671", "#3d638e", "#b15a30"


def numero(n, decimales=0):
    return f"{n:,.{decimales}f}".translate(str.maketrans(",.", ".,"))


def filas(archivo):
    with (TABLAS / archivo).open(encoding="utf-8", newline="") as fuente:
        return list(csv.DictReader(fuente))


def marco(titulo, subtitulo):
    fig, ejes = plt.subplots(1, 2, figsize=(10.8, 4.65))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=.16, right=.97, top=.68, bottom=.27, wspace=.67)
    fig.text(.035, .92, titulo, fontsize=19, weight="bold", color=TINTA)
    fig.text(.035, .82, subtitulo, fontsize=11, color="#36565e")
    for eje in ejes:
        eje.spines[["top", "right", "left"]].set_visible(False)
        eje.spines["bottom"].set_color("#acbdbd")
        eje.tick_params(axis="y", length=0, labelsize=10, pad=7)
        eje.tick_params(axis="x", colors="#536b70", labelsize=8, length=3)
        eje.grid(axis="x", color="#dce7e6", linewidth=.7)
        eje.set_axisbelow(True)
        eje.xaxis.set_major_locator(MaxNLocator(nbins=5))
    return fig, ejes


def barras(eje, etiquetas, valores, color, limite, formatear, rotulo_dentro=False):
    pos = range(len(valores))
    eje.barh(pos, valores, height=.56, color=color)
    eje.set(yticks=list(pos), yticklabels=etiquetas, xlim=(0, limite))
    eje.invert_yaxis()
    for y, v in enumerate(valores):
        x = v - .02 * limite if rotulo_dentro else v + .018 * limite
        eje.text(x, y, formatear(v), va="center", ha="right" if rotulo_dentro else "left",
                 fontsize=10, color="white" if rotulo_dentro else TINTA, weight="bold")


def guardar(fig, nombre, pie):
    fig.text(.035, .062, pie, fontsize=9.5, color="#36565e")
    FIGURAS.mkdir(exist_ok=True)
    fig.savefig(FIGURAS / nombre, dpi=135, facecolor="white")
    plt.close(fig)


def panel():
    d = json.loads((TABLAS / "ejecucion.json").read_text(encoding="utf-8"))
    fig, (a, b) = marco("¿Qué datos de Ishii se analizan?",
                         "De las señales del laboratorio a las concentraciones disponibles en los cultivos.")
    a.set_title("Selección de señales", loc="left", fontsize=12)
    barras(a, ["Registradas", "Con ≥10 medidas", "Incluidas"],
           [d["filas_originales"], d["con_diez"], d["objetivos"]], VERDE,
           1.25 * d["filas_originales"], numero)
    a.set_xlabel("Número de señales")
    b.set_title("Casillas del panel", loc="left", fontsize=12)
    barras(b, ["Con medición", "Sin medición"], [d["medidos"], d["huecos"]],
           [AZUL, OCRE], 1.27 * d["medidos"], numero)
    b.xaxis.set_major_formatter(FuncFormatter(lambda x, _: numero(x)))
    b.set_xlabel("Número de casillas")
    total = d["cultivos"] * d["objetivos"]
    assert total == d["medidos"] + d["huecos"]
    guardar(fig, "panel_ishii.png", f'{d["cultivos"]} cultivos × {d["objetivos"]} metabolitos = '
            f'{numero(total)} casillas. Un blanco indica que no se informó una concentración.')


def errores():
    datos = [f for f in filas("metricas.csv")
             if f["tipo"] == "individual" and f["alcance"] == "mismas_celdas"]
    metodos = {"A_amplio/Mediana": "Mediana", "A_amplio/KNN3": "KNN3",
               "A_amplio/Iterativo": "Iterativo", "Regla": "Regla original"}
    fig, ejes = marco("Error al reconstruir una concentración medida",
                      "Error absoluto medio (MAE) en mM. Una barra más corta indica menor error.")
    limite = max(float(f["MAE_mM"]) for f in datos if f["metodo"] in metodos) * 1.28
    for eje, (prueba, etiqueta) in zip(ejes, [("perturbacion", "Por perturbación"), ("serie", "Por serie")]):
        t = {f["metodo"]: f for f in datos if f["prueba"] == prueba}
        conteos = {int(t[m]["celdas_distintas"]) for m in metodos}
        assert len(conteos) == 1
        eje.set_title(f"{etiqueta} · {numero(conteos.pop())} casillas", loc="left", fontsize=12)
        barras(eje, list(metodos.values()), [float(t[m]["MAE_mM"]) for m in metodos],
               [AZUL, AZUL, AZUL, VERDE], limite, lambda x: numero(x, 3))
        eje.xaxis.set_major_formatter(FuncFormatter(lambda x, _: numero(x, 2)))
        eje.set_xlabel("MAE (mM)")
    guardar(fig, "error_metodos.png", "Cuatro estrategias de referencia. Dentro de cada prueba, "
            "los cuatro métodos reconstruyen exactamente las mismas mediciones de Ishii.")


def intervalos():
    t = {f["prueba"]: f for f in filas("intervalos.csv") if f["tipo"] == "individual"}
    fig, ejes = marco("Disponibilidad y cobertura del intervalo del 95 %",
                      "Regla original · ocultamiento individual. Cada panel responde a una pregunta distinta.")
    for eje, titulo, num, den in [
        (ejes[0], "¿Cuántas recibieron intervalo?", "con_intervalo", "candidatas"),
        (ejes[1], "¿Cuántos incluyeron la medición?", "dentro", "con_intervalo"),
    ]:
        etiquetas, porcentajes = [], []
        for prueba, nombre in [("perturbacion", "Perturbación"), ("serie", "Serie")]:
            a, b = int(t[prueba][num]), int(t[prueba][den])
            etiquetas.append(f"{nombre}\n{numero(a)} de {numero(b)}")
            porcentajes.append(100 * a / b)
        eje.set_title(titulo, loc="left", fontsize=12)
        barras(eje, etiquetas, porcentajes, VERDE if num == "dentro" else AZUL,
               110, lambda x: numero(x, 2) + " %", rotulo_dentro=num == "dentro")
        eje.set_xticks([0, 25, 50, 75, 100])
        eje.set_xlabel("Porcentaje")
    ejes[1].axvline(95, color=OCRE, linewidth=1.3, linestyle="--", zorder=0)
    guardar(fig, "cobertura_intervalos.png", "Izquierda: de todas las casillas ocultas. Derecha: "
            "de los intervalos emitidos, cuántos contienen la medición. La línea indica el 95 % nominal.")


if __name__ == "__main__":
    panel()
    errores()
    intervalos()

