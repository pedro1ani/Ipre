"""Localización de las fuentes, semilla del protocolo e identificadores de señales de Ishii."""
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
SEMILLA = 20261109

# Son señales de Ishii. Ninguna de estas claves designa una medición de Gerosa.
EXTRAS = {
    174: ("ish_g6p", "Glucose 6-phosphate"),
    173: ("ish_f6p", "Fructose 6-phosphate"),
    164: ("ish_ru5p", "Ribulose 5-phosphate"),
    53: ("ish_mal", "Malate"),
    124: ("ish_3pg", "3-Phosphoglycerate"),
    172: ("ish_g1p", "Glucose 1-phosphate"),
    165: ("ish_r5p", "Ribose 5-phosphate"),
}
