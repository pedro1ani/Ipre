"""Compara las salidas regeneradas con los CSV previos a la refactorización."""

import hashlib
import json
import unittest
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]


class EquivalenciaTests(unittest.TestCase):
    def test_resultados_conservados(self):
        informe = json.loads(
            (RAIZ / "docs/VERIFICACION_REFACTOR.json").read_text()
        )
        for nombre, esperado in informe["expected_csv_sha256"].items():
            with self.subTest(archivo=nombre):
                actual = hashlib.sha256(
                    (RAIZ / nombre).read_bytes()
                ).hexdigest()
                self.assertEqual(actual, esperado)


if __name__ == "__main__":
    unittest.main()
