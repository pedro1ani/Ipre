"""Verifica esta distribución en solo lectura; no normaliza ni imputa."""
from pathlib import Path
import csv
import hashlib
import json
import math
import sqlite3

ROOT = Path(__file__).resolve().parent

def require(condition, message):
    if not condition:
        raise ValueError(message)

def connect(path):
    return sqlite3.connect(path.resolve().as_uri() + '?mode=ro', uri=True)

def main():
    inventory = json.loads((ROOT / 'MANIFEST_PROYECTO.json').read_text(encoding='utf-8'))
    for name, expected in inventory['files_sha256'].items():
        path = ROOT / name
        require(path.is_file(), 'Falta: ' + name)
        require(hashlib.sha256(path.read_bytes()).hexdigest() == expected, 'Archivo modificado: ' + name)
    original = json.loads((ROOT / 'manifest.json').read_text(encoding='utf-8'))
    for key, folder in [('source_sha256', 'sources'), ('output_sha256', 'data')]:
        for name, expected in original[key].items():
            require(hashlib.sha256((ROOT / folder / name).read_bytes()).hexdigest() == expected, 'Diferencia respecto a conversión original: ' + name)
    with connect(ROOT / 'data/concentraciones_mM.db') as db:
        require(db.execute('PRAGMA integrity_check').fetchone()[0] == 'ok', 'Integridad SQLite')
        require(not db.execute('PRAGMA foreign_key_check').fetchall(), 'Claves foráneas')
        rows = db.execute('SELECT dataset_id, value_original, value_mM, uncertainty_original, uncertainty_mM, conversion_divisor, is_imputed FROM metabolite_measurements_mM').fetchall()
        require(len(rows) == 20483, 'Conteo inesperado')
        for dataset, value, converted, sd, sd_mM, divisor, imputed in rows:
            require(divisor == (2.3 if dataset in ('MMC2', 'MMC3') else 1.0), 'Divisor incorrecto')
            require(imputed == 0, 'Imputación inesperada')
            for a, b in [(value, converted), (sd, sd_mM)]:
                require((a is None) == (b is None), 'Patrón de faltantes modificado')
                if a is not None:
                    require(math.isclose(a / divisor, b, rel_tol=1e-12, abs_tol=1e-14), 'Conversión incorrecta')
        counts = db.execute('SELECT dataset_id, compartment, count(*), count(value_mM) FROM metabolite_measurements_mM GROUP BY 1,2').fetchall()
        sample_keys = {x[0] for x in db.execute("SELECT DISTINCT sample_key FROM metabolite_measurements_mM WHERE dataset_id='ISHII2007'")}
    with connect(ROOT / 'sources/ecoli_auditada_v3.db') as source:
        require(source.execute('PRAGMA integrity_check').fetchone()[0] == 'ok', 'Integridad fuente')
        query = source.execute('SELECT * FROM ishii_samples ORDER BY sample_id')
        expected_rows = [[x[0] for x in query.description]] + [[str(v) if v is not None else '' for v in row] for row in query]
    with (ROOT / 'data/samples_ishii.csv').open(encoding='utf-8-sig', newline='') as file:
        actual = list(csv.reader(file))
    require(actual == expected_rows, 'Metadatos diferentes de la fuente')
    ids = {row[actual[0].index('sample_id')] for row in actual[1:]}
    missing_metadata = sample_keys - ids
    require(missing_metadata <= {'KO05x'}, 'Muestras inesperadas sin metadatos: ' + str(missing_metadata))
    if missing_metadata:
        print('PENDIENTE DOCUMENTADO: KO05x tiene mediciones, pero no fila en ishii_samples; no se fusiona con KO05.')
    print('OK: hashes, SQLite, valores, incertidumbres, faltantes y metadatos.')
    print('Conjunto / compartimento / posiciones / valores observados')
    for row in counts:
        print(' / '.join(map(str, row)))
    print('Estado: conversión verificada. Consultar normalizacion/manifest_normalizacion.json para las matrices; evaluación e imputación separadas en resultados/<dataset>/; consultar docs/VALIDACION_Y_RESULTADOS.md.')

if __name__ == '__main__':
    main()
