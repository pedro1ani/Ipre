# IPre — Validación e imputación de Ishii

Se mantienen las 30 muestras Ishii, siete metabolitos y cuatro objetivos: fumarato, PEP, s7p y succinato. Se comparan dos escenarios separados: entrenamiento Ishii y entrenamiento Ishii + ocho condiciones admisibles de MMC2. Las verdades de evaluación y las estadísticas finales pertenecen exclusivamente a Ishii. Incorporar MMC2 no presupone equivalencia entre estudios.

## Flujo principal reproducible

```bash
python -m pip install -r requirements_imputacion.txt
python ejecutar_estudio.py --output resultados_nuevos
python verificar_estudio.py --results resultados_nuevos
python -m pip install -r requirements_descriptivos.txt
python graficar_estudio.py --results resultados_nuevos
```

La carpeta de salida debe ser nueva. Se reutiliza `validacion_comun/` solo si la huella de las entradas y la colección reconstruida coinciden exactamente. Para una matriz modificada debe seleccionarse una carpeta nueva mediante `--masks mascaras_nuevas`. Se guardan primero todas las evaluaciones y después se estiman los huecos reales. No se ajustan hiperparámetros con los resultados de evaluación.

Los resultados entregados están en `resultados_estudio/`. Leer **docs/VALIDACION_Y_RESULTADOS.md** antes de interpretar las matrices. Contiene diseño, cifras, límites y guía de archivos. **docs/IMPUTACION.md** documenta parámetros, sustituciones y reglas de soporte. No hay KNN propio ni PCA.

## Organización

- `validacion_comun/`: colección única de máscaras, verdades originales Ishii y series analíticas.
- `resultados_estudio/validacion/`: predicciones, MAE en mM y R, cobertura, motivos y comparaciones sobre las mismas celdas.
- `resultados_estudio/imputacion/`: candidatos y matrices por método, escenario y política de aceptación; observaciones originales intactas.
- `resultados_estudio/estadisticas/`: antes/después sobre las 30 muestras, con subconjuntos propios y comunes.
- `resultados_estudio/figuras/`: tres gráficos y PDF.
- `normalizacion/`: matrices fuente separadas Ishii (30 × 7, 138 medidos, 72 ausencias) y MMC2 (8 × 7, 54 medidos, 2 ausencias), referencias positivas de glucosa Gerosa con SD y trazabilidad.
- `data/` y `sources/`: concentraciones, conversiones, metadatos y fuentes conservadas; no se escriben imputaciones en ellas.
- `descriptivos/`: análisis exploratorio de observados por estudio, previo a imputación.

La referencia Gerosa es fija por metabolito para todas las muestras. KNN e Iterative estandarizan exclusivamente con entrenamiento, sin logaritmos. KO05x se conserva con metadatos incompletos; se excluye junto con KO05 al evaluar su perturbación.

## Aplicación independiente posterior a MMC2

El comando auxiliar `python ejecutar_imputacion.py --dataset mmc2 --output resultados_mmc2_nuevos` aplica los mismos métodos exclusivamente a MMC2. No se ejecutó en esta entrega y no sustituye una validación específica de MMC2. Su configuración independiente es `config_imputacion_mmc2.json`; piruvato sigue siendo predictor, no objetivo. El comando auxiliar también admite `--dataset ishii`, pero para evaluar antes de imputar debe utilizarse el flujo principal anterior.

## Comprobaciones y reconstrucción

`python verificar_proyecto.py` verifica hashes, integridad y conversiones. `python verificar_estudio.py` comprueba particiones, escalado, verdades, máscaras, errores, observaciones intactas y comparaciones comunes. `python -m unittest discover -s tests -v` ejecuta 15 pruebas.

`python normalizar.py --output normalizacion_nueva` reconstruye matrices separadas. `python analizar_descriptivos.py --output descriptivos_nuevos` reconstruye descriptivos observados. `convertir_unidades.py` reconstruye la capa mM en una copia separada con data/ vacía, usando requirements.txt; no regenera todos los manifiestos. No hace falta ejecutarlo sobre esta distribución ni usar Python -O.

Semilla, orden, versiones y hashes quedan registrados en configuraciones y manifiestos. Dos ejecuciones completas produjeron 24 archivos CSV/JSON idénticos byte a byte en el mismo entorno. Los gráficos no están incluidos en esa comparación binaria.
