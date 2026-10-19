# Diseño del flujo simplificado

Un único ejecutor, ejecutar_imputacion.py, analiza una base seleccionada y procesa dos secciones consecutivas: evaluación con máscaras y aplicación a candidatos reales. No mezcla estudios ni ejecuta una comparación Ishii + MMC2. Los resultados combinados del flujo anterior se retiraron de la distribución activa; la versión anterior conserva ese análisis.

## Evaluación

La colección determinista de máscaras se genera una vez por ejecución y se comparte entre los tres métodos. Guarda orden de filas/columnas, fuente de patrones, exclusiones y huellas de entradas en validacion/mascaras.json. La semilla es 20261019. Se guardan las verdades originales y nunca se usan los valores ocultados para entrenar o transformar la consulta. Una nueva ejecución reconstruye las máscaras, sin cargar colecciones de otra matriz o estudio.

| Base | Prueba por grupo/condición | Patrones reales | Serie analítica | Total de celdas por método |
|---|---:|---:|---:|---:|
| Ishii | 101 | 233 | 101 | 435 |
| MMC2 | 27 | 5 | No aplicable | 32 |

Ishii conserva las tres pruebas. La prueba de serie excluye la unión de familia de serie y grupo de perturbación; 1 y 1* forman una familia conservadora. La prueba de patrones reproduce simultáneamente todas las ausencias de las siete variables y evalúa solo objetivos originalmente medidos. Requiere que todas las variables visibles estén originalmente observadas; puede favorecer muestras más completas. Una medición puede repetirse en pruebas distintas: los intentos no son muestras independientes.

MMC2 se evalúa usando exclusivamente mediciones MMC2 y excluyendo la condición consultada del entrenamiento. No se inventa una serie analítica equivalente a Ishii. Se excluye MMC2_GLUCOSE como consulta evaluada porque sus concentraciones están reveladas en el denominador fijo de normalización. La referencia sigue siendo externa/fija para las otras condiciones. Los 32 casos de MMC2 constituyen una evaluación pequeña y exploratoria; no deben equipararse a los casos Ishii.

predicciones.csv contiene cada valor medido, imputado y error en mM/R. metricas.csv contiene MAE por metabolito y prueba, cobertura propia y comparación sobre las mismas celdas. cobertura_y_motivos.csv explica ausencias de predicción o rechazos por soporte; casos_comunes.csv identifica las intersecciones. Un MAE sin predicciones permanece vacío.

## Candidatos reales

Candidatos = celdas originalmente ausentes en los cuatro objetivos. Ishii tiene 19 y MMC2 tiene 1 (succinato en condición succinato); el faltante de piruvato de MMC2 no es objetivo. Se generan 57 y 3 intentos, respectivamente, porque cada candidato se procesa con los tres métodos. Las observaciones se conservan intactas. Todos los entrenamientos usan exclusivamente la base seleccionada, excluyendo el grupo/condición de la consulta.

Las matrices available incluyen todos los candidatos finitos. Las matrices supported_plus_baseline aceptan solo supported de KNN e Iterative y mantienen la mediana como referencia general. Más cobertura no demuestra más exactitud. No hay sustitución silenciosa entre métodos: la falta de convergencia conserva una estimación etiquetada como limitada.

estadisticas_antes_despues.csv informa, sobre todas las muestras de la base seleccionada, observados, imputados, pendientes, media, mediana, SD muestral, IQR y CV en mM. Incluye resultados propios y sobre los mismos huecos aceptados por los tres métodos. Los predictores conservan sus faltantes originales. CV no interpretable si n<2, media no positiva o magnitud ≤1e-12 mM. Tratar los imputados como puntos no incorpora incertidumbre de imputación; reducir el CV no demuestra mejor calidad.

## Simplificación y comprobación

Se retiraron el ejecutor redundante ejecutar_estudio.py, la configuración combinada y las salidas mezcladas de la carpeta activa. Se mantienen las fuentes y preparación por trazabilidad. El motor sklearn no cambió. La comparación con el flujo anterior verificó igualdad de las 1.305 predicciones y etiquetas Ishii. Las dos bases se ejecutaron independientemente y se comprobaron entrenamiento exclusivo, máscaras comunes, unidades, observados intactos y conteos de estadísticas.

Los gráficos real_vs_imputado muestran únicamente ocultaciones de valores conocidos; la diagonal es predicción exacta. Los colores distinguen candidatos con soporte de generales/limitados. No se conoce el error verdadero de los huecos reales.
