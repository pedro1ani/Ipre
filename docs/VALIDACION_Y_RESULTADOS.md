# Validación común y resultados

## Diseño fijado

Cohorte de evaluación y estadísticas: 30 muestras originales Ishii, siete metabolitos, cuatro objetivos (fum, pep, s7p, succ). Dos entrenamientos separados: Ishii e Ishii + las ocho condiciones MMC2 admisibles. Todas las verdades proceden exclusivamente de mediciones originales Ishii; ninguna estimación ni medición MMC2 sirve de verdad evaluada.

La referencia positiva Gerosa glucosa es idéntica entre muestras y escenarios, en mM, con incertidumbre documentada. Se transforma R a mM multiplicando por esa referencia. No se propaga la incertidumbre del denominador: falta información de covarianzas. La fila Gerosa glucosa puede entrar como entrenamiento auxiliar, pero no se evalúa como verdad ocultada.

Hay una única colección `validacion_comun/masks.json`, compartida por los seis resultados. La huella incluye matrices, referencia, configuraciones y base fuente; además se exige igualdad exacta de la colección reconstruida para reutilizarla. Las máscaras anteriores de observado no constituían máscaras artificiales de evaluación y se creó esta colección nueva.

| Prueba | Ensayos | Celdas evaluadas por método y escenario | Exclusión del entrenamiento |
|---|---:|---:|---|
| Perturbación | 101 | 101 | Grupo completo de la consulta; se oculta un objetivo medido |
| Patrones reales | 97 | 233 | Grupo completo de la consulta; se oculta simultáneamente el patrón |
| Serie analítica | 101 | 101 | Unión de la familia de serie y el grupo de perturbación; se oculta un objetivo medido |

Total: 299 ensayos y 435 celdas por método/escenario, 2.610 intentos. Una medición puede reaparecer en distintos ensayos: no son 435 muestras independientes. No se calculan intervalos ni pruebas de significación asumiendo independencia.

Los seis patrones únicos se extraen de las siete variables en filas con algún objetivo ausente. Para aplicar un patrón a otra fila, todos los objetivos puntuados y todas las variables que quedarán visibles deben estar originalmente medidos. Así, la ausencia resultante coincide exactamente con el patrón. Los predictores ocultados no se puntúan. El patrón totalmente vacío se conserva y permite medir el comportamiento de sustitución general. La elegibilidad favorece filas más completas; los errores no representan sin sesgo todos los huecos reales.

La serie proviene de las columnas analíticas originales de Ishii y queda exportada con trazabilidad. Se agrupan conservadoramente 1 y 1* en la familia 1. KO05x se conserva como fila independiente; KO05 y KO05x se excluyen juntos por perturbación. Todos los controles WT se excluyen juntos. Esta prueba por serie es más estricta al excluir también el grupo de perturbación.

Máscaras exhaustivas y deterministas, semilla 20261019, orden fijo y versiones fijadas en requirements_imputacion.txt. No hay logaritmos en la comparación. StandardScaler y los modelos se ajustan solo con entrenamiento; la consulta nunca participa. No hay búsqueda de hiperparámetros posterior a observar el error. Cada predicción referencia un fit_id con las filas de entrenamiento y parámetros aprendidos. Se finaliza la evaluación antes de producir imputaciones reales.

## Cobertura y comparabilidad

Se entregan dos políticas:

- `available`: todos los números finitos; incluye sustituciones basadas en entrenamiento y candidatos sin convergencia.
- `supported_plus_baseline`: KNN e Iterative requieren etiqueta supported; SimpleImputer se conserva como referencia de mediana general, explícitamente exceptuada del requisito de soporte individual.

Los motivos de rechazo operativo o abstención figuran en cobertura_y_motivos.csv; los diagnósticos completos, en predicciones.csv. La cobertura numérica fue 100% en los seis resultados, pero eso no implica precisión ni soporte suficiente.

Los ámbitos `own`, `common_methods` y `common_all_six` comparan, respectivamente, celdas disponibles propias, intersección entre tres métodos por escenario e intersección entre seis resultados. Un caso incluye el ensayo y la celda: no se empareja solamente por muestra/metabolito. efecto_MMC2_pareado.csv compara ambos escenarios sobre los mismos casos para cada método; delta MAE positivo significa empeoramiento al añadir MMC2. Los MAE sin celdas comunes quedan vacíos, nunca cero.

En la política operativa, la intersección de los seis resultados es muy pequeña: por perturbación solo 1 fumarato y 1 PEP; por serie, 5 y 5; por patrones, 1 y 1. No hay s7p ni succinato comunes a los seis bajo esta política. No se debe construir una clasificación general de métodos a partir de esas intersecciones pequeñas.

## MAE por perturbación: mismos casos, todos los candidatos disponibles

Esta tabla incluye sustituciones generales y candidatos sin convergencia; debe leerse junto con cobertura y diagnósticos.

| Escenario | Método | Metabolito | n | MAE mM | MAE R |
|---|---|---|---:|---:|---:|
| ishii | IterativeImputer | fum[c] | 23 | 0.022466 | 0.040978 |
| ishii | IterativeImputer | pep[c] | 22 | 0.056839 | 0.160946 |
| ishii | IterativeImputer | s7p[c] | 29 | 0.577989 | 0.160717 |
| ishii | IterativeImputer | succ[c] | 27 | 0.032808 | 0.027678 |
| ishii | KNNImputer | fum[c] | 23 | 0.022097 | 0.040305 |
| ishii | KNNImputer | pep[c] | 22 | 0.064333 | 0.182166 |
| ishii | KNNImputer | s7p[c] | 29 | 0.234675 | 0.065254 |
| ishii | KNNImputer | succ[c] | 27 | 0.024347 | 0.020540 |
| ishii | SimpleImputer | fum[c] | 23 | 0.017250 | 0.031464 |
| ishii | SimpleImputer | pep[c] | 22 | 0.048723 | 0.137963 |
| ishii | SimpleImputer | s7p[c] | 29 | 0.220048 | 0.061187 |
| ishii | SimpleImputer | succ[c] | 27 | 0.022065 | 0.018615 |
| ishii_mmc2 | IterativeImputer | fum[c] | 23 | 0.024886 | 0.045391 |
| ishii_mmc2 | IterativeImputer | pep[c] | 22 | 0.046141 | 0.130654 |
| ishii_mmc2 | IterativeImputer | s7p[c] | 29 | 0.399906 | 0.111199 |
| ishii_mmc2 | IterativeImputer | succ[c] | 27 | 0.049609 | 0.041853 |
| ishii_mmc2 | KNNImputer | fum[c] | 23 | 0.017390 | 0.031719 |
| ishii_mmc2 | KNNImputer | pep[c] | 22 | 0.058771 | 0.166416 |
| ishii_mmc2 | KNNImputer | s7p[c] | 29 | 0.251977 | 0.070065 |
| ishii_mmc2 | KNNImputer | succ[c] | 27 | 0.039045 | 0.032940 |
| ishii_mmc2 | SimpleImputer | fum[c] | 23 | 0.020270 | 0.036972 |
| ishii_mmc2 | SimpleImputer | pep[c] | 22 | 0.050395 | 0.142700 |
| ishii_mmc2 | SimpleImputer | s7p[c] | 29 | 0.229652 | 0.063857 |
| ishii_mmc2 | SimpleImputer | succ[c] | 27 | 0.024533 | 0.020698 |

En Ishii, la mediana tiene el menor MAE de los tres métodos para los cuatro objetivos en esta prueba. Esto no demuestra que sea el mejor método para cada hueco real. Añadir MMC2 empeora la mediana en los cuatro objetivos y las tres pruebas. KNN mejora fumarato y PEP en perturbación y serie, pero empeora otros objetivos; bajo patrones reales empeora en los cuatro. Iterative presenta cambios mixtos y convergencia problemática.

Entre los 37 entrenamientos distintos por escenario, Iterative converge en 31 con Ishii y solo 4 al incorporar MMC2. La mejora numérica de algunos MAE combinados no basta para recomendar ese resultado. Normalizar no elimina las diferencias de distribución entre estudios.

## Imputación posterior y estadísticas

Se estiman por separado los 19 huecos objetivo de Ishii para cada método y escenario (114 candidatos). Se excluye el grupo completo de perturbación antes de cada consulta real. Las mediciones originales se conservan exactamente y los tres predictores no se completan en las matrices finales. No se han imputado las muestras MMC2.

| Escenario | Mediana: referencia general | KNN con soporte | Iterative con soporte | Celdas comunes a los tres, política operativa |
|---|---:|---:|---:|---:|
| Ishii | 19 | 5 | 12 | 5 |
| Ishii + MMC2 | 19 | 9 | 0 | 0 |

La política available entrega 19 candidatos en cada resultado, por lo que su intersección común contiene los 19 huecos. La intersección operativa entre los seis resultados está vacía. Las estadísticas comunes operativas del escenario combinado y de los seis resultados coinciden con los observados: eso indica que no hay imputaciones comparables aceptadas, no ausencia de efecto del método.

`estadisticas/antes_despues.csv` contiene 504 filas, en mM, sobre las mismas 30 muestras: antes/después, dos escenarios, tres métodos, dos políticas, tres ámbitos y siete metabolitos. Incluye observados, imputados, pendientes, n disponible, media, mediana, SD muestral (ddof=1), Q1, Q3, IQR y CV porcentual. Para cada metabolito, observados + imputados + pendientes = 30. Los tres predictores conservan las estadísticas originales. `casillas_comunes.csv` identifica la pertenencia de cada hueco a las comparaciones comunes.

El CV se deja ausente y se documenta como no interpretable si n < 2, la media es no positiva o su magnitud es ≤ 1e-12 mM. En estos resultados no se activó ninguno de esos casos. Los CV altos, como s7p, describen dispersión; una caída del CV tras sustituir por medianas no demuestra mejor precisión. Las estimaciones se tratan como puntos al calcular estadísticas, sin incorporar incertidumbre de imputación. El patrón de faltantes puede sesgar los descriptivos antes y después.

## Archivos y reproducción

- `validacion/predicciones.csv`: cada intento con verdad original, estimación R/mM, errores, soporte y fit_id.
- `validacion/metricas.csv`: MAE R/mM, sesgo y cobertura por objetivo, prueba, política, escenario y ámbito.
- `validacion/cobertura_y_motivos.csv`, `casos_comunes.csv`, `efecto_MMC2_pareado.csv`: denominadores, razones e intersecciones auditables.
- `imputacion/candidatos_y_soporte.csv`: los 114 intentos reales; doce matrices mM por política/escenario/método.
- `ajustes_y_particiones.json`: entrenamientos, escaladores y convergencia; `manifest_estudio.json`: entradas, código, versiones y salidas de cálculo.
- `figuras/`: MAE por perturbación, cobertura operativa y distribución antes/después; PDF conjunto. La figura antes/después muestra todos los candidatos, incluidos los no convergentes o sin soporte.
- `reproducibilidad.json`: dos ejecuciones completas con 24 CSV/JSON idénticos en el mismo entorno. Las figuras y este registro suplementario no forman parte de esa comparación binaria; quedan inventariados en MANIFEST_PROYECTO.json.

Se superaron 15 pruebas y las comprobaciones de integridad, particiones, entrenamiento exclusivo, máscaras, errores, estadísticas y preservación de observaciones. Consulte README.md para comandos. Estos resultados evalúan ocultaciones de observados, no permiten conocer el error verdadero de las 19 ausencias reales.
