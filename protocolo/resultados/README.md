# Qué archivo abrir

Los CSV en esta carpeta resumen la regla original, los controles Ishii/MMC2, convergencia y trazabilidad del panel. La ejecución completa escribe **las comparaciones entre estrategias** en `salida_protocolo/revision/resultados/` (o en la carpeta elegida mediante `IPRE_DETALLE`). Estos resultados se calculan a partir de las fuentes y no se incorporan a la base experimental.

| Pregunta | Archivo generado en `revision/resultados/` |
| --- | --- |
| ¿Qué método tiene menos error y con qué cobertura? | `comparacion_estrategias.csv`, `diferencias_pareadas.csv`, `error_por_metabolito.csv.gz`. Filtre prueba, máscara y alcance antes de comparar. |
| ¿Cuánto del error de la regla corresponde a la mediana o a cada decisión? | `descomposicion_regla.csv`, `alternativas_mismas_casillas.csv`, `frecuencia_eleccion.csv` y `soporte_regla_original.csv`. |
| ¿Qué complejidad fue elegida sin consultar la prueba externa? | `seleccion_interna.csv`: fila por configuración y partición. `F` es el entrenamiento para las ausencias reales. |
| ¿Cómo eran las concentraciones antes y cómo cambian? | `estadistica_original.csv`, `estadistica_por_estrategia.csv.gz`, `efecto_distribucion.csv`, `efecto_correlaciones.csv`. |
| ¿Se recupera la estructura cuando conocemos la verdad? | `recuperacion_experimental.csv`, `influencia_extremos.csv`, `error_por_concentracion.csv` y `dificultad_por_metabolito.csv`. |
| ¿El intervalo contiene la medición? | `calibracion_posterior.csv` (cuando lo elige la regla original), `calibracion_iterativo_fijo.csv` (iterativo fijo) y `ancho_y_error.csv`. |
| ¿Qué respaldo tiene un metabolito o un hueco concreto? | `evidencia_por_metabolito.csv`, `evidencia_huecos_por_objetivo.csv`, `abstenciones_fijas.csv` y `estimaciones_por_casilla.csv.gz`. Esta última tiene concentración en mM, método, procedencia y casos comparables. |
| ¿Dónde están las visualizaciones y las huellas de esta ejecución? | `figuras/` y `reproducibilidad.json`. |

En la carpeta superior, `revision/matrices/` contiene una matriz por estrategia. `revision/validacion.csv.gz` conserva los intentos, incluidas las abstenciones; `revision/todas_las_alternativas_pareadas.csv.gz` compara las alternativas disponibles por casilla; `revision/correlaciones_por_par.csv.gz` recoge las diferencias de correlación y el número de cultivos por par; y `revision/extremos_trazados.csv.gz` permite localizar los valores extremos en el Excel. Los archivos `.csv.gz` se abren con `pandas.read_csv` sin descomprimirlos.

Los sorteos se conservan por partición, objetivo y opción bajo `salida_protocolo/posteriores/` (elección original) y `posteriores_fijos/` (iterativo fijo). Cada archivo NPZ incluye los identificadores de casilla y una matriz **sorteos × casillas**, en mM.

