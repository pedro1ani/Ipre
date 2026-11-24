# IPre · Estimación de concentraciones de metabolitos en *Escherichia coli*

Un **metabolito** es una sustancia que participa en las reacciones de una célula. El estudio de Ishii mide la concentración de estas sustancias dentro de la bacteria *Escherichia coli*, pero para algunos cultivos no informa un valor. Este proyecto estudia **cuándo se puede estimar esa concentración a partir de las mediciones disponibles y con qué error**. Ese cálculo se llama *imputación*. Una estimación se registra por separado de lo que mide el laboratorio.

| Cuaderno | Qué explica |
| --- | --- |
| [01 · Datos y selección de metabolitos](01_auditoria_preparacion.ipynb) | Qué mide cada experimento, por qué algunas señales entran en el análisis y cómo encontrar su celda de origen. |
| [02 · Evaluación y estimación de concentraciones](02_imputacion_evaluacion.ipynb) | Cómo se prueban los métodos con mediciones conocidas y qué decisión se toma cuando falta un valor. |

Los dos cuadernos se pueden leer en GitHub sin instalar programas.

## Qué contiene la tabla de datos

Cada fila representa uno de los **30 cultivos de Ishii**. Cada columna representa una señal del laboratorio y cada casilla contiene una concentración en **mM** (milimoles por litro de volumen intracelular) o queda vacía.

El [inventario](protocolo/resultados/inventario.csv) contiene 579 señales, que no necesariamente son 579 sustancias distintas. De ellas, 149 tienen mediciones en al menos diez cultivos. En 14 de esas señales la identidad es incierta o se registra una mezcla. **Las 135 restantes forman la tabla de análisis.** Para 17 de ellas también hay mediciones comparables en un segundo estudio, el de Gerosa. Su archivo **MMC2** aporta ocho condiciones de cultivo que se prueban como información adicional. Las otras 118 señales se estudian solo con Ishii.

## Dónde están las fuentes y los resultados

| Ubicación | Qué contiene |
| --- | --- |
| [data](data/) | La **base de trabajo permanente**, `concentraciones_mM.db`, sus exportaciones por estudio y la evidencia de conversión. Reúne 20.483 posiciones de concentración, incluidas las que carecen de medición. |
| [sources](sources/) | Los Excel y documentos experimentales de Ishii y Gerosa y la base multi-ómica `ecoli_auditada_v3.db`. Esta última también contiene datos de proteínas, fisiología y flujos. |
| [auditoria](auditoria/) | Mediciones con sus coordenadas originales, homologación química y metadatos de cultivos. Son las entradas tabulares que la preparación coteja con los Excel. |
| [protocolo](protocolo/METODO.md) | Preparación, diseño de validación, modelos, elección por casilla, intervalos y síntesis. Las **13 tablas y resúmenes científicos vigentes** se publican en [resultados](protocolo/resultados/). |
| [estudios_complementarios](estudios_complementarios/) | Documentación de tasas y pares proteína–metabolito que permite interpretar los predictores auxiliares. |
| [tests](tests/) | Comprobaciones de las fuentes, el aislamiento de cultivos, las máscaras, los controles y los productos por cultivo. |

MMC2 aporta las ocho condiciones estacionarias de Gerosa. MMC3 contiene mediciones a lo largo de transiciones y se conserva para consulta. MMC4 documenta reacciones y ayuda a interpretar los vínculos con proteínas. El protocolo prueba el aporte de MMC2 utilizando mediciones de Ishii como referencia de evaluación.

## Cómo se evalúa la estimación

Se oculta una concentración conocida de Ishii, se la estima y se compara el resultado con la medición original. Los cultivos con la misma perturbación experimental, por ejemplo la inactivación del mismo gen, se mantienen juntos. Cuando se evalúa uno de esos grupos, el método aprende solo de los demás.

| Estrategia | Cómo obtiene el valor | Error medio por medición (mM) |
| --- | --- | ---: |
| Mediana | Usa el valor central del metabolito en los cultivos de entrenamiento. | 0,246 |
| Vecinos (KNN), con tres vecinos | Promedia las concentraciones del metabolito en cultivos parecidos según las otras sustancias medidas. | 0,163 |
| Método iterativo | Estima la concentración a partir de relaciones entre metabolitos y refina el cálculo en varias rondas. | 0,168 |
| Regla de elección | Elige un método cuando las pruebas lo respaldan para ese metabolito y ese patrón de datos faltantes. | 0,204 |

Los cuatro resultados corresponden a **las mismas 3.302 mediciones ocultas** de 25 grupos. El error es la distancia promedio a la medición, tanto si el método se queda corto como si se excede. KNN con tres vecinos tiene el menor error promedio en esta prueba; el [resultado por metabolito](protocolo/resultados/por_objetivo.csv) permite revisar una sustancia en particular.

## Qué se obtiene para los valores faltantes

La tabla contiene **3.402 mediciones originales y 648 casillas vacías**. Ishii usa el mismo blanco cuando una sustancia no se detecta o queda por debajo del límite de detección. Por eso un blanco no equivale a cero y su estimación se interpreta como exploratoria.

| De las 648 casillas vacías | Casillas | Interpretación |
| --- | ---: | --- |
| Método elegido en las pruebas | 221 | Hay evidencia para elegir el método en ese caso; la concentración real sigue siendo desconocida. |
| Mediana de respaldo | 418 | Se informa un valor central del entrenamiento sin evidencia suficiente para preferir un método para ese caso. |
| Sin estimación | 9 | Los datos disponibles no permiten aplicar la regla. La casilla queda vacía. |

**Ejemplo: fumarato.** Ishii informa 23 concentraciones en 30 cultivos. Para dos de las siete casillas vacías se elige KNN, para tres se elige el método iterativo y para dos se usa la mediana de respaldo. Las 23 mediciones conservan su valor original.

El [protocolo de análisis](protocolo/METODO.md) explica los criterios y la organización del código. También se pueden consultar las tablas de [errores](protocolo/resultados/metricas.csv), [aporte de Gerosa](protocolo/resultados/MMC2_pareado.csv), [intervalos de predicción](protocolo/resultados/intervalos.csv) y [estadísticas por metabolito](protocolo/resultados/descriptivos.csv).

<details>
<summary><strong>Cómo ejecutar el análisis y dónde encontrar los archivos por cultivo</strong></summary>

Desde la carpeta principal del repositorio, con Python 3.12.14 en un entorno virtual, estos comandos preparan los datos y evalúan los métodos sobre mediciones conocidas:

```bash
python -m pip install -r requirements.txt
python -m auditoria.verificar
python -m protocolo.datos
python -m protocolo.validar --procesos 4
python -m protocolo.informe --solo-validacion
```

Los siguientes comandos aplican la regla a los valores faltantes, calculan los intervalos y generan las tablas del informe:

```bash
python -m protocolo.aplicar
python -m protocolo.intervalos --procesos 4
python -m protocolo.informe
python -m unittest discover -s tests
```

Los resúmenes se guardan en `protocolo/resultados/`. Los datos por cultivo se generan en `salida_protocolo/`:

| Archivo generado | Contenido |
| --- | --- |
| `matriz_solo_con_eleccion_validada.csv` | Las 3.402 mediciones y las 221 estimaciones con elección de método. Las otras 427 casillas siguen vacías. |
| `matriz_regla.csv` | Las mismas mediciones y estimaciones, más las 418 medianas de respaldo. Quedan nueve casillas vacías. |
| `auditoria_huecos.csv` | Una fila por casilla vacía: cultivo, metabolito, celda original, decisión, método y evidencia disponible. |
| `original_30_completa.csv` y `resultado_30_completo.csv` | Las tablas con las 579 señales del inventario, antes y después de aplicar la regla a las 135 seleccionadas. |

Para guardar el detalle en otra carpeta, establece la variable de entorno `IPRE_DETALLE`. Para trabajar en Jupyter, selecciona el mismo entorno de Python como kernel. Ejecuta el cuaderno 01 desde la raíz y luego el 02. Este último abre los resúmenes publicados por defecto y, si ya se generó el detalle, muestra también ejemplos por cultivo.

</details>
