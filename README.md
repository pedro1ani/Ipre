# IPre · Imputación de concentraciones intracelulares de *E. coli*

La pregunta es **qué se puede estimar para un cultivo de Ishii, con qué método y con qué evidencia**. Primero ocultamos mediciones originales y retiramos sus grupos del entrenamiento. Solo después aplicamos la regla a los huecos reales. El resultado se expresa en mM y conserva cada medición original.

| Notebook | Para qué abrirlo |
| --- | --- |
| [01 · ¿Qué midió Ishii y qué entra al análisis?](01_auditoria_preparacion.ipynb) | Revisar las 579 filas analíticas, el filtro de identidad y disponibilidad y la celda de origen de cada concentración. |
| [02 · ¿Qué método usar para este hueco?](02_imputacion_evaluacion.ipynb) | Comparar métodos sobre las mismas mediciones ocultas, ver cuándo aporta MMC2 y revisar por separado las estimaciones reales y las de respaldo. |

**Resultado de esta versión.** De 579 señales de Ishii, 149 tienen diez o más mediciones en los 30 cultivos. Se retuvieron 135 señales nominales: 17 compartidas con Gerosa y 118 exclusivas de Ishii. Los 648 huecos del panel se distribuyen en **221 estimaciones con elección validada internamente, 418 medianas de respaldo y 9 pendientes**. La matriz conserva las 3.402 mediciones originales del panel. El fichero con *solo* la elección validada deja expresamente en blanco los respaldos y los pendientes.

| Método | MAE micro en mM | MAE macro en mM |
| --- | ---: | ---: |
| Mediana fija | 0,246 | 0,220 |
| KNN fijo, 3 vecinos | 0,163 | 0,147 |
| Iterativo fijo | 0,168 | 0,151 |
| Regla por objetivo y perfil | 0,204 | 0,184 |

Comparación externa individual por perturbación: **las mismas 3.302 mediciones de 25 grupos**. La regla reduce el error frente a la mediana, pero los otros dos controles fijos tuvieron menor error global. Selección interna no equivale a una garantía para un hueco real. En el Excel de Ishii un blanco puede indicar «no detectable» o «por debajo del límite» sin especificar cuál: toda estimación de esos blancos es exploratoria.

**Para revisar los resultados:** [método y límites](protocolo/METODO.md), [inventario de las 579 señales](protocolo/resultados/inventario.csv), [comparación global](protocolo/resultados/metricas.csv), [error por objetivo en celdas comunes](protocolo/resultados/por_objetivo.csv), [aportación de MMC2](protocolo/resultados/MMC2_pareado.csv), [intervalos](protocolo/resultados/intervalos.csv) y [estadísticas antes/después](protocolo/resultados/descriptivos.csv).

## Reproducir

Desde la raíz, con Python 3.12 y un entorno virtual (para abrir los notebooks, seleccionarlo como kernel de Jupyter):

```bash
python -m pip install -r requirements.txt
python -m protocolo.datos
python -m protocolo.validar --procesos 4
python -m protocolo.informe --solo-validacion
```

Aquí se puede detener la ejecución y revisar el error externo **antes de imputar**. Para completar los huecos autorizados:

```bash
python -m protocolo.aplicar
python -m protocolo.intervalos --procesos 4
python -m protocolo.informe
python -m unittest discover -s tests
```

Los agregados quedan en `protocolo/resultados/`. En `salida_protocolo/` se generan las matrices, `auditoria_huecos.csv`, las máscaras por partición, los ajustes y la regla congelada. Ese directorio no se sube a Git. `original_30_completa.csv` y `resultado_30_completo.csv` conservan las 579 columnas; `matriz_solo_con_eleccion_validada.csv` contiene el panel de 135. Cada ejecución puede dirigirse a un directorio distinto con la variable de entorno `IPRE_DETALLE`.

Los [bloques 1–4](ampliacion/BLOQUES_1_4.md), [5–6](ampliacion/BLOQUES_5_6.md) y [7–8](ampliacion/BLOQUES_7_8.md) documentan los ensayos anteriores. Los [notebooks anteriores a este cambio](https://github.com/pedro1ani/Ipre/tree/2fac38dca143f872a6b9d8c8db7377e9d6690d13) y la [auditoría original](historico/) siguen disponibles. Sus cifras no son resultados del panel de 135 señales.
