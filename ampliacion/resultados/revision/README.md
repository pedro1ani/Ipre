# Tablas para leer la revisión

| Archivo | Pregunta que responde |
| --- | --- |
| [Resumen de comparaciones](resumen_comparaciones.csv) | ¿Cómo cambia el error global con la escala, la captación, el flujo o la dilución? |
| [Error del panel ampliado](error_panel_ampliado.csv) | ¿Qué MAE, RMSE y sesgo tiene cada método para cada metabolito, en las mismas celdas? |
| [Aporte de la captación](aporte_captacion.csv) | ¿Qué metabolitos mejoran en sus mismas celdas si se añade la captación como predictor? Un delta negativo indica menor error. |
| [Equivalencia de escalas](equivalencia.csv) | ¿Cuánto cambian las predicciones en mM al dividir por Gerosa antes de estandarizar? |

En el resumen, `comunes_familia` usa las mismas celdas para todos los métodos, estrategias y escenarios de una misma prueba; `comunes_todas` usa la intersección de *todos* los diseños. El denominador original son 457 mediciones ocultas. Las tablas enlazadas se muestran con seis cifras significativas. No contienen identificadores de cultivo.

Para reproducir el análisis desde la raíz del repositorio, ejecute `python -m ampliacion.revisar_escalas` o la sección 1.1 del [notebook de evaluación](../../../02_imputacion_evaluacion.ipynb). Se necesitan los resultados versionados de la validación base; si han cambiado sus entradas, ejecute antes `python -m ampliacion.validar`.

En este mismo directorio se generan localmente los datos completos: `condiciones.csv` y `fuentes_fisiologicas.csv.gz` identifican la celda Excel y distinguen ausencias; `inventario.csv` cuenta los datos disponibles por fila; `ajustes.csv` registra entrenamiento y convergencia; `predicciones.csv.gz` contiene los intentos de esta revisión; `global.csv`, `por_metabolito.csv`, `por_grupo.csv.gz` y `pesos_Gerosa.csv` permiten repetir las comparaciones sin redondeo de presentación. Para el panel ampliado de referencia, las predicciones provienen de la validación base, en el directorio superior. Las tablas por cultivo se generan al ejecutar, sin necesidad de publicarlas.
