# Evidencia para interpretar predictores auxiliares

Una variable adicional puede aportar información sobre un cultivo y estar disponible solo para algunas muestras. Estos dos catálogos permiten identificar qué mide cada variable y en qué cultivos se puede usar. Las comparaciones de error del protocolo se publican en [protocolo/resultados](../protocolo/resultados/).

| Pregunta | Evidencia y alcance |
| --- | --- |
| ¿Qué tasa se conoce en cada cohorte? | [Tasas documentadas](tasas_documentadas.csv) identifica hoja, fila, unidad, cobertura y rango observado. La captación de carbono se calcula con el número de carbonos del sustrato. El protocolo ensaya su aporte predictivo frente a un control con las mismas filas. |
| ¿Qué proteína puede vincularse a un metabolito? | [Pares proteína–metabolito](pares_de_red.csv) documenta seis hipótesis de reacción, las filas de Ishii y MMC4 y los cultivos con medición. Un enlace funcional permite proponer un predictor. La mejora predictiva se mide por separado para cada objetivo. |
| ¿Cómo se evalúan los auxiliares sobre los mismos casos? | [Contrastes por variable](../protocolo/resultados/auxiliares_pareados.csv) y [comparación de MMC2](../protocolo/resultados/MMC2_pareado.csv) dan los tamaños de muestra y errores pareados. Los [criterios de elección](../protocolo/METODO.md) exigen evidencia dentro del entrenamiento y disponibilidad en la aplicación. |

Las posiciones por cultivo y las comprobaciones contra los libros se reconstruyen mediante [fisiologia.py](../protocolo/fisiologia.py) y [proteinas.py](../protocolo/proteinas.py).
