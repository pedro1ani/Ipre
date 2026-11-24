# Origen de las concentraciones

Cada posición del registro representa una señal, un cultivo o condición y una celda de la fuente. Incluye las posiciones para las que el laboratorio dejó un blanco. El identificador `measurement_uid` permite seguir la misma posición en las dos bases y en el CSV. La columna `analyte_id` identifica la señal analítica y `sample_key` el cultivo o condición.

| Recurso | Contenido y uso |
| --- | --- |
| [Base de trabajo en mM](../data/concentraciones_mM.db) | La tabla `metabolite_measurements_mM` contiene 20.483 posiciones. Las vistas `intracellular_mM` y `extracellular_mM` reúnen 20.219 y 264. La tabla `conversion_factors` documenta el volumen de conversión. La base conserva el valor y unidad originales, la fuente y el valor en mM. |
| [Base experimental auditada](../sources/ecoli_auditada_v3.db) | La vista `audit_measurements_readable` reúne 22.761 posiciones de metabolitos y proteínas. `audit_analytes` identifica analitos, `audit_sources` identifica fuentes y `audit_reason_codes` explica estados de ausencia. Las tablas de corridas, proteínas, condiciones, fisiología y reacciones aportan datos complementarios. |
| [Mediciones auditadas](mediciones_auditadas.csv) | Entrada del protocolo con las 20.483 posiciones de metabolitos de Ishii, MMC2 y MMC3. Registra valores originales y convertidos, motivo del blanco, incertidumbre publicada, hoja y coordenada. Los campos de panel del inventario de origen se conservan por trazabilidad. La selección vigente se calcula en [datos.py](../protocolo/datos.py). |
| [Homologación maestra](homologacion_maestra.csv) | Identidad, fila, compartimento, protocolo, correspondencia química y justificación. Contiene registros de señales y entradas documentales del modelo de Gerosa. |
| [Metadatos](metadatos.csv) | Estudio, cultivo, perturbación, tasa de dilución, serie de medición y pertenencia a la cohorte. Identifica también controles y registros agregados. |
| [Conversiones](conversiones.csv) | Unidades, factor aplicado a cada fuente, evidencia y limitaciones. Los archivos de [data](../data/) conservan además las exportaciones por estudio. |

La preparación lee los tres CSV de mediciones, identidades y metadatos. Coteja cada posición de concentración utilizada con su celda de Excel. Lee directamente la base experimental para la fisiología, las proteínas y sus corridas. La base de trabajo permite consultar el conjunto completo en mM y cotejarlo antes de preparar las matrices. El ajuste de los modelos usa estas matrices. Las imputaciones se escriben en un directorio de salida independiente.

Para auditar el contenido desde la raíz del repositorio:

```bash
python -m auditoria.verificar
sha256sum -c auditoria/integridad.sha256
```

El primer comando abre ambas bases en modo de lectura y compara identificadores, etiquetas de ausencia, coordenadas, valores, incertidumbres y factores de conversión. El segundo comprueba los bytes de los ocho archivos de trabajo, las siete fuentes y las cuatro tablas documentales. El [cuaderno de procedencia de las tablas](https://github.com/pedro1ani/Ipre/blob/793cba71d181ae1ba4748e39652bf098efc71db7/historico/01_auditoria_preparacion.ipynb) conserva las decisiones por fila y su exportación.
