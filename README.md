# IPre — carpeta de trabajo para normalización e imputación

Estado: matrices en mM y normalizadas construidas; imputación pendiente.

Los resultados actuales están en `normalizacion/`. Consulte `docs/NORMALIZACION.md` para decisiones, conteos y límites. Reproduzca con `python normalizar.py --output normalizacion_nueva` (carpeta nueva). `config_normalizacion.json` identifica las filas analíticas, los objetivos y las condiciones incluidas.
Esta carpeta sustituye al conjunto mezclado de programas anteriores como punto de partida. No contiene sus métodos, configuraciones, matrices normalizadas ni resultados históricos.

## Qué utilizar

- `data/concentraciones_mM.db`: entrada canónica del análisis. Tabla `metabolite_measurements_mM`; vistas `intracellular_mM` y `extracellular_mM`.
- `data/ISHII2007_mM.csv`, `MMC2_mM.csv`, `MMC3_mM.csv`: exportaciones equivalentes para trabajar sin SQL. No se concatenan con la base: representan las mismas mediciones.
- `data/samples_ishii.csv`: metadatos de muestras, genotipo, dilución y series analíticas, exportados de la fuente auditada. No define por sí solo la partición de validación. KO05x tiene mediciones pero no fila en esta tabla de metadatos; se conserva como pendiente explícito, sin inventar una identidad ni fusionarlo con KO05.
- `data/conversion_factors.csv`, `factor_evidence.txt`, `summary.csv`: factor, evidencia y conteos.
- `sources/`: Ishii, MMC2, MMC3 y MMC4 originales, métodos de Ishii/Gerosa y una única base auditada de origen.
- `docs/CONVERSION_mM.md`: justificación y límites de la conversión.
- `convertir_unidades.py`: reconstrucción de la capa en mM desde la fuente auditada. No ejecuta métodos antiguos.
- `verificar_proyecto.py`: comprobación de integridad del paquete, conversiones y metadatos; solo lectura.
- `requirements.txt` y `environment.json`: versiones del entorno original de conversión. No constituyen todavía el entorno de imputación.
- `manifest.json`: manifiesto original de la conversión, conservado como evidencia de esa etapa.
- `MANIFEST_PROYECTO.json`: inventario y hashes de esta carpeta reorganizada.

La base `sources/ecoli_auditada_v3.db` se conserva íntegra porque es la entrada de `convertir_unidades.py` y la fuente de los metadatos. Incluye tablas fuera del alcance actual; no se usan como matrices de entrenamiento ni se importan resultados antiguos al nuevo análisis. Recortarla impediría mantener su identidad y hash de procedencia.

## Comprobar la carpeta

Extraer el ZIP completo. En una terminal situada dentro de la carpeta:

```powershell
python verificar_proyecto.py
```

Esta comprobación usa únicamente la biblioteca estándar de Python y no requiere instalar dependencias.

Para reconstruir la conversión, trabajar en una copia separada de la carpeta, instalar `python -m pip install -r requirements.txt`, y dejar `data/` vacía antes de ejecutar `python convertir_unidades.py`. No usar Python con `-O`: el convertidor original usa aserciones. Esa reconstrucción regenera los archivos de conversión, pero no `samples_ishii.csv` ni el manifiesto de esta distribución; estos deben exportarse/actualizarse antes de entregar otra versión. No hace falta reconstruir para continuar: las entradas actuales ya están disponibles y verificadas.

## Criterios del flujo (preparación y normalización ya ejecutadas)

1. Seleccionar los metabolitos intracelulares y documentar identidad, señales compuestas y elegibilidad. Mantener los faltantes como ausencias; un knockout no implica concentración cero.
2. Normalizar cada metabolito con su propia concentración Gerosa de referencia en mM:
   `ratio = concentración_mM / referencia_del_mismo_metabolito_mM`.
   La referencia inicial prevista es la condición de crecimiento con glucosa de MMC2; no es dividir todos los metabolitos por la concentración del metabolito glucosa. No usar referencias nulas, cero o ambiguas. Cualquier cambio de condición de referencia debe quedar explícito.
3. Definir cohorte, grupos de exclusión y tratamiento de MMC3 antes de entrenar. La conversión no justifica mezclar automáticamente estudios o tiempos. KO05 y KO05x no se fusionan.
4. Implementar el nuevo flujo con `KNNImputer` e `IterativeImputer`, registrar parámetros y versiones, y ajustar las transformaciones solo con entrenamiento. No se incluye aquí código antiguo de KNN/PCA ni se presupone su reutilización.
5. Fijar semilla, guardar máscaras comunes y particiones; evaluar error y cobertura sobre valores observados ocultados, sin filtración hacia entrenamiento o referencias.
6. Comparar media, mediana, desviación estándar, IQR y CV antes/después por metabolito. Identificar siempre observaciones e imputaciones y devolver estimaciones a mM para interpretación.

La preparación se ejecuta con `normalizar.py` y `config_normalizacion.json`. No hay todavía un ejecutor de imputación; se creará con el nuevo flujo. No es necesario recuperar los paquetes históricos para proceder.

## Alcance y pendientes

Ishii intracelular: 19.686 posiciones, 4.813 valores y 14.873 faltantes. Ishii extracelular: 264 valores. MMC2: 344 posiciones, 342 valores y 2 faltantes. MMC3: 189 valores. Ninguno está imputado en la capa en mM.

Ishii se mantiene en mM. Gerosa se convierte mediante `valor_original / 2,3` (mL/gCDW). La SD se divide por el mismo factor nominal, sin añadir incertidumbre no reportada del volumen.

La homologación incompleta de Ishii, las señales compuestas, la identidad/documentación de KO05x y la falta de límites de detección individuales siguen abiertas. Los campos de trazabilidad y homologación se preservan dentro de la base en mM. MMC4 se conserva para contrastar identidades y la red de referencia; no obliga a usar predictores de red.

Se excluyeron los lanzadores antiguos, scripts de auditoría/benchmark anteriores, configuraciones incompatibles, exportadores redundantes, journal SQLite, CSV de proteínas y resultados históricos. Los archivos originales entregados por el usuario no se modificaron.

## Exploración estadística

`descriptivos/` contiene estadísticos en mM y ratios, cuatro gráficos PNG, un PDF conjunto y su interpretación. Se usan solo valores observados. Ejecutar `python analizar_descriptivos.py --output descriptivos_nuevos`; dependencias en `requirements_descriptivos.txt`. No se han imputado datos.
