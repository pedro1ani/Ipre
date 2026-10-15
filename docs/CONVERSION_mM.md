# Paso 1: fuentes y concentraciones en mM

Se reconstruyó una capa de concentraciones a partir de la base auditada v3 y se contrastaron todos sus valores numéricos originales de metabolitos con sus celdas en los libros Ishii, MMC2 y MMC3. Las fuentes permanecen intactas. No se incorporaron resultados imputados ni se recalcularon los modelos anteriores.

## Factor verificado

**V = 2,3 mL/gCDW**, entendido como volumen citosólico intracelular por gramo de biomasa seca.

- Fuente numérica: `sources/mmc2.xlsx`, hoja `Conversion Factors`, celda **B6**.
- Definición: **A5**, “Conversion factor from gCDW to intracellular cytosolic volume”.
- Unidad: **A6**, mL/gCDW.
- Alcance indicado: **B5**, “All conditions”.
- Respaldo metodológico: `sources/Gerosa_methods.pdf`, página PDF **9**, sección “Intracellular metabolite concentrations by LC-MS/MS”. El apartado describe tanto cultivos estacionarios como transición diáuxica y declara el uso del factor de Dataset S1 para convertir µmol/gCDW a µmol/mL.

Esto respalda aplicar el mismo factor nominal a las concentraciones intracelulares de MMC2 y MMC3. No demuestra que el volumen celular real permanezca constante durante la transición. La fuente lo describe como previamente determinado, pero el material revisado no especifica aquí el ensayo original de determinación ni su incertidumbre. La procedencia está documentada hasta el factor publicado por Gerosa; no se inventa una SD del factor.

No confundir 2,3 mL/gCDW con los factores OD600–biomasa de la fila 2 de esa misma hoja, que tienen otra definición y dependen del sustrato. Tampoco se aplica a Ishii el factor de Gerosa.

## Conversión y alcance

`C_Gerosa_mM = C_Gerosa_umol_gCDW / 2.3`

La equivalencia dimensional es µmol/mL = mmol/L = mM. No hace falta multiplicar por 1.000.

`SD_Gerosa_mM = SD_Gerosa_umol_gCDW / 2.3`

La SD convertida está condicionada al factor nominal fijo y no incluye incertidumbre del volumen. Ishii conserva sus concentraciones originales en mM, sin conversión previa a µmol/g. Se conserva cualquier incertidumbre original junto con su tipo y evidencia; no se crea incertidumbre donde no existe.

| Conjunto | Posiciones | Valores conservados/convertidos | Faltantes conservados |
|---|---:|---:|---:|
| Ishii intracelular, sin conversión | 19.686 | 4.813 | 14.873 |
| Ishii extracelular, sin conversión | 264 | 264 | 0 |
| Gerosa MMC2, dividido por 2,3 | 344 | 342 | 2 |
| Gerosa MMC3, dividido por 2,3 | 189 | 189 | 0 |

Las filas extracelulares de Ishii se conservan identificadas y separadas por compartimento. No se mezclan automáticamente con concentraciones intracelulares. Las filas analíticas no equivalen necesariamente a moléculas únicas: se preservan protocolos, mediciones compuestas y estado de homologación.

Los dos faltantes de MMC2 corresponden a cuantificación intracelular no fiable con abundancia externa del mismo metabolito, según la nota original; siguen siendo NULL. Los vacíos Ishii conservan sus códigos previos de ausencia. Un knockout no se convierte en cero.

## Trazabilidad por medición

Cada registro contiene identificador único, conjunto, analito, condición/muestra, tiempo, compartimento, protocolo, estado de homologación, valor y unidad originales, incertidumbre original y tipo, estado de ausencia, archivo/hoja/celda originales y hash del archivo. Los campos añadidos incluyen `value_mM`, `uncertainty_mM`, `conversion_divisor`, `intracellular_volume_mL_gCDW`, `factor_id`, `conversion_formula` y la nota de incertidumbre.

La información heredada de la auditoría se conserva, incluyendo sus ambigüedades pendientes. La conversión de unidades no resuelve por sí sola la comparabilidad biológica, la identidad de compuestos o la causa de todos los faltantes.

## Archivos

- `data/concentraciones_mM.db`: nueva base derivada de concentraciones. Tabla principal `metabolite_measurements_mM`, tabla `conversion_factors`, vistas `intracellular_mM` y `extracellular_mM`.
- `data/ISHII2007_mM.csv`, `MMC2_mM.csv`, `MMC3_mM.csv`: exportaciones con originales y conversiones, sin ocultar faltantes.
- `data/conversion_factors.csv`: definición, valor, alcance y procedencia del factor.
- `data/factor_evidence.txt`: evidencia extraída de la hoja y página metodológica.
- `data/summary.csv`: conteos verificados.
- `sources/`: libros y métodos originales más la base auditada usada como entrada. Esta base se preserva únicamente como fuente trazable, sin ejecutar sus métodos o reutilizar imputaciones previas.
- `convertir_unidades.py`: reconstrucción reproducible desde las fuentes.
- `manifest.json`: hashes y controles de integridad y conversión.

Para reconstruir, instalar pandas, numpy, openpyxl y pypdf en un entorno Python y ejecutar `python convertir_unidades.py` en una copia con la carpeta `data` vacía. El programa protege la base de salida contra sobrescritura. Las versiones exactas del entorno ejecutado se incluyen en `environment.json`.

## Controles y siguiente paso

Se verificaron valores numéricos contra celdas originales, SD de Gerosa, unidades, unicidad, conservación de faltantes, identidad de Ishii, reversibilidad de conversión y la integridad SQLite. Se comprobaron hashes para confirmar que las fuentes no cambiaron.

Este paso no normaliza por glucosa ni imputa. Cuando se retome la normalización, debe calcularse desde estas concentraciones en mM. Los cocientes históricos Ishii/Gerosa calculados convirtiendo Ishii con su propio factor de volumen no deben copiarse como resultados de este nuevo esquema: el cambio de base puede modificar esos cocientes. Para Gerosa/Gerosa, el factor común 2,3 se cancela en el cociente.
