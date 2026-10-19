# Matrices independientes por estudio

| Estudio | Filas × metabolitos | Medidos | Faltantes | Medidos en cuatro objetivos | Faltantes en cuatro objetivos |
|---|---:|---:|---:|---:|---:|
| Ishii | 30 × 7 | 138 | 72 | 101 | 19 |
| MMC2 | 8 × 7 | 54 | 2 | 31 | 1 |

Cada ejecución entrena exclusivamente con una matriz seleccionada: Ishii o MMC2. `normalizacion/` entrega por separado `matriz_ishii_*` y `matriz_mmc2_*` en mM, ratios R y máscara de observado (1 medido, 0 faltante). Las celdas vacías son ausencias, no ceros. Orden: fumarato, PEP, s7p, succinato, piruvato, α-cetoglutarato y citrato. Los primeros cuatro son objetivos actuales.

## Ishii

29 muestras con dilución documentada de 0,2 h⁻¹ y KO05x, conservada para mantener la cohorte previa de 30. KO05x/pfkA_1 tiene mediciones, pero no fila inequívoca en S2; no se inventa su dilución o fecha. Se agrupa con KO05/pfkA_2 únicamente para excluir ambos del entrenamiento al consultar ese grupo. Se excluyen otras cuatro muestras con diferentes diluciones. Los controles WT forman un grupo de exclusión conjunta.

Para citrato se usa la fila analítica ISHII:metabolites:523 (Nucleotide), con mediciones. No se combina con la fila 133 (Anion), sin valores. Las otras seis variables son Anion. `config_normalizacion.json` fija estas identidades.

## MMC2, aplicación independiente

Ocho condiciones estacionarias: acetato, fructosa, galactosa, gluconato, glucosa, glicerol, piruvato y succinato. Identidades resueltas, compartimento intracelular y señales no compuestas. Se conservan dos ausencias por cuantificación no fiable: piruvato en condición piruvato y succinato en condición succinato. Se entrega una ejecución independiente de evaluación e imputación MMC2.

Los metadatos de muestras y condiciones mantienen dataset_id. Compartir el formato y el motor de cálculo permite aplicar la metodología por separado sin suponer equivalencia biológica entre estudios.

## Referencia fija

R_ij = concentración_ij [mM] / concentración_Gerosa_glucosa,j [mM].

Se conserva la misma referencia positiva y medida por metabolito, con SD y coordenadas de origen en `referencias_glucosa.csv`. La normalización de Ishii no introduce filas MMC2 en su entrenamiento. No se divide por el metabolito glucosa.

La SD Gerosa corresponde al factor de volumen nominal 2,3 mL/gCDW, sin incertidumbre del propio factor. Se conservan las incertidumbres de origen, sin propagar SD del ratio: el denominador es compartido y faltan covarianzas. La fila MMC2_GLUCOSE tiene siete unos porque se divide por sí misma. No se debe ocultar ese objetivo para validación mientras el denominador continúa revelando su valor.

No se aplican logaritmos en estas matrices. Volver a mM consiste en multiplicar por la referencia. Los escaladores aprendidos de los imputadores se ajustan exclusivamente con el entrenamiento del estudio seleccionado.

## Reproducir

`python normalizar.py --output normalizacion_nueva` genera una carpeta nueva con ambos juegos de matrices separados, metadatos, trazabilidad y manifiesto. No entrena ni imputa. Las máscaras de observado no son máscaras de validación artificial.
