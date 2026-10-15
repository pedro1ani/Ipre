# Matrices normalizadas

## Resultado

| Escenario | Filas × metabolitos | Medidos | Faltantes | Medidos en cuatro objetivos | Faltantes en cuatro objetivos |
|---|---:|---:|---:|---:|---:|
| Ishii | 30 × 7 | 138 | 72 | 101 | 19 |
| Ishii + MMC2 | 38 × 7 | 192 | 74 | 132 | 20 |

Orden de columnas: fumarato, PEP, s7p, succinato, piruvato, α-cetoglutarato y citrato. Los cuatro primeros son objetivos. `sample_id` identifica filas y nunca debe utilizarse como predictor numérico.

`normalizacion/` contiene, para cada escenario, una matriz `*_mM.csv`, una matriz de cocientes `*_R.csv` y una máscara `*_observado.csv`. Las celdas vacías representan ausencias; la máscara indica 1 observado y 0 ausente. Los ceros numéricos reales, si existen, no se transforman en ausencias. Son matrices de entrada candidatas; no son modelos entrenados ni particiones de validación.

## Cohorte Ishii

Se conservan las 29 muestras metabolómicas con dilución documentada de 0,2 h⁻¹ más KO05x para mantener el panel previo de 30 muestras. KO05x es pfkA_1 en la columna F del libro original, pero no tiene fila inequívoca en la tabla S2 de metadatos. Su dilución se deja vacía; no se inventa una fecha ni se fusiona con KO05/pfkA_2. Ambas se etiquetan en el grupo pfkA para una futura exclusión conjunta. Los controles WT comparten grupo WT. Cuatro muestras con otras diluciones se excluyen y quedan registradas.

Se seleccionan explícitamente las filas analíticas en `config_normalizacion.json`. Para citrato se usa `ISHII:metabolites:523`, protocolo Nucleotide, que tiene mediciones. La fila `ISHII:metabolites:133`, protocolo Anion y sin valores, no se combina con ella. Los demás seis analitos usan Anion. La selección mantiene identidades y protocolos; no estima equivalencias analíticas entre métodos.

## Admisibilidad MMC2

Se incluyen ocho condiciones estacionarias (acetato, fructosa, galactosa, gluconato, glucosa, glicerol, piruvato y succinato). Sus siete posiciones tienen identidad resuelta, compartimento intracelular, medición no compuesta y referencia positiva. Se acepta una condición incompleta: no se exige completar artificialmente su panel.

Piruvato en la condición piruvato y succinato en la condición succinato permanecen ausentes por cuantificación no fiable. Las otras seis condiciones tienen los siete valores. Se incluyen 54 valores MMC2 y dos ausencias. `condiciones_mmc2.csv` registra sustrato, tasa de crecimiento, SD de tasa, cobertura y decisión.

La admisión significa que pueden ser datos auxiliares para comparar los dos escenarios. No demuestra intercambiabilidad biológica entre cultivos Ishii y Gerosa: las tasas de crecimiento y condiciones difieren. La evaluación primaria posterior deberá ocultar las mismas celdas de Ishii y usar las mismas exclusiones en ambos escenarios; MMC2 se añade solamente al entrenamiento. No se ha ejecutado esa evaluación.

## Referencias y cálculo

Para cada metabolito j y cada muestra i:

R_ij = C_ij [mM] / C_Gerosa_glucosa,j [mM].

Los siete denominadores provienen de MMC2_GLUCOSE. Son positivos, reportados, no imputados, no compuestos y tienen SD y coordenadas de origen documentadas en `referencias_glucosa.csv`. Cada denominador es idéntico en ambos escenarios y para todas las muestras del mismo metabolito. No se divide por la concentración del metabolito glucosa.

Se conservan las concentraciones en mM junto con los ratios, sin aplicar todavía logaritmo ni estandarización. Deshacer la normalización consiste en multiplicar el ratio por su referencia. La fila MMC2_GLUCOSE tiene siete ratios iguales a 1 porque sus numeradores son las propias referencias; no representa siete evidencias independientes del denominador.

La SD de referencia está expresada en mM y condicionada al volumen nominal fijo 2,3 mL/gCDW. No incluye incertidumbre del volumen. Se conservan las incertidumbres de origen, pero no se calcula una SD del ratio: la referencia es compartida entre filas y no se dispone de toda la información de incertidumbre/covarianza. En particular, no debe aplicarse una fórmula de errores independientes al cociente de la referencia consigo misma.

## Evitar filtración en el análisis posterior

No validar MMC2_GLUCOSE ocultando sus valores mientras permanecen disponibles en los denominadores: revelaría el objetivo. En la comparación primaria propuesta se validan únicamente valores de Ishii; las referencias externas se mantienen fijas. La fila glucosa queda marcada como referencia en los metadatos. El efecto de incluirla como donante puede estudiarse posteriormente, sin presentarlo como validación independiente.

Todos los escalados, transformaciones aprendidas, selección de vecinos y ajustes de imputadores deben hacerse dentro del entrenamiento de cada partición. Excluir el grupo de perturbación retenido de Ishii, incluyendo conjuntamente KO05/KO05x. Las máscaras de observado actuales no son máscaras aleatorias de validación; estas aún deben generarse con semilla fija y reutilizarse entre métodos y escenarios.

## Reproducir

Desde la carpeta raíz:

    python normalizar.py --output normalizacion_nueva

La carpeta destino debe ser nueva. El script usa solo la biblioteca estándar de Python. Las entradas se abren en solo lectura, se validan las referencias y unicidad, se comprueba reversibilidad y se comparan hashes de entrada antes y después. `manifest_normalizacion.json` registra las entradas, el código y los hashes de las salidas. Dos ejecuciones se contrastaron byte a byte durante la preparación de esta entrega. No se imputó ningún valor.
