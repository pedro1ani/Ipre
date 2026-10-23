# Diez metabolitos, parámetros e incertidumbre

## Alcance y cambio respecto de la versión anterior

Esta ampliación añade ATP, ADP y AMP a fumarato, PEP, s7p, succinato, piruvato,
α-cetoglutarato y citrato. Los diez se evalúan como objetivos; solo se imputan
sus ausencias reales. Ishii y MMC2 se analizan por separado. Se mantienen mM,
ratios respecto de glucosa Gerosa y el escalado aprendido con entrenamiento.

La evaluación nueva es **validación cruzada anidada por grupos**: cinco bloques
externos, y tres bloques internos dentro de cada entrenamiento. Los grupos de
perturbación se mantienen íntegros, incluidos KO05/KO05x y WT. En MMC2 cada
condición constituye un grupo. La fila glucosa no se enmascara para evaluación,
porque es la referencia conocida. Puede estar en entrenamiento.

En cada consulta de evaluación se oculta un metabolito medido. Los otros faltantes
naturales permanecen ausentes. Se eligen parámetros por el promedio del MAE_R por
metabolito sobre las mismas celdas respondidas por todas las configuraciones.
El nivel externo nunca participa en esta elección ni en el escalado. No se
combinan los resultados internos con los externos para afirmar precisión.

Este protocolo permite evaluar la selección de parámetros sin usar los mismos
casos para elegirlos y juzgarlos. **Cambia el tamaño del entrenamiento respecto
del protocolo anterior**: se excluye todo un bloque de grupos. Por ello no se
atribuye cualquier diferencia de resultados exclusivamente a nuevos parámetros
o metabolitos. No es un experimento causal de ablación.

Las pruebas anteriores por serie analítica y patrones reales siguen disponibles
en `ejecutar_imputacion.py`. Sus resultados de siete metabolitos se conservan en
`resultados/`; no deben confundirse con `resultados_ampliados/`. La evaluación
nueva de una celda oculta no establece precisión ante todos los patrones de
faltantes múltiples. Los candidatos reales sí se completan conjuntamente.

## Homologación y trazabilidad

Se agregan equivalencias explícitas por etiqueta y compartimento:

| Identificador | Ishii, hoja metabolites | Etiqueta original | MMC2 |
|---|---|---|---|
| amp[c] | analyte 538 | AMP | amp |
| adp[c] | analyte 554 | ADP | adp |
| atp[c] | analyte 566 | ATP | atp |

Las 102 mediciones originales de esos tres analitos Ishii y las 24 celdas MMC2
se contrastaron con los Excel incluidos, antes de restringir a la cohorte de
30 muestras. Los valores Ishii se conservan en mM; las referencias MMC2 proceden
de la conversión auditada con 2,3 mL/g seco. Se conservan fuente, hoja, celda,
unidad y UID en `trazabilidad.csv` y las referencias completas en `referencias.csv`.
Las equivalencias nuevas se registran en `homologacion.csv`; la base SQLite
original no se modifica. No se incorporan mediciones compuestas ni isómeros
ambiguos para alcanzar diez columnas. ATP, ADP y AMP tienen datos completos en
la cohorte Ishii; esto no autoriza inventar faltantes reales.

## Por qué estos parámetros

| Decisión | Razón y alcance |
|---|---|
| Mediana | Referencia robusta y simple; no usa el perfil de la consulta. |
| KNN: 3 o 5 vecinos | Sensibilidad local pequeña frente al escaso número de muestras. Cinco suaviza más; tres conserva mayor localización. |
| Pesos uniformes o distancia | Compara promediar vecinos con dar más peso a los más cercanos. |
| BayesianRidge | Regresión lineal regularizada para pocos datos y predictores relacionados; admite desviación predictiva y muestreo posterior. Es también el estimador predeterminado de Iterative. No demuestra que la relación biológica sea lineal. |
| 300 iteraciones, tolerancia 0,001 del regresor | Configuración base de sklearn 1.8; control numérico, no un parámetro biológico. Se compara con 1000 y 0,0001. |
| 50 o 100 rondas de Iterative | Sensibilidad al criterio de convergencia, sin confundirlo con las iteraciones internas de cada regresión. |
| Priors Gamma: cuatro parámetros 1e-6 o 0,01 | Contraste conjunto de sensibilidad de los hiperpriors de precisión del ruido y coeficientes. Cambiarlos juntos no identifica el efecto de cada uno ni equivale a fijar una penalización ridge. |
| Inicialización por mediana | Punto de partida reproducible, coherente con la referencia sencilla. |
| Orden ascending | Empieza por las variables con menos ausencias. Se mantiene fijo para limitar la búsqueda. |
| Todas las variables activas | Solo hay diez metabolitos; no se añade una selección aleatoria de predictores. |
| skip_complete=False | Aprende también modelos para columnas completas en train que pueden estar ocultadas en una consulta. Es relevante para ATP, ADP y AMP. |
| sample_posterior=False en selección y predicción puntual | Compara estimaciones deterministas sin añadir variación Monte Carlo a la selección. |
| sample_posterior=True en incertidumbre | Genera imputaciones múltiples; esta fase se reporta por separado. |
| Límite inferior -media/escala | Es R=0 escrito en escala estandarizada. No significa permitir concentraciones negativas. |
| Sin límite superior específico | No se dispone de un máximo biológico validado; no se inventa. |

La rejilla contiene una mediana, tres KNN y tres Iterative. Las tres opciones de
Iterative son: base, convergencia más exigente y sensibilidad a priors. No es una
búsqueda exhaustiva ni identifica el óptimo universal. Se registra el elegido
por partición en `particiones.json`, todos los errores internos en
`validacion_interna.csv` y el resumen en `hiperparametros.csv`.

La estandarización hace que dividir previamente por una referencia positiva no
cambie z respecto de estandarizar mM directamente (salvo redondeos). Se mantiene
R por continuidad y para comunicar errores relativos a Gerosa. El macro MAE_R
pondera cada metabolito por su referencia: se muestran también los MAE en mM.

## Incertidumbre: qué se calcula y qué no

**Iterative:** 30 reajustes completos con semillas distintas y muestreo posterior,
con el mismo entrenamiento y configuración elegida. Así también varían los
rellenos del entrenamiento. Se guardan las 30 muestras, su media, SD muestral y
cuantiles 2,5–97,5 %. La estimación puntual principal sigue siendo determinista;
no se le atribuyen los resultados de la media posterior sin calcularlos.
Con muestreo posterior sklearn no aplica la parada temprana por tolerancia:
se ejecuta el número fijado de rondas. No se interpreta esto como convergencia
estadística de una cadena.

**Los tres métodos:** rangos empíricos construidos con los residuos internos del
método y configuración elegidos, separados por metabolito. Se requieren al menos
cinco muestras con residuo. Se suman los cuantiles firmados a la predicción y se
trunca en cero el extremo negativo. No son intervalos conformales ni tienen una
garantía del 95 %. La selección del parámetro usa estos mismos datos internos,
por lo que el rango puede ser optimista; por eso se verifica externamente.

Se informan la fracción de valores ocultados dentro de cada intervalo y su
anchura. Un intervalo amplio puede cubrir más y ser poco útil. Los 30 sorteos
son una estimación Monte Carlo inicial: los extremos tienen precisión limitada.
En MMC2 los residuos disponibles por metabolito pueden ser demasiado pocos;
se deja el intervalo vacío en vez de inventar precisión.

No se propaga la incertidumbre del volumen celular, del denominador Gerosa,
ni toda la incertidumbre experimental. La SD de imputación no sustituye la SD
de medición. No se aplica la regla de Rubin a cada concentración aislada: esta
combina estimaciones de un análisis posterior y sus varianzas dentro/entre
imputaciones, que aquí no se ha definido.

## Distribuciones y supuestos

Antes/después comparten mediciones, por lo que `ttest_ind`, Mann–Whitney y el
KS de dos muestras independientes no proporcionan aquí una comparación válida
sin un diseño adicional. Además, una prueba t evalúa medias; Mann–Whitney no es
una prueba general de igualdad de medianas sin supuestos de forma. `kstwo` es
una distribución de referencia del estadístico KS, no la función para comparar
dos muestras.

Se entregan histogramas con bins comunes, media, mediana, SD, IQR, CV y distancia
de Wasserstein. Para los pares medido–estimado externos se añade distancia
máxima entre ECDF con intercambio pareado de etiquetas por grupo, y una prueba
de signos del sesgo por grupo. Se corrigen los p-valores con Holm dentro de
cada familia de 30 comparaciones por base. Para pocos grupos se enumeran todos
los intercambios; en los demás casos se usan 1999 permutaciones y semilla fija.

Estas pruebas son **exploratorias**: la permutación supone intercambiabilidad de
medido y estimado bajo la hipótesis nula; la prueba de signos exige grupos
independientes con probabilidad de error positivo 1/2. Los modelos de los folds
comparten entrenamiento, de modo que la independencia no está garantizada.
No se usa el p-valor como criterio de selección del método, ni se afirma
preservación de distribuciones por falta de significación. La evaluación
principal es el error externo y el diagnóstico de cobertura.

## Referencias técnicas consultadas

- [IterativeImputer, versión 1.8](https://scikit-learn.org/1.8/modules/generated/sklearn.impute.IterativeImputer.html)
- [BayesianRidge, versión 1.8](https://scikit-learn.org/1.8/modules/generated/sklearn.linear_model.BayesianRidge.html)
- [Validación anidada](https://scikit-learn.org/stable/auto_examples/model_selection/plot_nested_cross_validation_iris.html)
- [Permutación y muestras pareadas](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html)
- [KS de dos muestras: independencia](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ks_2samp.html)
- [Wilcoxon: simetría de las diferencias](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html)
