IPre | Imputación de concentraciones de metabolitos en E. coli
1. Descripción del proyecto
En este proyecto evaluamos métodos de imputación para estimar concentraciones faltantes de metabolitos en Escherichia coli, utilizando datos experimentales de Ishii y Gerosa.

Comparamos tres métodos de scikit-learn:
- Mediana: utiliza la mediana de las concentraciones observadas de cada metabolito.
- KNNImputer: estima valores a partir de muestras con perfiles metabólicos similares.
- IterativeImputer: estima cada metabolito mediante regresiones que utilizan las concentraciones de otros metabolitos. Se emplea BayesianRidge como estimador.
Evaluamos el error, la cobertura y la incertidumbre de las estimaciones. También analizamos si incorporar datos de Gerosa MMC2 al entrenamiento modifica el desempeño de los métodos sobre las mediciones de Ishii. Los valores faltantes se completan únicamente cuando cumplen los criterios de elegibilidad y evaluación establecidos.
2. Organización del proyecto
El análisis se desarrolla en dos notebooks, que se ejecutan en el siguiente orden.
01_auditoria_preparacion.ipynb
Prepara las bases experimentales y construye las matrices utilizadas en la evaluación. Incluye:
- Revisión de las fuentes de Ishii y Gerosa.
- Identificación y homologación de metabolitos.
- Revisión de la disponibilidad de mediciones.
- Conversión de concentraciones a mM.
- Definición de las referencias de normalización.
- Revisión preliminar de elegibilidad por metabolito.
- Construcción de las matrices de entrada.
02_imputacion_evaluacion.ipynb
Evalúa y aplica los métodos de imputación. Incluye:
- Preparación de máscaras para ocultar concentraciones medidas.
- Validación externa e interna.
- Selección de parámetros de KNN e Iterative.
- Cálculo de errores, cobertura y abstenciones.
- Estimación de incertidumbre de IterativeImputer.
- Generación de matrices imputadas.
- Cálculo de estadísticas descriptivas y elaboración de figuras.
3. Datos experimentales
Utilizamos 30 muestras de Ishii correspondientes a una tasa de dilución de $D = 0.2\ \mathrm{h}^{-1}$. Las muestras incluyen diferentes grupos de perturbación metabólica. KO05 y KO05x se consideran cultivos separados dentro de su grupo, mientras que las cinco muestras WT se mantienen agrupadas durante la validación. Se excluyen cuatro muestras correspondientes a otras tasas de dilución.

En el escenario que incorpora Gerosa utilizamos ocho condiciones experimentales admisibles de MMC2. Estas mediciones provienen de otro estudio y de condiciones de cultivo diferentes. Su aporte se evalúa mediante la validación sobre Ishii. Los datos temporales de MMC3 y las anotaciones de MMC4 no se utilizan como muestras de entrenamiento.
Homologación de unidades
Todas las concentraciones utilizadas en el análisis se expresan en mM. Ishii entrega las concentraciones en esas unidades. Las mediciones de Gerosa se convierten desde µmol/gCDW utilizando un volumen intracelular de 2.3 mL/gCDW, documentado en los métodos experimentales.

La conversión corresponde a:
\[
C_{\mathrm{Gerosa}}[\mathrm{mM}] = \frac{C_{\mathrm{Gerosa}}[\mu\mathrm{mol}/\mathrm{gCDW}]}{2.3[\mathrm{mL}/\mathrm{gCDW}]}
\]
También calculamos una concentración normalizada para cada metabolito. Como referencia utilizamos la concentración del mismo metabolito reportada por Gerosa bajo la condición de glucosa, expresada en mM:
\[
R_{ij} = \frac{C_{ij}[\mathrm{mM}]}{C_{\mathrm{ref},j}[\mathrm{mM}]}
\]
En esta expresión, $C_{ij}$ es la concentración del metabolito $j$ en la muestra $i$ y $C_{\mathrm{ref},j}$ es su concentración de referencia en Gerosa. Utilizamos la misma referencia por metabolito en todas las muestras y comprobamos su disponibilidad y positividad. No propagamos la incertidumbre del volumen intracelular ni de las concentraciones de referencia.
Metabolitos considerados
El panel contiene 20 metabolitos:

pyr, fum, succ, akg, pep, dhap, s7p, ade, cit, cAMP, amp, gmp, adp, gdp, atp, gtp, nad, nadh, nadp y nadph.

La disponibilidad experimental varía entre metabolitos. DHAP, NADH y NADPH tienen 8, 6 y 6 observaciones en Ishii, respectivamente. Se conservan sus concentraciones medidas, pero no se completan sus faltantes. Los otros 17 metabolitos cumplen la revisión preliminar de disponibilidad; su participación efectiva depende de los criterios aplicados en cada entrenamiento y de los resultados de evaluación.
4. Métodos de imputación
Mediana
Sustituye un valor faltante por la mediana de las concentraciones observadas del metabolito en el conjunto de entrenamiento. No utiliza las concentraciones de otros metabolitos para calcular la estimación.
KNNImputer
Identifica muestras similares utilizando las concentraciones disponibles y calcula la estimación a partir de las muestras vecinas. Se evalúan cuatro configuraciones:
- 3 vecinos con ponderación uniforme.
- 5 vecinos con ponderación uniforme.
- 7 vecinos con ponderación uniforme.
- 3 vecinos con ponderación por distancia.
IterativeImputer
Estima las concentraciones faltantes mediante regresiones sucesivas que utilizan otros metabolitos como predictores. El estimador es BayesianRidge, que incorpora regularización y permite realizar imputaciones estocásticas mediante muestreo de su distribución posterior.

Se evalúan dos configuraciones:
- max_iter=50, tol=0.001.
- max_iter=100, tol=0.0001.
BayesianRidge utiliza los parámetros predeterminados de scikit-learn. Las configuraciones de KNN e Iterative se seleccionan mediante la validación interna. La mediana no requiere una búsqueda de parámetros.
5. Escenarios de análisis
Se utilizan dos conjuntos de entrenamiento:
- Escenario A — Ishii: incluye únicamente las muestras de Ishii.
- Escenario B — Ishii + MMC2: incorpora ocho condiciones experimentales admisibles de Gerosa MMC2.
Los métodos se evalúan sobre las mismas concentraciones observadas de Ishii en ambos escenarios. Para mantener la comparabilidad, A y B utilizan la intersección de metabolitos elegibles en cada entrenamiento.
6. Validación de los métodos
Para evaluar las predicciones, ocultamos concentraciones experimentales conocidas y comparamos los valores estimados con sus mediciones originales.

El protocolo utiliza cinco bloques externos definidos por grupos de perturbación y tres bloques internos para seleccionar parámetros. Cada consulta oculta una concentración individual y mantiene visibles las demás mediciones disponibles de la muestra. El grupo experimental de esa muestra se excluye del entrenamiento.

Las 457 concentraciones documentadas de Ishii definen los intentos de evaluación comunes. Los valores faltantes originales no se utilizan como valores de referencia para calcular errores.
Criterios de elegibilidad
La elegibilidad se determina dentro de cada conjunto de entrenamiento. Para incluir un metabolito se exige:
- Al menos 10 mediciones disponibles.
- Mediciones distribuidas en al menos 5 grupos experimentales.
- Al menos 2 concentraciones distintas.
- Al menos 2 predictores con 5 pares de mediciones disponibles.
KNN e Iterative requieren, además, al menos dos predictores medidos en la muestra evaluada. KNN necesita al menos tres muestras donantes con una distancia calculable. Si no existe soporte suficiente, se registra una abstención y no se genera una predicción.

El escalamiento se ajusta exclusivamente con los datos de entrenamiento. No se aplica transformación logarítmica. La selección interna utiliza el MAE normalizado macro, calculado sobre las mismas celdas entre las configuraciones evaluadas.
Métricas de evaluación
- MAE en mM: error absoluto medio entre concentraciones estimadas y observadas.
- MAE normalizado (MAE_R): error absoluto medio calculado sobre las concentraciones normalizadas.
- Cobertura: proporción de intentos de evaluación que generan una predicción.
- Número de predicciones: cantidad de concentraciones efectivamente estimadas y evaluadas.
La comparación principal utiliza las mismas celdas para los tres métodos y ambos escenarios. El MAE macro se obtiene calculando primero el error por metabolito y después promediándolo entre metabolitos. MAE_R no corresponde a un porcentaje de error.
7. Estimación de incertidumbre
La incertidumbre de IterativeImputer se analiza mediante 30 imputaciones estocásticas por escenario y bloque externo, además de la aplicación final. Se utiliza sample_posterior=True con semillas distintas, desde 1100 hasta 1129. En esta etapa se conservan los parámetros seleccionados durante la evaluación determinista.

Para cada concentración estimada se guardan:
- Los 30 valores generados.
- La media y la desviación estándar muestral.
- Los cuantiles de 2.5 % y 97.5 %.
- La amplitud del intervalo entre cuantiles.
- La inclusión o exclusión del valor experimental en el intervalo de validación.
La desviación estándar describe la variación de la estimación entre las ejecuciones estocásticas. Los intervalos dependen de los supuestos del modelo y no incluyen todas las fuentes de incertidumbre experimental y biológica. Una dispersión baja no garantiza un error bajo.
8. Aplicación de la imputación
Aplicamos los métodos evaluados a las matrices originales. Para incorporar una estimación a la matriz final, el metabolito debe cumplir los criterios de elegibilidad y contar con al menos 10 predicciones externas distribuidas en 3 bloques.

Se genera una matriz por método y escenario, para un total de seis matrices de 30 muestras por 20 metabolitos. Todas conservan las concentraciones medidas. Las posiciones que no cumplen los criterios permanecen sin imputar y se registra el motivo.

Las matrices resultantes se utilizan para análisis exploratorios. La elegibilidad y el soporte de evaluación no equivalen a una garantía de exactitud biológica.
9. Resultados del análisis
Los resultados corresponden a la ejecución del 9 de noviembre de 2026 y se encuentran en evaluacion/.

La comparación principal incluye 416 concentraciones predichas en común, de 457 intentos de evaluación, con una cobertura de 91.0 %. Se obtienen predicciones externas para 16 metabolitos.

El MAE normalizado macro es:
Método	Ishii (A)	Ishii + MMC2 (B)
Mediana	0.6587	0.7009
KNN	0.4801	0.4835
Iterative	0.4713	0.4956


Iterative en el escenario A presenta el menor MAE normalizado macro, seguido de KNN en el mismo escenario. La diferencia entre ambos no demuestra superioridad estadística. Al incorporar MMC2, el MAE normalizado macro aumenta en los tres métodos. Los resultados también deben revisarse en mM y por metabolito, porque las conclusiones pueden depender de la métrica utilizada.
Concentraciones imputadas
Las matrices originales contienen 457 mediciones y 143 posiciones faltantes. Cada matriz final completa 36 concentraciones y conserva 107 posiciones faltantes.
Metabolito	Concentraciones imputadas
Piruvato (pyr)	16
Fumarato (fum)	7
Succinato (succ)	3
PEP (pep)	8
S7P (s7p)	1
cAMP	1
Total	36


De los 107 valores pendientes, 70 no cumplen las condiciones de elegibilidad en el entrenamiento y 37 no cuentan con evaluación externa suficiente.
Incertidumbre
Las imputaciones estocásticas generan 416 intervalos externos por escenario. Los intervalos nominales del 95 % contienen 371 de 416 valores experimentales (89.2 %) en el escenario A y 361 de 416 (86.8 %) en el escenario B. Esta cobertura observada es inferior al nivel nominal.

Durante la selección interna determinista se registraron tres advertencias de convergencia de IterativeImputer en evaluacion/ajustes.csv. Se conservan para revisar su posible efecto en la selección de parámetros.
10. Estadísticas y visualización
Calculamos estadísticas descriptivas de las concentraciones en mM, considerando valores observados, imputados y pendientes:
- Media y mediana.
- Desviación estándar.
- Rango intercuartílico (IQR).
- Coeficiente de variación (CV).
- Cantidad de mediciones e imputaciones.
También se calculan estadísticas sobre los mismos faltantes completados por los seis análisis. El CV se deja indefinido cuando hay menos de dos observaciones o la media no es positiva o está próxima a cero.

Seleccionamos cinco metabolitos para revisar situaciones específicas:
- Succinato (succ): caso con buen desempeño.
- Citrato (cit): caso con cobertura limitada.
- Adenina (ade): caso utilizado para analizar incertidumbre.
- NADP (nadp): caso con error elevado.
- DHAP (dhap): metabolito excluido de la imputación por disponibilidad insuficiente.
En el caso del citrato solo existen dos predicciones externas, por lo que su error no permite establecer una conclusión general sobre el desempeño del método.
11. Archivos y resultados
sources/
Contiene las fuentes experimentales, los archivos de Ishii, MMC2–MMC4 y la documentación de métodos y unidades.
auditoria/
Contiene la revisión de identidad de metabolitos, las referencias de normalización, la disponibilidad, los criterios preliminares de elegibilidad y las matrices preparadas.
evaluacion/
Contiene las métricas, los resultados de validación, las imputaciones, las estadísticas y las figuras.
Archivo	Contenido
resumen_global.csv	Comparación de errores y cobertura.
metricas.csv	Métricas por metabolito y bloque.
validacion_celdas.csv	Concentraciones observadas y predicciones externas.
busqueda_interna.csv	Resultados de selección de parámetros.
parametros_seleccionados.csv	Configuraciones seleccionadas.
ajustes.csv	Ajustes realizados y advertencias registradas.
elegibilidad_entrenamientos.csv	Elegibilidad en cada entrenamiento.
sorteos_posteriores.csv.gz	Valores de las imputaciones estocásticas.
incertidumbre.csv	Resúmenes de incertidumbre por concentración.
cobertura_intervalos.csv	Cobertura de los intervalos estimados.
imputaciones.csv	Valores estimados y decisiones de incorporación.
soporte_aplicacion.csv	Soporte de las imputaciones finales.
estadisticos.csv	Estadísticas descriptivas.
casos_representativos.csv	Metabolitos seleccionados para las figuras.
mascaras.json	Máscaras utilizadas en la evaluación.
verificacion.json	Semillas, versiones y verificaciones de la ejecución.


Las seis matrices de salida se identifican mediante archivos matriz_A_*_mM.csv y matriz_B_*_mM.csv. Los gráficos se guardan en evaluacion/figuras/.
12. Cómo ejecutar el proyecto
Para reproducir la ejecución documentada se utiliza Python 3.13.14 y se instalan las dependencias desde la carpeta principal:
pip install -r requirements.txt
Si Jupyter no está instalado en el entorno:
pip install jupyterlab ipykernel
Para abrir los notebooks:
jupyter lab
Ejecutar en orden:
1. 01_auditoria_preparacion.ipynb
2. 02_imputacion_evaluacion.ipynb
Las máscaras existentes se reutilizan únicamente si coinciden con las matrices y la configuración de validación. Para realizar otro estudio, se deben conservar por separado los resultados generados y evitar mezclar archivos de ejecuciones distintas. Los registros de verificación documentan el entorno y las semillas utilizados.
13. Consideraciones del análisis
La evaluación se realiza ocultando mediciones experimentales individuales. Los errores obtenidos no describen necesariamente el comportamiento de los métodos frente a faltantes naturales asociados, por ejemplo, a límites de detección.

Las muestras de Ishii y Gerosa pertenecen a estudios y condiciones diferentes. La incorporación de MMC2 se evalúa como transferencia de información, sin asumir equivalencia experimental.

El MAE debe interpretarse junto con el número de predicciones, la cobertura y los resultados por metabolito. Los bloques externos comparten parte de los datos de entrenamiento, por lo que la desviación estándar entre bloques no constituye un intervalo de confianza.

Los intervalos estocásticos de IterativeImputer reflejan la variabilidad estimada bajo el modelo y no toda la incertidumbre biológica. Las concentraciones imputadas se consideran resultados para análisis exploratorio y no sustituyen las mediciones experimentales.