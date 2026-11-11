# IPre | Imputación de concentraciones de metabolitos en *E. coli*

## 1. Descripción del proyecto

Este proyecto busca estudiar distintos métodos de imputación para estimar concentraciones faltantes de metabolitos en *Escherichia coli*, utilizando datos experimentales de Ishii y Gerosa.

La idea es evaluar qué tan bien podemos reconstruir concentraciones que no conocemos a partir de la información disponible de otros metabolitos, y determinar en qué condiciones estas estimaciones son suficientemente respaldadas por los datos.

Para esto, trabajamos con tres métodos de `scikit-learn`:

- **Mediana:** reemplaza los valores faltantes utilizando la mediana de las concentraciones disponibles de cada metabolito.
- **KNNImputer:** estima concentraciones a partir de muestras que presentan perfiles metabólicos similares.
- **IterativeImputer:** utiliza las relaciones entre los metabolitos para estimar sus concentraciones mediante modelos de regresión. En este caso, utilizamos `BayesianRidge`.

También buscamos estudiar si incorporar información experimental de Gerosa permite mejorar las estimaciones obtenidas con los datos de Ishii.

El objetivo no es necesariamente completar todos los valores faltantes, sino identificar cuáles podemos estimar, qué error presentan los métodos y qué tan confiables son estas predicciones.

## 2. Organización del proyecto

El análisis está organizado en dos notebooks principales, que deben ejecutarse en el siguiente orden.

### `01_auditoria_preparacion.ipynb`

Este notebook se encarga de revisar y preparar los datos experimentales.

Sus principales funciones son:

- Revisar las fuentes de Ishii y Gerosa.
- Identificar los metabolitos presentes en las distintas bases de datos.
- Analizar la disponibilidad de mediciones para cada metabolito.
- Homologar las concentraciones a mM.
- Establecer las referencias utilizadas para normalizar las concentraciones.
- Evaluar qué metabolitos cuentan con información suficiente para participar en la imputación.
- Construir las matrices de datos que utilizaremos en el análisis.

### `02_imputacion_evaluacion.ipynb`

Este notebook desarrolla la evaluación y aplicación de los métodos de imputación.

Incluye:

- Generación de máscaras para ocultar concentraciones conocidas.
- Validación de los tres métodos.
- Selección de parámetros para KNN e Iterative.
- Cálculo de errores y cobertura.
- Estimación de incertidumbre de IterativeImputer.
- Imputación de las concentraciones faltantes.
- Análisis estadístico y generación de figuras.

## 3. Datos experimentales

Trabajamos principalmente con concentraciones de metabolitos medidas en *E. coli* a partir de los estudios de Ishii y Gerosa.

La base de Ishii contiene 30 muestras seleccionadas correspondientes a una tasa de dilución de 0.2 h⁻¹, incluyendo diferentes perturbaciones metabólicas.

Para complementar esta información utilizamos las condiciones experimentales admisibles de Gerosa MMC2.

Como ambos estudios utilizan condiciones experimentales diferentes, los datos de Gerosa se incorporan como información complementaria. Su utilidad se determina a partir de los resultados de validación.

### Homologación de unidades

Para poder trabajar con ambas fuentes, expresamos las concentraciones en mM.

Las concentraciones de Ishii ya se encuentran en estas unidades, mientras que las de Gerosa se convierten desde µmol/gCDW utilizando un volumen intracelular de 2.3 mL/gCDW, respaldado por la documentación experimental.

La conversión corresponde a:

\[
C_{\mathrm{Gerosa}}[\mathrm{mM}]
=
\frac{C_{\mathrm{Gerosa}}[\mu\mathrm{mol/gCDW}]}{2.3\ [\mathrm{mL/gCDW}]}
\]

También trabajamos con concentraciones normalizadas respecto a la concentración de referencia de cada metabolito en Gerosa bajo la condición de glucosa:

\[
R_{ij}=\frac{C_{ij}[\mathrm{mM}]}{C_{\mathrm{ref},j}[\mathrm{mM}]}
\]

Esto nos permite analizar los errores tanto en concentraciones absolutas como en valores relativos a una referencia común.

### Metabolitos considerados

El panel de estudio contiene 20 metabolitos:

`pyr`, `fum`, `succ`, `akg`, `pep`, `dhap`, `s7p`, `ade`, `cit`, `cAMP`, `amp`, `gmp`, `adp`, `gdp`, `atp`, `gtp`, `nad`, `nadh`, `nadp` y `nadph`.

No todos presentan la misma cantidad de mediciones. Por esta razón, evaluamos su disponibilidad antes de decidir cuáles pueden participar en los modelos.

En particular, DHAP, NADH y NADPH presentan una cantidad limitada de observaciones en Ishii, por lo que conservamos sus valores experimentales sin completar sus faltantes.

## 4. Métodos de imputación

Evaluamos tres métodos con el objetivo de estudiar diferentes formas de utilizar la información disponible.

**Mediana**

Utiliza la mediana de las concentraciones observadas de cada metabolito dentro del conjunto de entrenamiento. Funciona como un método de referencia, ya que no considera directamente las relaciones entre los metabolitos.

**KNNImputer**

Busca muestras similares a partir de sus metabolitos observados y utiliza sus concentraciones para estimar los valores faltantes.

Evaluamos cuatro configuraciones, con 3, 5 y 7 vecinos con ponderación uniforme, además de 3 vecinos ponderados por distancia.

**IterativeImputer**

Estima las concentraciones de cada metabolito mediante un modelo que utiliza los demás metabolitos como predictores.

Utilizamos `BayesianRidge` como estimador, ya que permite incorporar regularización a los modelos de regresión y también generar estimaciones estocásticas para estudiar su incertidumbre.

Evaluamos dos configuraciones de IterativeImputer:

- 50 iteraciones, con tolerancia de 0.001.
- 100 iteraciones, con tolerancia de 0.0001.

Los parámetros de `BayesianRidge` corresponden a los valores predeterminados de `scikit-learn`.

## 5. Escenarios de análisis

Realizamos la evaluación considerando dos escenarios.

**Escenario A: Ishii**

Los métodos se entrenan utilizando únicamente las concentraciones experimentales de Ishii.

**Escenario B: Ishii + MMC2**

Incorporamos ocho condiciones experimentales admisibles de Gerosa MMC2 al conjunto de entrenamiento.

En ambos escenarios, la evaluación se realiza sobre las mismas concentraciones observadas de Ishii. Esto permite estudiar cómo influye incorporar los datos de Gerosa en el desempeño de los métodos.

## 6. Validación de los métodos

Para evaluar la imputación, ocultamos concentraciones que conocemos experimentalmente y pedimos a cada método que las estime.

Después comparamos las concentraciones estimadas con las originales para calcular el error.

La validación utiliza cinco bloques externos organizados por grupos de perturbación y tres bloques internos para seleccionar parámetros.

En cada evaluación se oculta una concentración, manteniendo visibles los demás metabolitos medidos de esa muestra. El grupo experimental al que pertenece se excluye del entrenamiento.

Esto permite evaluar las predicciones sin utilizar las mediciones del mismo grupo experimental como información de entrenamiento.

### Criterios de elegibilidad

Antes de entrenar los modelos, verificamos que cada metabolito cuente con suficiente información.

Los criterios principales son:

- Al menos 10 mediciones disponibles.
- Observaciones distribuidas en al menos 5 grupos.
- Al menos 2 valores de concentración distintos.
- Al menos 2 metabolitos predictores con 5 pares de mediciones disponibles.

Estos criterios se aplican dentro de cada entrenamiento para evitar utilizar información de las observaciones reservadas para evaluación.

Además, KNN e Iterative requieren al menos dos predictores medidos en la muestra que se desea imputar. KNN también necesita al menos tres muestras donantes con distancias calculables.

Cuando no existe suficiente información, no realizamos la estimación y registramos el motivo.

### Métricas de evaluación

Utilizamos principalmente:

- **MAE en mM:** error absoluto medio entre las concentraciones experimentales y estimadas.
- **MAE_R:** error absoluto medio calculado sobre las concentraciones normalizadas.
- **Cobertura:** proporción de concentraciones evaluables para las cuales se obtiene una predicción.
- **Cantidad de predicciones:** número de observaciones efectivamente evaluadas para cada metabolito y método.

La comparación principal considera las mismas celdas para los tres métodos y ambos escenarios.

Esto es importante porque un método podría presentar un error menor simplemente por haber realizado predicciones sobre un conjunto de observaciones diferente.

## 7. Estimación de incertidumbre

Además del error de las predicciones, nos interesa estudiar cuánto pueden variar las estimaciones de IterativeImputer.

Para esto realizamos 30 imputaciones estocásticas utilizando `BayesianRidge` con `sample_posterior=True` y diferentes semillas.

A partir de estas repeticiones calculamos:

- Media de las concentraciones estimadas.
- Desviación estándar de las imputaciones.
- Cuantiles de 2.5 % y 97.5 %.
- Amplitud de los intervalos.
- Proporción de concentraciones experimentales que quedan dentro de los intervalos estimados.

La desviación estándar permite describir cuánto cambia la estimación de una misma concentración entre las diferentes ejecuciones.

Es importante considerar que esta incertidumbre depende del modelo utilizado. Una desviación estándar pequeña no significa necesariamente que la predicción sea correcta, por lo que también evaluamos el error y la cobertura de los intervalos.

## 8. Aplicación de la imputación

Una vez evaluados los métodos, los aplicamos a las matrices de concentraciones originales.

Para incorporar una estimación a la matriz final exigimos que el metabolito cumpla los criterios de elegibilidad y cuente con al menos 10 predicciones externas distribuidas en 3 bloques de validación.

Las concentraciones experimentales se conservan y los valores que no cumplen estas condiciones permanecen faltantes.

Generamos una matriz para cada combinación de método y escenario, obteniendo seis matrices finales de 30 muestras y 20 metabolitos.

Estas matrices están destinadas al análisis exploratorio y deben interpretarse considerando la disponibilidad experimental, el error y la incertidumbre de las imputaciones.

## 9. Resultados del análisis

Los resultados de la ejecución del 9 de noviembre de 2026 se encuentran guardados en la carpeta `evaluacion/`.

La comparación principal considera 416 concentraciones evaluadas en común, de un total de 457 intentos, correspondientes a una cobertura del 91.0 %.

El MAE_R macro obtenido fue:

| Método | Ishii (A) | Ishii + MMC2 (B) |
|---|---:|---:|
| Mediana | 0.6587 | 0.7009 |
| KNN | 0.4801 | 0.4835 |
| Iterative | **0.4713** | 0.4956 |

Iterative entrenado únicamente con Ishii presenta el menor MAE_R macro, seguido de KNN en el mismo escenario.

Las diferencias entre estos dos métodos son relativamente pequeñas y no permiten establecer por sí solas una superioridad estadística.

Al incorporar MMC2, el MAE_R macro aumenta ligeramente para los tres métodos. Esto muestra que agregar información de otro estudio no necesariamente mejora la calidad de las estimaciones.

### Concentraciones imputadas

Las matrices iniciales contienen 457 concentraciones experimentales y 143 posiciones faltantes.

En cada una de las seis matrices finales se completan 36 concentraciones, mientras que 107 permanecen sin imputar.

Los metabolitos con concentraciones completadas son:

| Metabolito | Concentraciones imputadas |
|---|---:|
| Piruvato | 16 |
| Fumarato | 7 |
| Succinato | 3 |
| PEP | 8 |
| S7P | 1 |
| cAMP | 1 |
| **Total** | **36** |

En cuanto a la incertidumbre, los intervalos nominales de 95 % obtenidos mediante las imputaciones estocásticas contienen el 89.2 % de los valores experimentales evaluados en el escenario A y el 86.8 % en el escenario B.

Esto indica que los intervalos presentan una cobertura inferior a la nominal, por lo que es necesario considerar esta limitación al interpretar la incertidumbre del modelo.

## 10. Estadísticas y visualización

Para analizar cómo cambia la distribución de las concentraciones al incorporar las estimaciones, calculamos distintos estadísticos descriptivos:

- Media y mediana.
- Desviación estándar.
- Rango intercuartílico (IQR).
- Coeficiente de variación (CV).
- Cantidad de valores observados, imputados y pendientes.

También seleccionamos cinco metabolitos representativos para visualizar diferentes situaciones:

- **Succinato:** caso con buen desempeño.
- **Citrato:** caso con cobertura limitada.
- **Adenina:** caso para estudiar incertidumbre.
- **NADP:** caso con error importante.
- **DHAP:** ejemplo de metabolito que no cumple los criterios para imputación.

Estas figuras permiten estudiar el comportamiento de los métodos a nivel de metabolito, complementando las métricas generales.

## 11. Archivos y resultados

Las principales carpetas de trabajo son:

### `sources/`

Contiene las fuentes experimentales y la documentación utilizada para preparar las bases de datos.

### `auditoria/`

Contiene los resultados de la revisión de fuentes, la disponibilidad de los metabolitos, las referencias de normalización, los criterios de elegibilidad y las matrices preparadas.

### `evaluacion/`

Contiene los resultados de imputación y validación.

Entre los archivos principales se encuentran:

| Archivo | Contenido |
|---|---|
| `resumen_global.csv` | Resumen de errores y cobertura de los métodos. |
| `metricas.csv` | Métricas por metabolito y bloque. |
| `validacion_celdas.csv` | Valores observados y predicciones de validación. |
| `busqueda_interna.csv` | Evaluación de configuraciones de parámetros. |
| `parametros_seleccionados.csv` | Parámetros seleccionados para cada análisis. |
| `elegibilidad_entrenamientos.csv` | Disponibilidad y elegibilidad de los metabolitos. |
| `incertidumbre.csv` | Estadísticos de incertidumbre de Iterative. |
| `cobertura_intervalos.csv` | Cobertura de los intervalos estimados. |
| `imputaciones.csv` | Concentraciones imputadas y sus condiciones de aceptación. |
| `estadisticos.csv` | Estadísticos descriptivos de las concentraciones. |
| `casos_representativos.csv` | Metabolitos seleccionados para el análisis gráfico. |
| `mascaras.json` | Máscaras utilizadas durante la validación. |
| `verificacion.json` | Información de versiones, semillas y controles de ejecución. |

Las seis matrices finales se guardan con nombres que comienzan por `matriz_A_` o `matriz_B_`, según el escenario.

La carpeta `figuras/` contiene los gráficos utilizados para interpretar los resultados.

## 12. Cómo ejecutar el proyecto

Para reproducir el análisis se recomienda utilizar Python 3.13.14 y trabajar dentro de un entorno virtual.

Primero, instalar las dependencias:

```bash
pip install -r requirements.txt
```

En caso de no contar con Jupyter, instalar:

```bash
pip install jupyterlab ipykernel
```

Luego abrir Jupyter desde la carpeta principal:

```bash
jupyter lab
```

Ejecutar los notebooks en este orden:

1. `01_auditoria_preparacion.ipynb`
2. `02_imputacion_evaluacion.ipynb`

El primer notebook prepara las matrices de entrada y el segundo realiza la validación, imputación, análisis estadístico y generación de resultados.

Las ejecuciones utilizan semillas y máscaras documentadas para facilitar la reproducibilidad. Los archivos de verificación registran el entorno y las condiciones utilizadas.

## 13. Consideraciones del análisis

Los resultados deben interpretarse teniendo presentes algunas limitaciones.

Primero, la validación evalúa nuestra capacidad de reconstruir concentraciones experimentales que ocultamos deliberadamente. Esto no garantiza que obtengamos el mismo desempeño al estimar concentraciones que originalmente no fueron medidas.

Segundo, los datos de Ishii y Gerosa provienen de condiciones experimentales diferentes. Por lo tanto, incorporar MMC2 puede aportar información útil, pero también introducir diferencias que afecten las predicciones.

Tercero, el desempeño no debe interpretarse únicamente a partir del MAE. También es necesario considerar cuántas concentraciones se evaluaron, qué cobertura se obtuvo y cómo cambia el error entre metabolitos.

Finalmente, los intervalos de IterativeImputer representan incertidumbre bajo los supuestos del modelo, pero no incluyen todas las fuentes de variabilidad experimental y biológica.

Por estas razones, las concentraciones imputadas se consideran **estimaciones para el análisis exploratorio**, y no sustitutos equivalentes de nuevas mediciones experimentales.

El proyecto busca establecer una metodología reproducible para evaluar qué metabolitos pueden imputarse, qué métodos entregan mejores estimaciones y cuáles son las limitaciones de utilizar estas concentraciones en análisis posteriores.