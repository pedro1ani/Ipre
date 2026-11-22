# IPre · Imputación de concentraciones intracelulares en E. coli

**[Resultados y decisiones](ampliacion/INFORME.md)** · [Métricas, captación y flujos](ampliacion/REVISION_METRICAS_FISIOLOGIA.md) · [Diccionario de metabolitos](ampliacion/resultados/diccionario.csv) · [Error por metabolito](ampliacion/resultados/por_metabolito.csv)

| Notebook principal | Qué permite revisar |
| --- | --- |
| [01 · Qué se midió](01_auditoria_preparacion.ipynb) | Reconstruye la matriz, verifica las celdas en los Excel y separa las 20 correspondencias nominales de las siete señales exclusivas de Ishii. |
| [02 · Evaluación e imputación](02_imputacion_evaluacion.ipynb) | Primero oculta mediciones reales, retira sus grupos, compara Mediana, KNN e Iterative y muestra el error. Los bloques 7 y 8 comparan tasas y proteínas; la aplicación a huecos reales es opcional y está desactivada por defecto. |

El análisis principal compara los mismos 30 cultivos de Ishii en dos escenarios: A entrena con Ishii; B añade ocho condiciones de Gerosa MMC2. Las unidades son mM, no se usan logaritmos y los modelos multivariados estandarizan solo el entrenamiento. **El panel ampliado se evalúa como candidato:** los resultados no justifican sustituir por defecto el panel común ni usar MMC2 para todo metabolito.

**[Bloques de trabajo 1 a 4: datos, referencia, número de cultivos y faltantes conjuntos](ampliacion/BLOQUES_1_4.md).** Incluyen un control de ocho filas adicionales de cada fuente, patrones aprendidos solo en Ishii de entrenamiento y una prueba por serie analítica. Desde la raíz, `python -m ampliacion.bloques_1_4` reproduce estos bloques por separado y no completa ningún hueco real.

**[Bloques 5 y 6: elección del método por metabolito e intervalos de predicción](ampliacion/BLOQUES_5_6.md).** El segundo notebook permite revisar ambos resultados sin recalcularlos ni rellenar ausencias. La comparación enfrenta el selector y los tres métodos fijos sobre las mismas mediciones, y comprueba cuántas observaciones caen en los intervalos de 30, 100 y 200 sorteos. Para ejecutar solo estos análisis, desde la raíz:

```bash
python -m ampliacion.selector
python -m ampliacion.incertidumbre --procesos 4
python -m ampliacion.lectura_5_6
```

**[Bloques 7 y 8: tasas y proteínas](ampliacion/BLOQUES_7_8.md).** Compara los mismos cultivos con y sin captación, flujo o dilución. En Ishii contrasta añadir una proteína elegida por disponibilidad con añadir una propuesta por la reacción. El [notebook 02](02_imputacion_evaluacion.ipynb) lee los resultados guardados. Para reproducir la validación y sus intervalos desde la raíz: `python -m ampliacion.bloques_7_8 --procesos 4`. No completa huecos reales.

Desde la raíz del repositorio, con Python 3.12 y un entorno virtual:

```bash
python -m pip install -r ampliacion/requirements.txt
python -m ampliacion.validar
python -m ampliacion.revisar_escalas
python -m unittest discover -s tests
```

También se pueden ejecutar las celdas en orden con ese entorno como kernel de Jupyter. La primera ejecución utiliza las fuentes y las máscaras ya versionadas. En el segundo notebook la sección de aplicación puede volver a ejecutarse desde un kernel nuevo tras la celda de importación; comprueba que los datos y la validación coinciden. Los ocho auxiliares nunca se usan como verdades.

**Lectura del código.** [preparar.py](ampliacion/preparar.py) identifica las señales y comprueba su procedencia. [modelos.py](ampliacion/modelos.py) contiene un ajuste y explicita cuándo hay datos suficientes. [validar.py](ampliacion/validar.py) separa los grupos y las comparaciones controladas. [resumir.py](ampliacion/resumir.py) calcula los errores y aplica estimaciones donde la validación tiene suficiente cobertura. Los CSV, los parámetros, el registro de ejecución y el gráfico están en [ampliacion/resultados](ampliacion/resultados).

**Histórico del 16 de noviembre.** Los [dos notebooks y su documentación](historico/) se conservan íntegros. Para reproducir esa auditoría SQLite/Excel y sus 30 sorteos posteriores por ajuste use el directorio raíz como directorio de trabajo, las dependencias del [estudio histórico](requirements.txt) y, preferiblemente, un checkout separado de [3c6ab55](https://github.com/pedro1ani/Ipre/tree/3c6ab55a2bfcfc819b475916d530384e9cb78fb2). Sus resultados permanecen en [evaluacion](evaluacion); los intervalos históricos no describen la incertidumbre de las nuevas predicciones.
