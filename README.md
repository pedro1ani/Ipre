# IPre · Imputación de concentraciones de E. coli

Comparar Mediana, KNN e Iterative de scikit-learn para reconstruir concentraciones de Ishii, con y sin información de Gerosa MMC2, sin forzar la imputación de todos los faltantes.

## Dos notebooks, en este orden

1. **01_auditoria_preparacion.ipynb**: fuentes, identidad de metabolitos, disponibilidad, conversión a mM, referencia de glucosa, elegibilidad preliminar y matrices.
2. **02_imputacion_evaluacion.ipynb**: máscaras, validación anidada, selección de parámetros, error/cobertura, incertidumbre, aplicación final y cinco figuras.

No hay ejecutores externos, clases personalizadas, KNN propio ni PCA. Las pocas funciones del segundo notebook representan operaciones repetidas: elegibilidad, máscaras, ajuste/predicción y selección interna. La regla de elegibilidad se conserva literalmente en ambos notebooks para que puedan leerse sin cargar módulos auxiliares.

## Ejecutar

Con Python 3.12, crear un entorno e instalar `pip install -r requirements.txt`. Abrir Jupyter desde esta carpeta y ejecutar cada notebook de arriba abajo. Si no se dispone de Jupyter, instalar `jupyterlab ipykernel` en ese entorno y abrir `jupyter lab`. Las versiones numéricas efectivamente utilizadas están fijadas en `requirements.txt`; la interfaz Jupyter no interviene en los cálculos. Se incluyen las salidas ejecutadas mediante IPython. No es necesario ejecutar los scripts antiguos.

Las máscaras existentes solo se reutilizan si coinciden exactamente con los datos y el protocolo. Para un estudio deliberadamente distinto, guardar primero sus resultados y retirar `evaluacion/`; no mezclar ejecuciones con matrices diferentes. La ejecución completa incluye 360 ajustes posteriores y puede tardar varios minutos.

## Datos y reglas científicas

- `sources/`: Excel Ishii, MMC2–MMC4, métodos de Ishii/Gerosa y base auditada. `data/concentraciones_mM.db` se conserva únicamente como comparación de la conversión previa, requerida por el notebook 01.
- 30 muestras Ishii a D = 0.2 h⁻¹. KO05 y KO05x son cultivos separados dentro del mismo grupo; los cinco WT se mantienen juntos al validar. Se excluyen cuatro muestras de otras diluciones.
- Ishii ya está en mM. Gerosa se convierte desde µmol/gCDW dividiendo por **2.3 mL/gCDW**, documentado en MMC2 y métodos. No se utiliza 12.3.
- R = concentración mM / concentración Gerosa bajo glucosa mM, con la misma referencia por compuesto para todas las muestras. Se comprueban positividad, medición y SD documentada. No se propaga la incertidumbre de referencia ni del volumen celular.
- Panel documental de 20: pyr, fum, succ, akg, pep, dhap, s7p, ade, cit, cAMP, amp, gmp, adp, gdp, atp, gtp, nad, nadh, nadp, nadph. DHAP, NADH y NADPH tienen solo 8, 6 y 6 observaciones Ishii: se conservan sin completar. Los otros 17 tienen soporte preliminar, sujeto a entrenamiento y evaluación.
- A entrena con Ishii; B agrega las ocho condiciones admisibles de MMC2. MMC3 temporal y nodos de MMC4 no se tratan como muestras de entrenamiento. MMC2 es otro estudio/cultivo/carbono: se evalúa transferencia, no equivalencia experimental.

## Validación y aplicación

Se mantienen cinco bloques externos por grupos de perturbación, con tres internos para seleccionar parámetros. Cada consulta oculta una observación y conserva los demás valores visibles de esa muestra; todo su grupo se excluye del entrenamiento. Las 457 mediciones documentales Ishii definen los mismos intentos para los tres métodos y ambos escenarios. Los faltantes originales nunca son verdades de evaluación.

Dentro de cada entrenamiento: ≥10 mediciones, ≥5 grupos, ≥2 valores distintos y ≥2 predictores con ≥5 pares medidos. A/B usan la intersección de elegibles. Escalamiento ajustado solo al entrenamiento, sin logaritmos. Mediana no se ajusta por búsqueda; KNN compara cuatro configuraciones (3/5/7 uniforme y 3 por distancia); Iterative compara 50 rondas/tol 0.001 y 100/tol 0.0001. BayesianRidge mantiene sus parámetros predeterminados, registrados. El criterio interno es MAE_R macro sobre celdas comunes entre configuraciones.

Mediana es una sustitución general del entrenamiento. KNN/Iterative requieren dos predictores medidos en la consulta; KNN además tres donantes con distancia definida. Sin soporte se registra abstención, sin rescate silencioso. Las advertencias de convergencia se conservan.

La aplicación final exige además ≥10 predicciones externas en ≥3 bloques. Es una regla operativa de evaluabilidad, no un certificado de exactitud. No hay un umbral biológico de error aceptable establecido. Las matrices son estimaciones para análisis exploratorio: conservan las 30×20 posiciones y todas las mediciones originales.

## Incertidumbre y estadísticas

Después de la comparación determinista, se ejecutan 30 ajustes Iterative con posterior y semillas 1100–1129 por escenario/bloque externo, y por escenario en la aplicación final. No se repite la búsqueda interna. Se ejecutan explícitamente todas las rondas seleccionadas: no hay parada temprana posterior. En la aplicación, los huecos de cada muestra se generan conjuntamente en cada repetición.

Se guardan los 30 valores, media, SD muestral, cuantiles 2.5/97.5%, amplitud e inclusión de las verdades externas. Los intervalos dependen del modelo; no cuantifican toda la incertidumbre biológica. Baja SD no garantiza precisión. Los candidatos sin evaluación suficiente pueden tener sorteos exploratorios, sin incorporarse a la matriz final.

Los descriptivos en mM incluyen observados/imputados/pendientes, media, mediana, SD, IQR y CV. Se calculan tanto con todos los resultados disponibles como sobre los mismos huecos completados por los seis análisis. El CV queda indefinido si n<2 o la media no es positiva/se aproxima a cero.

## Dónde leer los resultados

- `auditoria/`: trazabilidad, referencias, disponibilidad, elegibilidad y matrices de entrada.
- `evaluacion/resumen_global.csv`, `metricas.csv`, `validacion_celdas.csv`: MAE macro/micro, conteos, cobertura, errores por compuesto/bloque y abstenciones. `comunes_seis` es la comparación principal.
- `busqueda_interna.csv`, `parametros_seleccionados.csv`, `ajustes.csv`, `elegibilidad_entrenamientos.csv`: selección, rondas, advertencias y soporte.
- `sorteos_posteriores.csv.gz`, `incertidumbre.csv`, `cobertura_intervalos.csv`: incertidumbre separada del análisis determinista.
- `matriz_A_*_mM.csv`, `matriz_B_*_mM.csv`, `imputaciones.csv`, `soporte_aplicacion.csv`: seis matrices, candidatos, estimaciones aceptadas y motivos para abstenerse.
- `estadisticos.csv`, `casos_representativos.csv`, `figuras/`: cambios descriptivos y cinco gráficos.
- `mascaras.json`, `verificacion.json`: máscaras, semillas, versiones, huellas y controles.

## Límites de la comparación

Los grupos externos comparten datos de entrenamiento, por lo que la SD entre bloques no es un intervalo de confianza. Algunos metabolitos tienen pocas predicciones: interpretar siempre MAE junto con n y cobertura. Ocultar mediciones detectadas no reproduce necesariamente faltantes por límite de detección. El protocolo actual evalúa ocultamiento individual por perturbación; no presenta como recalculadas las antiguas pruebas de patrones simultáneos o series analíticas.

Los MAE históricos del panel de siete no son directamente comparables: cambian máscaras, panel y reglas de soporte. Para atribuir una diferencia solo a ampliar el panel sería necesario un experimento controlado adicional. No se reutilizan esos MAE como resultados actuales.

## Resultado principal de esta ejecución

| Método | Ishii A | Ishii + MMC2 B |
|---|---:|---:|
| Mediana | 0.6682 | 0.6981 |
| KNN | 0.4964 | 0.4551 |
| Iterative | 0.4272 | 0.4715 |

MAE_R macro sobre 413 celdas comunes (90.4% de 457 intentos). Iterative/Ishii tiene el menor error global; añadir MMC2 favorece KNN pero no Iterative. Se completan 20 huecos y quedan 123 pendientes en cada matriz documental. La cobertura de intervalos es 85.2% en A y 88.4% en B entre intervalos emitidos. La interpretación por metabolito y sus limitaciones se encuentran al final del notebook 02.
