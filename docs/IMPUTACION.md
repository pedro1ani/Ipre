# Métodos y decisiones conservadas

El ejecutor usa una sola base por ejecución. `Imputador.fit` aprende con entrenamiento y `transform` aplica el modelo a una consulta. Una clase conserva ese estado y permite reutilizar ajustes cuando coinciden las exclusiones. Los dos diagnósticos privados separan las reglas específicas de KNN e Iterative; no calculan estimaciones propias.

| Método | Cálculo |
|---|---|
| SimpleImputer | Mediana del objetivo observado en entrenamiento, sobre R. |
| KNNImputer | Tres vecinos, pesos uniformes y distancia nan_euclidean, sobre R estandarizado. |
| IterativeImputer | Regresiones BayesianRidge, inicialización por mediana y actualizaciones sucesivas sobre R estandarizado. |

La estandarización aprende media y desviación exclusivamente con entrenamiento. No se aplican logaritmos. Se deshace la estandarización y se multiplica por la referencia para volver a mM. Los valores medidos se copian exactamente, evitando redondeos de ida y vuelta.

## Los seis parámetros visibles

| Parámetro | Valor | Motivo |
|---|---:|---|
| semilla | 20261019 | Registrar una ejecución reproducible. |
| vecinos | 3 | Mantener el vecindario del protocolo comparado. |
| rondas_maximas | 50 | Limitar las rondas de Iterative. |
| tolerancia | 0,001 | Conservar su criterio de parada. |
| minimo_donantes | 3 | Regla operativa de soporte. |
| minimo_predictores | 2 | Exigir información observada de la consulta. |

Son decisiones previas del protocolo, no valores optimizados con los errores de evaluación. Las configuraciones duplicadas por estudio se reemplazaron por config_imputacion.json. Los objetivos se leen de config_normalizacion.json.

Se mantienen fijos en el código: BayesianRidge(max_iter=300, tol=0.001), inicialización por mediana, orden ascending, todas las variables activas, skip_complete=False y sample_posterior=False. El límite inferior es R=0, trasladado a la escala z como -media/escala. Los demás valores predeterminados se conservan mediante scikit-learn 1.8.0 fijado en requirements.txt. No se modificaron los hiperparámetros efectivos para simplificar la interfaz.

## Soporte y casos con poca información

- **supported**: KNN requiere al menos tres vecinos con objetivo medido, dos predictores observados compartidos por cada vecino y dos predictores informativos en la consulta. Iterative requiere tres observaciones del objetivo, dos predictores informativos y convergencia.
- **limited_support**: información insuficiente según esos mínimos o falta de convergencia. Se conserva el candidato numérico.
- **training_only**: sustitución general; la mediana siempre pertenece a esta categoría. También se usa cuando la consulta carece de predictores informativos o el objetivo es constante en entrenamiento.
- **abstention**: columna objetivo vacía en entrenamiento o estimación inválida. Se conserva NaN; nunca se rellena cero por falta de datos.

Un predictor informativo está originalmente observado en la consulta y varía en entrenamiento. Los valores rellenados internamente no cuentan como observados. Estas reglas describen información disponible, no garantizan precisión individual.

KNN puede usar menos vecinos que los solicitados. Si no hay ninguna distancia definida, sklearn sustituye por la media de entrenamiento y se registra ese mecanismo. El diagnóstico reproduce la selección de vecinos únicamente para informar sus IDs. Iterative puede utilizar predictores completados internamente y devolver un perfil general incluso para una fila vacía. No hay sustitución silenciosa entre métodos.

Convergencia describe la estabilización de la matriz de entrenamiento, no el acierto de la predicción ni un criterio nuevo durante transform. Se conservan rondas y advertencias. Alcanzar una estimación no convergente no la vuelve supported. Se registran además las iteraciones internas de BayesianRidge.

Los valores negativos de redondeo entre -1e-12 y 0 se corrigen a cero y se marcan. Otros negativos o infinitos producen abstención. El límite no negativo de Iterative puede recortar predicciones; un cero estimado no es evidencia experimental de ausencia del metabolito.

## Estadísticas

`available` acepta todo candidato finito; `supported_plus_baseline` acepta supported de KNN/Iterative y mantiene explícitamente la mediana como referencia. Se comparan celdas propias y comunes a los tres métodos.

`apoyo.describir` calcula media, mediana, SD muestral, IQR lineal y CV. El CV queda sin valor si n<2, media no positiva o magnitud ≤1e-12 mM. El notebook reutiliza ese cálculo para las verdades ocultadas y lee las estadísticas finales del ejecutor. Una caída de dispersión tras imputar no demuestra mejor precisión.
