# Tres métodos de imputación y diagnóstico de soporte

Implementación: scikit-learn 1.8.0. Se usan exclusivamente SimpleImputer, KNNImputer e IterativeImputer. No hay un estimador KNN propio ni PCA en el código activo. El cálculo de distancias en el diagnóstico de KNN solo identifica los donantes usados por sklearn; no sustituye su estimación.

## Alcance y ejecución

Único ejecutor: ejecutar_imputacion.py --dataset ishii|mmc2. Primero evalúa mediciones ocultadas y después estima candidatos reales. Ambos bloques usan AuditedImputer de imputadores.py. Cada entrenamiento se obtiene exclusivamente de la matriz seleccionada; no hay concatenación de estudios.

Para Ishii se excluye el grupo completo (incluye KO05/KO05x y los controles WT conjuntos); para MMC2 se excluye la condición consultada. Los cuatro objetivos son fumarato, PEP, s7p y succinato; piruvato, akg y citrato actúan como predictores. Cambiar el estudio selecciona su configuración independiente.

## Parámetros y escalas

| Método | Parámetros y procesamiento |
|---|---|
| SimpleImputer | strategy=median; ratios R originales, sin transformar. Mediana calculada únicamente con valores observados del entrenamiento. |
| KNNImputer | n_neighbors=3, weights=uniform, metric=nan_euclidean. StandardScaler sobre R ajustado solo con entrenamiento; se conservan NaN. No se aplica logaritmo. |
| IterativeImputer | BayesianRidge(max_iter=300, tol=0,001), max_iter=50 rondas, tol=0,001, initial_strategy=median, imputation_order=ascending, n_nearest_features=None, skip_complete=False, sample_posterior=False, random_state=20261019. StandardScaler sobre R ajustado solo con entrenamiento, sin logaritmo. |

Se fijan add_indicator=False y keep_empty_features=False. Antes de ajustar se excluyen explícitamente las columnas completamente vacías del entrenamiento y se restaura su posición con ausencia: nunca se inventa un cero. Las mediciones originales de consulta se copian exactamente, evitando cambios por redondeo al escalar y desescalar.

Para Iterative, el límite inferior es cero en R, trasladado a la escala estandarizada como `-media_train / escala_train` por variable; el superior es infinito. sklearn recorta sus predicciones a ese límite. Esto es una restricción del modelo, no evidencia experimental de un cero. No se cuantifica cuántos recortes internos ocurrieron; los resultados cero deben interpretarse con esa limitación. Al desescalar solo se corrigen a cero negativos de redondeo en [-1e-12, 0), marcados en la salida. Otros negativos o valores no finitos provocan abstención por celda.

La predicción se devuelve a R y luego a mM multiplicando por la referencia fija Gerosa del metabolito. No se usan las unidades ni transformaciones de los antiguos métodos. La estandarización controla escalas numéricas; no corrige diferencias biológicas entre Ishii y MMC2.

## Qué significa soporte

| Etiqueta | Interpretación y uso |
|---|---|
| supported | Disponibilidad mínima de información observada de la muestra según las reglas siguientes. Es soporte operativo, no precisión validada ni intervalo de confianza. |
| limited_support | Hay información parcial, pocos donantes o falta convergencia. Se conserva el candidato, sin aceptarlo en la matriz de candidatos con soporte. |
| training_only | Sustitución basada en la distribución o modelo general del entrenamiento, sin soporte individual suficiente para distinguir la muestra. Se conserva como referencia explícita. |
| abstention | No se puede producir un candidato válido, por ejemplo una columna objetivo completamente vacía. El valor sigue ausente. |

SimpleImputer siempre es training_only: incluso con muchos datos estima la mediana general del objetivo, sin condicionar en los otros metabolitos de la muestra. Esto es el propósito de la referencia básica, no un fallo.

KNN obtiene supported solo si utiliza al menos tres donantes con objetivo medido, cada uno comparte al menos dos predictores observados con la consulta, y la consulta tiene al menos dos predictores observados que varían en entrenamiento. Los predictores internos rellenos no cuentan como observados. Estas reglas etiquetan la estimación; no alteran el KNN estándar ni descartan vecinos antes de su cálculo. La suficiencia estadística no está garantizada por estos umbrales.

Iterative requiere al menos tres valores objetivo de entrenamiento, dos predictores observados que varían en entrenamiento y convergencia externa. Sus modelos pueden usar otros predictores completados internamente. La etiqueta no asegura que cada valor observado tenga peso importante en la regresión ni prueba información causal. La regularización puede producir resultados muy próximos a la media incluso en filas etiquetadas con soporte. Una consulta sin ningún predictor informativo se etiqueta training_only, aunque sklearn devuelva un número mediante sus cadenas de regresión. Si el objetivo es constante en entrenamiento también se etiqueta training_only.

## Sustituciones automáticas y pocos datos

- SimpleImputer ignora NaN y calcula mediana; con un solo valor del objetivo utiliza ese valor, siempre etiquetado training_only y con n registrado.
- KNN puede utilizar menos de tres vecinos si no hay suficientes donantes comparables. Un único predictor compartido permite una distancia válida de sklearn, pero no satisface nuestra regla de soporte. Se registran el número utilizado, los identificadores y el mínimo de predictores compartidos.
- Si ninguna distancia a donantes del objetivo está definida, KNN usa automáticamente la media del objetivo en entrenamiento. Con estandarización afín y pesos uniformes, desescalar devuelve la media en R. No es una mediana ni una estimación por vecinos. Se marca `automatic_training_mean_no_defined_distance`.
- Si la selección de vecinos no tiene observaciones informativas variables, se etiqueta training_only. Una distancia numérica no equivale por sí sola a información individual útil.
- Iterative comienza completando internamente con medianas. Después ajusta regresiones por variable y reutiliza esa secuencia para nuevas consultas. `skip_complete=False` permite estimar una ausencia en consulta aunque la variable estuviera completa en entrenamiento. Una sola columna activa puede terminar sin regresiones y quedarse con la inicialización; se marca explícitamente.
- Una fila enteramente vacía no identifica su perfil metabólico. KNN devuelve medias; Iterative puede devolver un perfil derivado de la inicialización y del modelo general. Ninguno se considera supported en esa situación.
- Una columna vacía en entrenamiento produce abstención en los tres métodos. No se usan datos de consulta para aprender esa columna ni se rellena cero como si estuviera medido.
- No hay sustitución entre métodos: una falta de convergencia no activa mediana ni KNN silenciosamente. Errores inesperados detienen la ejecución.

## Convergencia

SimpleImputer y KNNImputer no tienen bucle de convergencia. Sus campos `iterations` y `converged` quedan vacíos, no se declara éxito iterativo.

Iterative detiene el ajuste por la regla de sklearn sobre el cambio máximo de la matriz completada, dividido por el máximo valor absoluto observado, menor que tol; aquí se evalúa en la escala estandarizada. Con sample_posterior=False la regla está activa. Se conserva `n_iter_` y cada advertencia de ajuste. Si sklearn emite ConvergenceWarning, se marca converged=False y no se acepta el candidato como supported. El número de rondas por sí solo no basta para determinar convergencia. Si no hubo rondas, converged queda vacío.

`transform` aplica la secuencia de regresiones aprendida, no realiza un nuevo ajuste ni un criterio de parada independiente para cada consulta. `converged` describe el ajuste del entrenamiento, no garantiza convergencia individual de la consulta ni exactitud.

Además, el registro informa el máximo número de iteraciones de BayesianRidge y cuántos ajustes internos alcanzaron su límite. No se confunde convergencia del regresor con convergencia del imputador externo.

## Resultados y diagnóstico

En cada carpeta resultados/<dataset>/, validacion/ contiene las máscaras y las comparaciones con mediciones conocidas; imputacion/ contiene las estimaciones de los faltantes originales y estadísticas. ajustes.json registra entrenamiento, escaladores, iteraciones y advertencias mediante fit_id. manifest.json registra parámetros, versiones y hashes.

available conserva todos los candidatos finitos, incluyendo sustituciones generales y ajustes no convergentes. supported_plus_baseline exige supported a KNN e Iterative, pero mantiene explícitamente la mediana como referencia general. No convertir esa excepción en evidencia de soporte individual. Las estadísticas se entregan sobre celdas propias y sobre la intersección común de los tres métodos.

Los métodos y parámetros permanecen iguales al flujo anterior. Las 1.305 predicciones de validación Ishii conservan sus valores y etiquetas de soporte tras la simplificación. Se verificaron las dos ejecuciones y las 15 pruebas del motor y las máscaras. Esto demuestra consistencia de implementación, no certeza científica de los valores imputados.

## Documentación oficial consultada

- https://scikit-learn.org/1.8/modules/generated/sklearn.impute.SimpleImputer.html
- https://scikit-learn.org/1.8/modules/generated/sklearn.impute.KNNImputer.html
- https://scikit-learn.org/1.8/modules/generated/sklearn.impute.IterativeImputer.html

IterativeImputer sigue siendo experimental y requiere habilitar explícitamente sklearn.experimental.enable_iterative_imputer. La versión queda fijada para evitar cambios silenciosos de comportamiento.
