# Tres métodos de imputación y diagnóstico de soporte

Implementación: scikit-learn 1.8.0. Se usan exclusivamente SimpleImputer, KNNImputer e IterativeImputer. No hay un estimador KNN propio ni PCA en el código activo. El cálculo de distancias en el diagnóstico de KNN solo identifica los donantes usados por sklearn; no sustituye su estimación.

## Alcance y ejecución

El flujo principal es `python ejecutar_estudio.py --output resultados_nuevos`: valida los tres métodos en Ishii y en Ishii + MMC2 con máscaras comunes, y después estima los 19 huecos objetivo de Ishii por método y escenario (114 candidatos). La imputación real excluye el grupo de perturbación de la consulta: KO05/KO05x juntos y todos los controles WT juntos. Las siete variables pueden completarse internamente, pero solo se exportan imputaciones de los cuatro objetivos; los otros tres predictores conservan sus ausencias originales.

Instalación: `python -m pip install -r requirements_imputacion.txt`. Diseño y resultados en `VALIDACION_Y_RESULTADOS.md`. Configuración del estudio: `config_validacion.json`, con parámetros de métodos en `config_imputacion_ishii.json` para ambos escenarios.

El ejecutor auxiliar `ejecutar_imputacion.py --dataset ishii|mmc2` permanece disponible para aplicaciones independientes. No ofrece opción combinada ni contiene una prohibición interna por origen. La comparación combinada se implementa en el nuevo ejecutor de estudio. No se han imputado muestras MMC2; su aplicación futura requiere validación independiente. En particular, no evaluar Gerosa glucosa como un valor ocultado mientras el denominador revela su concentración.

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

En el estudio completo hay 37 entrenamientos distintos por escenario para cada método, reutilizados cuando las exclusiones coinciden. Iterative converge en 31/37 ajustes Ishii y 4/37 combinados. En el subconjunto de siete ajustes usado para huecos reales converge en 7/7 Ishii y 1/7 combinados. La consulta del único ajuste combinado convergente carece de predictores observados; por tanto no aporta imputaciones reales con soporte.

## Salidas y políticas

Las salidas principales están en `resultados_estudio/`. `ajustes_y_particiones.json` permite vincular cada candidato a su entrenamiento mediante fit_id y contiene IDs, escaladores, iteraciones y advertencias. Los CSV de validación e imputación conservan candidatos, soporte y mecanismo.

`available` conserva todos los candidatos finitos, incluidas sustituciones generales y ajustes sin convergencia. `supported_plus_baseline` acepta solo supported para KNN e Iterative y mantiene explícitamente la mediana como referencia general. Esta excepción no convierte la mediana en una estimación individual con soporte. Las matrices incluyen esa política en el nombre; las bases y matrices fuente permanecen intactas.

Las estadísticas distinguen celdas propias, comunes a tres métodos dentro de cada escenario y comunes a los seis resultados. Una intersección vacía no permite evaluar cambios comparables por imputación.

## Verificación de implementación

Quince pruebas automatizadas cubren comportamientos del motor, entradas, máscaras, exclusión completa de grupos/series y ausencia de influencia de valores ocultados sobre la predicción. El verificador del estudio comprueba además el entrenamiento del escalador, las verdades originales y la preservación de las mediciones. Dos ejecuciones completas dieron los mismos 24 archivos CSV/JSON. Estas comprobaciones verifican implementación y reproducibilidad; no garantizan precisión científica en los huecos reales.

## Documentación oficial consultada

- https://scikit-learn.org/1.8/modules/generated/sklearn.impute.SimpleImputer.html
- https://scikit-learn.org/1.8/modules/generated/sklearn.impute.KNNImputer.html
- https://scikit-learn.org/1.8/modules/generated/sklearn.impute.IterativeImputer.html

IterativeImputer sigue siendo experimental y requiere habilitar explícitamente sklearn.experimental.enable_iterative_imputer. La versión queda fijada para evitar cambios silenciosos de comportamiento.
