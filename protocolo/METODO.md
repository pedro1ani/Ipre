# Cómo se estima una concentración que falta

En un cultivo de *Escherichia coli* pueden estar medidas las concentraciones de varias sustancias y faltar la de otra, por ejemplo la de fumarato. La sustancia cuya concentración se quiere estimar se llama **objetivo**. Las otras sustancias que se usan en el cálculo se llaman **predictores**. El protocolo responde dos preguntas por separado: ¿qué tan bien se recupera una medición cuando se la oculta al método?, y ¿qué evidencia hay para aplicar ese método cuando no se conoce el valor?

## Qué representa cada dato

La matriz es una tabla con una fila por cultivo y una columna por señal analítica. Esta disposición mantiene juntas las concentraciones medidas en el mismo cultivo. El [inventario](resultados/inventario.csv) identifica cada señal por su nombre en Ishii, fila del archivo original, protocolo de medición, compartimento y unidad. La etiqueta del laboratorio permite seguir una señal dentro del estudio; no siempre basta para asignarle un único compuesto o para identificarla con una señal de otro estudio.

Los 30 cultivos principales crecen con glucosa a la misma tasa de dilución (0,2 por hora). Se distribuyen en 25 grupos según la perturbación experimental: por ejemplo, dos cultivos en los que se inactiva el mismo gen se mantienen en el mismo grupo. Los cuatro cultivos de la cepa sin modificar que crecen a otras tasas de dilución se usan solo para comprobar el efecto de incluirlos; también se respeta su grupo. Medir varias veces una proteína en el mismo cultivo cuenta como una observación biológica: se usa la mediana de esas mediciones técnicas.

Las concentraciones de Ishii se leen en **mM** (milimoles por litro de volumen intracelular). Gerosa publica micromoles por gramo de peso seco. Para obtener concentraciones comparables se divide la cantidad por gramo por el volumen intracelular documentado de **2,3 mL por gramo de peso seco**. La unidad resultante, µmol/mL, equivale a mM. Se conserva la cantidad original y su procedencia.

La [base de trabajo en mM](../data/concentraciones_mM.db) conserva 19.950 posiciones de Ishii, 344 de MMC2 y 189 de MMC3. La [base experimental](../sources/ecoli_auditada_v3.db) reúne además proteínas, corridas técnicas y variables fisiológicas. [La auditoría](../auditoria/README.md) describe qué tablas se consultan, cómo se cotejan por identificador y de dónde se leen las matrices. Las ocho condiciones estacionarias de MMC2 entran solo en los ensayos de aporte de datos externos. MMC3 documenta transiciones y MMC4 aporta las anotaciones de reacciones.

## Qué señales tienen datos suficientes

De las 579 señales registradas, 149 tienen al menos diez mediciones. Entre esas 149 hay 11 con identidad incierta y tres que reúnen distintas sustancias en una misma señal. Las 135 restantes se consideran posibles objetivos. Ishii y Gerosa comparten una identidad documentada para 17 de ellas. Las otras 118 se pueden estimar dentro de Ishii sin suponer una equivalencia entre estudios.

En cada prueba se apartan algunos cultivos. En los que quedan para el **entrenamiento** (los que puede usar el método), el objetivo necesita al menos diez mediciones en cinco grupos y dos valores distintos. Para relacionarlo con otras sustancias se exigen además dos predictores con al menos cinco cultivos en los que ambos estén medidos: el predictor y el objetivo. Con estos requisitos hay entre 127 y 132 objetivos aptos en los cinco entrenamientos principales. La [tabla de umbrales](resultados/umbrales.csv) también muestra cuántos quedarían con mínimos de ocho o doce.

Por ejemplo, DHAP, NADH y NADPH solo tienen 8, 6 y 6 mediciones en los 30 cultivos: quedan fuera de la prueba principal. Otras señales tienen diez en total pero pueden quedar por debajo del mínimo cuando se aparta un grupo.

## Cómo se comprueba una predicción sin mostrarle la respuesta al método

Una **máscara** indica qué casillas medidas se ocultan durante la prueba. La concentración original se conserva para calcular el error al final. Todas las alternativas reciben las mismas casillas ocultas. Las únicas respuestas conocidas con las que se las evalúa son mediciones originales de Ishii.

| Prueba | Qué permite comprobar |
| --- | --- |
| Una medición oculta por vez, con grupos apartados | Se distribuyen los 25 grupos de perturbación en cinco bloques. Cada bloque sirve una vez como prueba y el método se ajusta con los grupos restantes. Así se evalúa un cultivo cuya perturbación no aportó datos al ajuste. |
| Varias mediciones ocultas a la vez | Se buscan combinaciones que realmente faltan juntas en los cultivos de entrenamiento y se ocultan donde las concentraciones sí son conocidas. Permite ver qué ocurre cuando también falta información que podría ayudar a la estimación. |
| Una serie analítica apartada | Se retiran todos los cultivos de una misma serie de medición y, del entrenamiento, cualquier otro cultivo que comparta su grupo de perturbación. Así se comprueba cuánto depende el resultado de la serie. |

El **perfil de ausencia** describe cuáles de cuatro sustancias están sin medir en el cultivo: fumarato, fosfoenolpiruvato (PEP), sedoheptulosa-7-fosfato (S7P) y succinato. Si una de las cuatro es el objetivo, se consideran las otras tres. Estas referencias sirven para describir qué falta; los predictores del modelo pueden ser otras sustancias. Ningún patrón conjunto de las 135 columnas se puede reproducir por completo con mediciones conocidas en los cultivos de prueba. Al considerar el objetivo y las cuatro referencias sí hay [2.201 casos de prueba conjunta por perturbación](resultados/mascaras_resumen.csv). Una medición puede participar en más de un caso, con diferentes acompañantes ocultos.

## Cómo calcula cada método

| Método | De dónde sale el valor |
| --- | --- |
| Mediana | Se ordenan las concentraciones conocidas del objetivo en el entrenamiento y se toma su valor central. No necesita conocer los otros metabolitos del cultivo. Se implementa con `SimpleImputer`. |
| Vecinos (KNN) | Se buscan cultivos con el objetivo medido y al menos dos predictores medidos en común con el cultivo que se quiere estimar. Se comparan las concentraciones disponibles y se promedian las del objetivo en los cultivos más parecidos. Se ensayan tres o cinco vecinos con el mismo peso y se requieren al menos tres. Si hay menos de cinco se registra cuántos se usan. Se implementa con `KNNImputer` y distancia euclídea adaptada a datos faltantes. |
| Iterativo | Se comienza con medianas y se estima cada variable a partir de las demás. El cálculo se repite para refinar las estimaciones. Se usa `IterativeImputer` con un modelo de regresión lineal (`BayesianRidge`), al menos dos predictores medidos en el cultivo y un límite inferior de cero mM. El ajuste determinista permite hasta cinco rondas, con tolerancia de 0,001; el [registro de convergencia](resultados/convergencia.csv) indica si se alcanza el máximo sin satisfacer el criterio de parada. |

KNN y el método iterativo usan hasta tres metabolitos predictores por objetivo. Se favorecen las sustancias que tienden a variar junto con el objetivo y están medidas con él en suficientes cultivos. Para ordenar las candidatas se multiplica la correlación absoluta de rangos de Spearman por la fracción de mediciones del objetivo que tienen pareja. Esta relación ayuda a predecir, aunque por sí sola no demuestra una conexión bioquímica.

Para que una sustancia no pese más solo por sus valores numéricos, cada variable se centra y se divide por su desviación estándar con `StandardScaler`. Tanto el escalado como la selección de predictores utilizan solo el entrenamiento. Tras el cálculo se recupera la escala en mM, que también se usa para expresar el error.

## Cuándo se elige un método y cuándo se usa información adicional

Dentro de cada entrenamiento se forman **tres bloques internos** para probar las opciones. Se compara el *error absoluto medio* (MAE): cuánto se aleja la predicción de la medición, en promedio y sin considerar el signo. Por ejemplo, estimar 0,25 mM cuando la medición es 0,30 mM da un error absoluto de 0,05 mM. Para elegir por objetivo y perfil se requieren al menos diez cultivos comparables, cinco grupos y dos bloques internos. La ventaja en MAE debe repetirse en al menos dos bloques. La elección se fija antes de evaluar el bloque externo, que no participa en esa decisión.

El estudio de Gerosa aporta ocho condiciones de cultivo en el archivo **MMC2**. Para medir su aporte se comparan dos entrenamientos que usan los mismos metabolitos compartidos: **A** usa solo Ishii y **B** usa Ishii y MMC2. Los cultivos de prueba son los mismos. Por ejemplo, con el método iterativo, el [MAE para el fumarato](resultados/MMC2_pareado.csv) es 0,019 mM con Ishii y 0,098 mM al añadir MMC2; para PEP los valores son 0,067 y 0,058 mM, respectivamente. El aporte depende de la sustancia y del método.

También se ensaya, una por vez, la captación de carbono (cantidad de carbono consumido por gramo de peso seco y hora) o la abundancia de una proteína. Para proteínas se estudian candidatas relacionadas con la reacción y candidatas elegidas por su disponibilidad en el entrenamiento. El control se ajusta sobre los mismos cultivos sin esa variable. Usar una variable adicional requiere que mejore frente a los controles correspondientes y que esté medida en el cultivo donde se quiere aplicar la regla.

| Estado en los archivos | Interpretación para una casilla sin medición |
| --- | --- |
| `validada` | La comparación interna respalda la **elección del método** para ese objetivo y perfil. La concentración real de esa casilla sigue siendo desconocida. La propia mediana puede ser el método elegido. |
| `respaldo_general` | Los datos del objetivo permiten calcular una mediana de entrenamiento, pero las pruebas no respaldan elegir un método específico para ese caso. |
| `pendiente` | La regla no puede producir una estimación con los datos disponibles. Si la opción elegida necesita una tasa o proteína ausente, solo se usa una alternativa sin ella si ya tenía respaldo para ese perfil. |

## Cómo interpretar el error y el intervalo

El [informe de errores](resultados/metricas.csv) distingue cuántas veces responde un método de cuán cerca está de la medición. Una comparación usa todos los valores que estima cada método; la otra usa **solo las mismas mediciones** que todos pueden estimar.

| Indicador | Qué pregunta responde |
| --- | --- |
| MAE en mM | ¿A qué distancia de la medición queda la estimación, en promedio? El promedio *micro* da el mismo peso a cada medición. El *macro* calcula primero el error de cada metabolito y luego da el mismo peso a cada metabolito. |
| RMSE en mM | ¿Hay errores grandes que el promedio absoluto disimula? Se elevan las diferencias al cuadrado antes de promediarlas y sacar la raíz. |
| Sesgo en mM | ¿El método tiende a estimar por encima o por debajo? Se promedia «estimación menos medición»; un valor positivo indica sobreestimación en promedio. |
| Error relativo a la variabilidad | ¿Es grande la diferencia frente a la variación habitual del metabolito? Se divide el error absoluto por el rango intercuartílico del entrenamiento Ishii: la distancia entre los percentiles 25 y 75. Se omite si ese rango es cero. |

En la regla original, solo las estimaciones del método iterativo elegido llevan **intervalo de predicción**. El modelo genera 200 concentraciones plausibles por casilla, o 100 cuando utiliza una proteína, mediante muestreo de su distribución posterior. Los percentiles 2,5 y 97,5 de esos valores definen los extremos del intervalo nominal del 95 %. Son límites del modelo, no límites del instrumento de medición.

Hay dos recuentos distintos. En la prueba individual por perturbación, de **3.402 mediciones ocultas**, la regla entrega 3.325 estimaciones y **809 intervalos**. De esos 809, **704 contienen la medición**: la cobertura observada es 87,0 %. Al apartar una serie completa contienen la medición 488 de 640 intervalos (76,25 %). Para una casilla sin medición no se puede saber si su intervalo contiene el valor real.

## Cómo seguir una estimación hasta su origen

Los archivos por cultivo se generan en `salida_protocolo/`. En `auditoria_huecos.csv` se puede buscar un cultivo y un metabolito para ver la celda del Excel, el motivo del blanco, el método, las fuentes y la estimación en mM. Para KNN también se registran los cultivos vecinos. Si hay casos comparables en la prueba externa, se informa su error y cuántas mediciones y grupos lo respaldan. **Ese error describe otros casos con respuesta conocida**, no el error individual de la casilla vacía.

La [tabla de estadísticas](resultados/descriptivos.csv) describe los mismos 30 cultivos: cuenta las mediciones, las estimaciones y las casillas pendientes, y resume las concentraciones mediante media, mediana, desviación estándar muestral y rango intercuartílico. El coeficiente de variación (CV) expresa la desviación estándar como porcentaje de la media; se deja sin valor si hay menos de dos datos o la media no es positiva. Se muestran tanto los resultados disponibles de cada método como los de las mismas casillas estimables por todos. El cambio en una media al añadir estimaciones describe el contenido de la tabla, no valida las estimaciones.

| Si se quiere revisar… | Archivo del programa |
| --- | --- |
| De qué hoja sale una concentración y cómo se convierte a mM | [datos.py](datos.py), [identificadores de señales](fuentes.py) y [verificación de fuentes](../auditoria/verificar.py) |
| Cómo se cotejan tasas, corridas de proteínas y candidatos de control | [fisiologia.py](fisiologia.py) y [proteinas.py](proteinas.py) |
| Qué cultivos se apartan, qué valores se ocultan y cómo se eligen predictores | [diseno.py](diseno.py) y [validar.py](validar.py) |
| Cómo funcionan los tres métodos y cómo se comparan las fuentes | [modelos.py](modelos.py) y [ensayos.py](ensayos.py) |
| Qué evidencia permite elegir un método y cómo se aplica a cada casilla | [regla.py](regla.py) y [aplicar.py](aplicar.py) |
| Cómo se obtienen los intervalos y las tablas de resultados | [intervalos.py](intervalos.py) e [informe.py](informe.py) |

Para reproducir la evaluación se fija la semilla 20261109; los sorteos de los intervalos usan semillas consecutivas desde 20262109. `reglas.json` guarda las decisiones y `congelado.json` guarda huellas digitales (SHA-256) que permiten comprobar que las entradas y las máscaras corresponden a la evaluación. Las [instrucciones de ejecución](../README.md) indican en qué orden usar los programas.

## Ocho estrategias, una comparación

**Métodos fijos.** Mediana, KNN3, KNN5 e iterativo usan el panel amplio de Ishii y el mismo diseño estadístico. Sus valores de entrenamiento y sus predictores se fijan sin ver la prueba. Si no se puede aplicar un método, se conserva la ausencia. La regla original puede usar además MMC2 o una covariable si supera sus controles en los mismos cultivos. El respaldo KNN conserva cada elección validada de la regla original; en los otros casos intenta KNN3 y, si no tiene soporte, la mediana. El respaldo interno funciona igual, pero aprende el método global en las pruebas internas. Una mediana elegida por evidencia local y una mediana que solo sustituye una estimación fallida llevan etiquetas distintas.

**Jerarquía.** El punto de partida es el método con menor MAE en mediciones individuales comunes a los cuatro métodos dentro del entrenamiento. Se intenta mejorar la elección con evidencia para metabolitos parecidos, después para el propio metabolito y finalmente para su perfil. Así tiene prioridad la información más específica que demuestre utilidad. La similitud se define con dos variables del entrenamiento: fracción ausente y CV por encima/debajo de sus respectivas medianas. No supone una conexión química. Se exigen al menos tres metabolitos distintos del objetivo antes de prestar evidencia. El nivel registrado es el que **cambió** la decisión. La evidencia global o prestada se marca como respaldo; no se confunde con validación local.

Para aceptar un cambio se comparan, en las mismas casillas, el mejor candidato y la elección del nivel anterior: ganancia de MAE superior al margen estipulado y ventaja en al menos dos bloques. Se ensayan mínimos de 8, 10 y 12 casillas distintas y un margen de cero o un error estándar agrupado por perturbación. Se mantienen cinco grupos y dos bloques. Para el nivel objetivo/perfil cada cultivo aporta una vez; para niveles compartidos se usan ocultamientos individuales, por lo que cada casilla cultivo–objetivo aporta una vez. El nivel global también requiere esos mínimos. Si faltan pruebas, se conserva el nivel anterior. Una mediana general se etiqueta como tal y, si no se puede calcular, se deja pendiente.

**Cómo se decide la complejidad.** Cada entrenamiento externo tiene tres bloques internos. Para ensayar las reglas sin entrenarlas con el bloque que las puntúa, cada entrenamiento interno se vuelve a separar en tres por perturbación. En conjuntos pequeños puede no quedar ningún metabolito elegible: no se inventa evidencia. En las casillas individuales comunes de la validación interna se elige la primera configuración, ordenada de menor a mayor complejidad, cuya pérdida frente a la mejor no exceda un error estándar pareado por grupo. Se vuelve a aprender la elección con todo el entrenamiento externo y solo entonces se compara con la verdad externa. En la aplicación final se repite el procedimiento usando los 30 cultivos. Los resultados externos de una variante más flexible sirven para estudiarla, no para instalarla retrospectivamente.

## Qué significan las comparaciones y figuras

- **Error individual:** MAE micro y macro, RMSE, sesgo (estimación menos medición), mediana y percentiles 90, 95 y 99 del error absoluto. Cada escala de IQR proviene del entrenamiento Ishii de su bloque. Se informa el número de denominadores utilizables. La SD del MAE entre grupos y entre bloques describe estabilidad, no la precisión de un valor individual.
- **Casillas:** «disponibles» incluye todos los intentos y expresa cobertura; «comunes» restringe el error a las consultas que responden las ocho estrategias. Para un contraste entre dos estrategias se usa su intersección y se declara el denominador, que puede diferir del común a ocho. Los IC del 95 % corresponden a 2.000 remuestreos de perturbaciones completas, conservando juntos cultivos y metabolitos del mismo grupo. Describen la variación de la diferencia usando las predicciones externas calculadas; no reajustan los modelos en cada réplica. El límite ideal selecciona retrospectivamente la mejor de las cuatro predicciones fijas disponibles y no integra ninguna regla operativa.
- **De dónde viene la pérdida:** para un subconjunto se suma `|regla − observado| − |KNN3 − observado|` y se divide por el número total de casillas comunes. Los aportes por estado se suman; no se deben sumar dimensiones que describen las mismas casillas, como método y perfil. El archivo de alternativas permite contrastar además cada opción disponible, incluida la información auxiliar, sobre la misma medición.
- **Distribución:** la SD y la varianza son muestrales; el CV es 100 × SD/media. La asimetría requiere ocho valores y dispersión. Los porcentajes de cambio se omiten si el denominador es casi cero (máximo de 10⁻¹² y 10⁻⁶ veces la escala de la magnitud). Se presentan las casillas disponibles y las mismas casillas completables por las ocho estrategias; las mediciones originales son idénticas en ambas. Los posibles extremos se rastrean a su hoja, celda y protocolo usando las cercas de Tukey; se compara su influencia en el error sin borrarlos.
- **Correlaciones y PCA:** Pearson y Spearman usan al menos ocho pares y guardan el número de cultivos antes y después. El cambio relativo se omite si |r original| < 0,1. En los huecos reales pueden incorporarse nuevos pares; en el ocultamiento conocido se calcula la dispersión de las predicciones externas y las correlaciones de todos los métodos sobre los mismos cultivos y pares con coeficiente definido. Son reconstrucciones individuales, no un borrado simultáneo de los 135 metabolitos. El PCA solo visualiza las cuatro anclas: un escalador y dos componentes ajustados a cultivos originalmente completos, con proyección de los cultivos que se completan después. No imputa datos.
- **Incertidumbre:** además del método iterativo elegido por la regla, se generan sorteos del iterativo fijo para todas las mediciones individuales aplicables por perturbación y para los huecos reales estimables. Se conserva la SD posterior por casilla, los percentiles centrales del 50, 80, 90 y 95 %, el ancho frente al error y la calibración empírica. El ancho depende en parte de la escala del metabolito: se informa también la relación ajustada al IQR y por objetivo en la elección original. Para un hueco real no hay verdad con la cual calibrar.
- **Uso biológico:** la covariación y el error de reconstrucción son evidencias distintas. Una concentración no mide un flujo ni prueba una interacción enzimática. Un balance o una correlación que cambia al sustituir los respaldos por valores ausentes es sensible a la imputación y requiere evidencia independiente.

El [índice de resultados](resultados/README.md) relaciona cada pregunta con su tabla. [estrategias.py](estrategias.py) aprende y aplica las decisiones; [comparacion.py](comparacion.py) comprueba las predicciones contra verdades conocidas; [distribuciones.py](distribuciones.py) calcula estadísticos y recuperación; [figuras.py](figuras.py) crea las figuras y permite consultar cualquier metabolito. El punto de entrada es [revision.py](revision.py).

