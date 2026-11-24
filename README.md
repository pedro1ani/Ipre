# IPre · Imputación de concentraciones de metabolitos en *Escherichia coli*

Ishii midió concentraciones intracelulares en 30 cultivos. La tabla de análisis contiene **135 metabolitos, 3.402 concentraciones medidas y 648 casillas sin medición**. El proyecto estudia dos preguntas: ¿qué error comete un método al reconstruir una concentración conocida? Y, cuando se completan los huecos reales, ¿qué cambia en la dispersión y en las relaciones entre metabolitos? Las concentraciones se conservan en **mM**.

| Dónde empezar | Qué explica |
| --- | --- |
| [01 · Auditoría, preparación y análisis exploratorio](01_auditoria_preparacion.ipynb) | Qué midió cada estudio, cómo se eligieron las señales, qué falta y cómo se distribuyen las mediciones originales. |
| [02 · Evaluación, imputación y análisis estadístico](02_imputacion_evaluacion.ipynb) | Cómo se ocultan mediciones, se comparan métodos, se evalúan intervalos y se estudia el efecto de completar las casillas vacías. |

Para consultar una sustancia en los cuadernos, cambie `METABOLITO`, por ejemplo a `fum` (fumarato), `cit` (citrato) o `atp`. El [inventario de señales](protocolo/resultados/inventario.csv) permite buscar su nombre y protocolo analítico.

## Qué se compara

Se oculta una concentración conocida de Ishii y se entrena sin los cultivos de su grupo de perturbación. También se ocultan varias concentraciones a la vez y se apartan series analíticas completas. Todos los métodos reciben las mismas máscaras. El [resumen de la regla original](protocolo/resultados/metricas.csv) muestra por qué hace falta evaluarla: en 3.302 mediciones comunes, su error medio absoluto (MAE) es 0,204 mM, frente a 0,163 mM de KNN3, 0,168 mM del iterativo y 0,246 mM de la mediana.

El cuaderno 02 compara **ocho estrategias**: cuatro métodos fijos (mediana, KNN3, KNN5 e iterativo), la regla original, dos reglas que cambian el respaldo y una regla jerárquica. Muestra el error en las mismas casillas, la cobertura sobre todos los intentos, la contribución de cada decisión y el intervalo de la diferencia entre métodos. El referente ideal por casilla conoce la medición: solo sirve para cuantificar el margen de mejora. Nunca se aplica a un dato faltante.

## Cómo interpretar una estimación

- **KNN** utiliza cultivos parecidos. Requiere al menos dos predictores compartidos y tres donantes que hayan medido el objetivo. Los vecinos se eligen de nuevo dentro de cada entrenamiento.
- **El método iterativo** estima mediante relaciones entre metabolitos y puede producir un intervalo por casilla. Su calibración se mide sobre concentraciones ocultas. La SD de los sorteos de una predicción y la SD experimental entre cultivos responden a preguntas distintas.
- **La mediana de respaldo** representa el centro del entrenamiento. Si se necesita información específica del cultivo y esta falta, se puede dejar la casilla pendiente. Las salidas distinguen método fijo, elección con evidencia, respaldo y ausencia de estimación.
- **Un blanco de Ishii** puede significar «no detectado» o «por debajo del límite de detección». El error medido al ocultar valores conocidos no demuestra la exactitud de los realmente ausentes. Tampoco basta con conservar la distribución para dar por correcta una predicción individual.

KNN es una referencia sencilla para reconstruir concentraciones individuales con soporte. La regla adaptativa debe demostrar una mejora estable frente a ese referente en ambas pruebas externas. Para estudiar varianzas, balances y relaciones con proteínas o flujos conviene contrastar la matriz observada con varias matrices imputadas, manteniendo identificadas las estimaciones y sus respaldos.

## Fuentes y ejecución

La [base permanente en mM](data/concentraciones_mM.db), los [archivos experimentales](sources/) y la [auditoría](auditoria/) conservan las mediciones y su procedencia. De las 135 señales, 17 tienen equivalencias que permiten probar las ocho condiciones de Gerosa (MMC2); las otras 118 se estudian dentro de Ishii. Gerosa aporta entrenamiento auxiliar, nunca la medición que verifica una estimación. El [método](protocolo/METODO.md) documenta los criterios y los controles.

Desde la raíz del repositorio, con **Python 3.12.14** en un entorno virtual:

```bash
python -m pip install -r requirements.txt
python -m auditoria.verificar
python -m protocolo.revision --procesos 4
python -m unittest discover -s tests
```

En Jupyter seleccione el mismo entorno y ejecute los cuadernos de arriba abajo. El cuaderno 01 solo necesita la preparación de fuentes. El 02 reutiliza los resultados si están completos; en su primera ejecución puede tardar varios minutos porque ajusta modelos y genera sorteos posteriores.

Los [resúmenes de la validación de referencia](protocolo/resultados/) están en el repositorio. Las tablas detalladas, matrices, sorteos y figuras de la comparación de estrategias se generan en su equipo dentro de `salida_protocolo/`, separadas de SQLite y de los archivos versionados. El [índice de salidas](protocolo/resultados/README.md) indica qué abrir para responder cada pregunta. Para cambiar la carpeta de trabajo, defina la variable de entorno `IPRE_DETALLE`.

