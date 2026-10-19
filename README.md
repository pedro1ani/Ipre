# IPre — Imputación de una base por ejecución

## Los dos códigos principales

**ejecutar_imputacion.py** prepara los datos y contiene las dos secciones del análisis. **imputadores.py** implementa la interfaz de los tres métodos de scikit-learn. El ejecutor usa una sola matriz de entrada, elegida con --dataset; no combina Ishii y MMC2.

### Sección 1: evaluar con mediciones conocidas

Oculta valores originales, entrena sin la muestra ni su grupo y estima los valores ocultados. Usa exactamente las mismas máscaras para mediana, KNN e Iterative. Compara estimación con medición mediante MAE en mM y ratios, cobertura y celdas comunes. Guarda resultados en validacion/ antes de pasar a la segunda sección.

### Sección 2: imputar candidatos reales

Identifica como candidatos los faltantes de los cuatro objetivos (fumarato, PEP, s7p y succinato). Aplica cada método por separado, conservando todas las mediciones originales. Guarda candidatos, diagnósticos, matrices y estadísticas antes/después en imputacion/. Los otros tres metabolitos son predictores: sus ausencias no se rellenan en las matrices finales.

## Ejecutar

```bash
python -m pip install -r requirements_imputacion.txt
python ejecutar_imputacion.py --dataset ishii --output resultados_ishii_nuevos
```

Para analizar MMC2 de forma independiente:

```bash
python ejecutar_imputacion.py --dataset mmc2 --output resultados_mmc2_nuevos
```

La carpeta de salida debe ser nueva. Sin argumentos se selecciona Ishii y resultados/ishii. Se entregan ejecuciones independientes verificadas en resultados/ishii y resultados/mmc2. No es necesario ejecutar otro código para obtener evaluación e imputación; se retiró ejecutar_estudio.py.

Cada salida tiene solo dos carpetas de análisis: validacion/ e imputacion/, más ajustes.json y manifest.json para trazabilidad. Los gráficos opcionales se guardan en figuras/.

## Ver y comprobar

```bash
python verificar_estudio.py --results resultados/ishii
python -m pip install -r requirements_descriptivos.txt
python graficar_estudio.py --results resultados/ishii
python -m unittest discover -s tests -v
```

Los gráficos comparan directamente valor medido ocultado (eje horizontal) con valor imputado (eje vertical). La diagonal representa acierto exacto. Son gráficos de validación, no de huecos reales cuya verdad es desconocida.

## Archivos de apoyo

| Archivo | Función |
|---|---|
| normalizar.py | Construye las matrices en mM y ratios a partir de los datos preparados. |
| convertir_unidades.py | Reconstruye la conversión desde fuentes; no hace falta repetirla para imputar. |
| analizar_descriptivos.py | Describe únicamente mediciones originales. |
| graficar_estudio.py | Grafica resultados existentes; no entrena. |
| verificar_estudio.py | Revisa cada ejecución independiente. |
| verificar_proyecto.py | Revisa integridad y conversiones de la distribución. |
| config_imputacion_ishii.json / config_imputacion_mmc2.json | Parámetros y objetivos del estudio seleccionado. |

sources/ y data/ conservan las fuentes y conversiones. normalizacion/ contiene las matrices preparadas. docs/IMPUTACION.md describe el motor y docs/VALIDACION_Y_RESULTADOS.md el diseño y las limitaciones. Los registros técnicos permiten reproducir los resultados, pero no son entradas adicionales que el usuario deba seleccionar.

Se mantienen las 30 muestras Ishii: 29 con dilución 0,2 h⁻¹ documentada y KO05x con metadatos incompletos. Todos los escaladores y modelos se ajustan exclusivamente con entrenamiento; no se usan logaritmos. La referencia Gerosa glucosa sigue fija, con incertidumbre documentada; no se propaga su incertidumbre a los ratios.
