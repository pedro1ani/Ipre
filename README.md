# IPre: evaluar e imputar una base

## Por dónde empezar

1. **Explorador_Imputacion_Simplificado.ipynb** muestra datos, errores y candidatos.
2. **ejecutar_imputacion.py** ejecuta el análisis: prepara una base, evalúa mediciones ocultadas y completa huecos reales.
3. **imputadores.py** ajusta mediana, KNN e Iterative y registra su soporte.

Cada ejecución usa Ishii o MMC2 por separado. Los datos preparados y resultados de ambas bases ya están incluidos.

## Ejecutar

```bash
python -m pip install -r requirements.txt
python ejecutar_imputacion.py --dataset ishii --output resultados_ishii_nuevos
```

Para MMC2, cambia `--dataset mmc2` y usa otra carpeta de salida. El ejecutor exige una carpeta nueva para proteger resultados anteriores. Las dos únicas opciones de ejecución son base y carpeta de salida.

Para explorar lo entregado, abre el notebook, selecciona el Python instalado y pulsa **Ejecutar todo**. El notebook lee `resultados/ishii` o `resultados/mmc2`; puedes cambiar `RES` si quieres explorar una nueva ejecución. No exporta archivos.

## Dos partes del análisis

- **validacion/**: máscaras compartidas, valores medidos ocultados, estimaciones, MAE en mM/R, cobertura y casos comunes.
- **imputacion/**: estimaciones de los faltantes reales, matrices completadas y estadísticas antes/después.

`available` incluye todos los candidatos finitos. `supported_plus_baseline` exige soporte a KNN e Iterative y conserva la mediana como referencia general. Más cobertura no implica mejor precisión.

## Preparación y comprobación

| Archivo | Propósito |
|---|---|
| convertir_unidades.py | Verificar fuentes y convertir Gerosa a mM. |
| normalizar.py | Seleccionar muestras y construir matrices con referencia Gerosa glucosa. |
| analizar_descriptivos.py | Describir mediciones originales. |
| graficar_estudio.py | Graficar valores ocultados frente a estimados. |
| apoyo.py | Compartir escritura de archivos, huellas y estadísticos. |
| verificar_estudio.py | Comprobar particiones, verdades y conservación de observados. |
| verificar_proyecto.py | Comprobar integridad, unidades y fuentes. |

```bash
python verificar_proyecto.py
python verificar_estudio.py --results resultados/ishii
python verificar_estudio.py --results resultados/mmc2
python -m unittest discover -s tests -v
```

## Configuración

`config_imputacion.json` contiene seis parámetros comunes a ambas bases. `config_normalizacion.json` fija la selección científica y los objetivos; es la única definición de esos objetivos. Los parámetros se justifican en docs/IMPUTACION.md.

Las fuentes y datos originales permanecen en sources/ y data/. Los manifiestos registran versiones y huellas. docs/REFACTORIZACION.md explica los cambios y su verificación; docs/VALIDACION_Y_RESULTADOS.md describe las pruebas y sus limitaciones.
