# IPre: análisis ampliado a diez metabolitos

Abre **Explorador_Ampliado.ipynb** para ver el análisis nuevo. Selecciona `ishii`
o `mmc2` y ejecuta todas las celdas. Los resultados ya están incluidos.

## Reproducir

```bash
python -m pip install -r requirements.txt
python estudio_ampliado.py --dataset ishii --output resultados_recalculados/ishii
python estudio_ampliado.py --dataset mmc2 --output resultados_recalculados/mmc2
```

El directorio de salida debe ser nuevo. Para regenerar tablas y figuras del
análisis entregado:

```bash
python analizar_ampliacion.py --dataset ishii
python analizar_ampliacion.py --dataset mmc2
```

El segundo script lee `resultados_ampliados/`. Para inspeccionar otro resultado,
colócalo en esa ruta después de conservar la ejecución anterior. Las ejecuciones
posteriores de Iterative pueden tardar varios minutos. No hace falta repetirlas
para abrir el notebook.

## Archivos principales

- `estudio_ampliado.py`: prepara diez metabolitos, selecciona parámetros dentro
  del entrenamiento, evalúa grupos apartados y completa candidatos con incertidumbre.
- `analizar_ampliacion.py`: estadísticas, pruebas exploratorias, gráficos y conclusiones.
- `Explorador_Ampliado.ipynb`: revisión visual sin exportaciones ni entrenamiento.
- `docs/METODOLOGIA_AMPLIADA.md`: decisiones, supuestos y referencias.
- `resultados_ampliados/ishii/` y `resultados_ampliados/mmc2/`: análisis nuevos separados.

El ranking se basa en macro MAE_R externo. Consulta siempre errores por metabolito,
cobertura y anchura de intervalos. No hay garantía de precisión en huecos reales.

## Relación con la versión anterior

Esta entrega amplía la metodología; no es otra refactorización con resultados
idénticos. Las fuentes y el flujo de siete metabolitos se conservan para
trazabilidad y para las pruebas previas por patrones/series. La documentación
anterior describe ese flujo, no el análisis ampliado. `MANIFEST_PROYECTO.json`
y `manifest.json` son registros de la versión anterior; el manifiesto vigente
de esta entrega es `MANIFEST_AMPLIACION.json`.

Para actualizar tu repositorio puedes agregar los tres archivos principales,
`tests/test_ampliacion.py`, `verificar_ampliacion.py`, la documentación nueva,
`resultados_ampliados/` y el manifiesto nuevo. Mantén `.git` y tus trabajos personales.
Las bases originales no fueron modificadas. No necesitas sustituir tus fuentes.

## Verificación

```bash
python -m unittest discover -s tests -p test_ampliacion.py -v
python verificar_ampliacion.py
```
