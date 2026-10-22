# Refactorización integral — 22 de octubre de 2026

## Diagnóstico y propuesta aplicada

| Problema encontrado | Cambio aplicado | Razón |
|---|---|---|
| Líneas con varias operaciones y nombres de una letra | Reescritura con nombres descriptivos, pasos explícitos y formato de 79 columnas | Permitir seguir el cálculo sin descifrar expresiones compactas. |
| Dos configuraciones casi iguales y opciones que el motor rechazaba cambiar | Una configuración común de seis parámetros; opciones fijas junto a cada constructor sklearn | Reducir decisiones expuestas sin modificar el método efectivo. |
| Objetivos repetidos en configuraciones | Una definición en config_normalizacion.json | Evitar selecciones divergentes entre preparación e imputación. |
| Generación de máscaras con función interna que alteraba listas externas | Construcción explícita de ensayos y posterior extracción de verdades | Separar qué se oculta de lo que se usa para evaluar. |
| Predicción y diagnósticos en un bloque largo | Clase Imputador con ajuste, transformación y dos diagnósticos específicos | Conservar el estado entrenado y distinguir cálculo de soporte. |
| Escritura y descriptivos repetidos | apoyo.py con cuatro operaciones compartidas | Reutilizar solo las tareas realmente comunes. |
| Selección de estudios por prefijos en descriptivos | Lectura explícita de las matrices de cada estudio | Mantener identificable el origen de los datos. |
| Reglas de aceptación y estadísticas recalculadas en el notebook | Lectura de casos comunes y estadísticas del ejecutor | Evitar implementaciones paralelas que puedan discrepar. |
| Tres listas de dependencias solapadas | Un requirements.txt con versiones fijadas | Una instalación para ejecutar y explorar todo el proyecto. |
| Comprobaciones de conversión dependientes de assert | Condiciones explícitas con mensajes breves | Mantener controles incluso si Python se ejecuta con -O. |
| Conversión ligada a una carpeta fija y metadatos preparados manualmente | --output y exportación reproducible de metadatos Ishii | Permitir reconstruir en una carpeta nueva y verificar equivalencia. |

El flujo principal sigue en ejecutar_imputacion.py: preparación → evaluación con máscaras → candidatos reales → estadísticas. Se conserva una clase porque sus modelos y escaladores deben reutilizarse. No se agregó un framework, jerarquía de clases ni sistema de plugins.

Se reescribieron todos los scripts de producción y se ajustaron las pruebas y el notebook. Algunas partes ocupan más líneas: separar operaciones mejora su lectura, aunque no reduzca el tamaño del archivo. Los nombres de columnas de salida se conservaron para mantener compatibilidad.

## Verificación frente a la versión 7

- 29 CSV de resultados (15 Ishii, 14 MMC2) idénticos byte a byte: predicciones, candidatos, MAE, cobertura, matrices y estadísticas.
- Máscaras y registros de ajuste idénticos, exceptuando las huellas de archivos que necesariamente cambiaron.
- 10 CSV de normalización y 4 CSV descriptivos numéricamente equivalentes (rtol=1e-12, atol=1e-14).
- 20.483 filas de concentraciones y la tabla de conversión exactamente iguales tras reconstruir desde las fuentes.
- Pasaron 16 pruebas: las 15 previas del motor, exclusiones y máscaras, más una comprobación de los CSV contra las huellas anteriores.
- Las celdas del notebook se ejecutaron secuencialmente con Python para Ishii, MMC2 y una selección sin casos comunes aceptados. Se validó el formato ipynb. El entorno impidió iniciar un kernel Jupyter por sus permisos de sockets; no se afirma haber probado la interfaz de VS Code.
- Ruff comprueba reglas E/F/I en los scripts, y E/F en el notebook; Black comprueba formato de 79 columnas.

El detalle verificable queda en VERIFICACION_REFACTOR.json. Esta comprobación garantiza continuidad del comportamiento implementado, no exactitud biológica de los huecos reales. No se optimizaron métodos ni se alteraron reglas de soporte durante esta refactorización.

## Actualizar una carpeta anterior

Reemplazar los scripts Python, tests/, docs/, README.md, requirements.txt, MANIFEST_PROYECTO.json y los resultados entregados. Agregar apoyo.py, config_imputacion.json y el notebook incluido.

Eliminar config_imputacion_ishii.json, config_imputacion_mmc2.json, requirements_imputacion.txt y requirements_descriptivos.txt. Si quedan de versiones más antiguas, retirar ejecutar_estudio.py, config_validacion.json, resultados_estudio/ y validacion_comun/.

Conservar .git y trabajos personales ajenos a esta distribución. Las fuentes, valores originales y selección de muestras no cambiaron. El manifiesto raíz manifest.json conserva la procedencia de la conversión original; no se presenta como una nueva ejecución.
