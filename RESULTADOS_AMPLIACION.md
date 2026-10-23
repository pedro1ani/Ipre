# Resumen de la ampliación

Se implementaron los cuatro puntos solicitados: diez metabolitos con trazabilidad,
selección interna de parámetros, incertidumbre y análisis de distribuciones.

| Base | Mediana: macro MAE_R | KNN: macro MAE_R | Iterative: macro MAE_R | Menor error |
|---|---:|---:|---:|---|
| Ishii | 0,12045 | 0,10229 | 0,13331 | KNN |
| MMC2 | 0,54706 | 0,62937 | 0,76357 | Mediana |

El error se calcula sobre 228 mediciones ocultadas de Ishii y 68 de MMC2;
cada método responde los mismos casos. El promedio da igual peso a cada
metabolito, en ratios respecto de su referencia Gerosa. No se comparan
magnitudes de error entre bases como si fueran la misma población.

KNN es la primera opción global para Ishii bajo este protocolo. No gana en
todos los metabolitos: la mediana tiene menor MAE_R en fumarato, PEP, s7p y
α-cetoglutarato; Iterative lo tiene en citrato. Estas diferencias por metabolito
son exploratorias, no una segunda selección validada independientemente.
En MMC2, con solo ocho condiciones, la mediana es la referencia más defendible
globalmente en esta evaluación. No hay un ganador universal.

Se generaron candidatos para 72 huecos reales Ishii y dos MMC2, por cada uno de
los tres métodos. Las mediciones no fueron alteradas. ATP, ADP y AMP amplían el
panel a diez y están completos en las 30 muestras Ishii.

## Parámetros

La selección anidada mantiene fuera los grupos de evaluación. En Ishii se eligió
KNN con cinco vecinos y pesos uniformes en cuatro de cinco particiones; en una,
tres vecinos con pesos por distancia. Iterative eligió la configuración base en
tres particiones, la de convergencia más exigente en una y priors 0,01 en una.
En MMC2 las cinco particiones de Iterative escogieron priors 0,01. Esto no prueba
que esos valores sean óptimos para futuras bases.

Iterative solo alcanzó su criterio de parada en los ajustes que generan el
65,4 % de las predicciones externas Ishii; en MMC2 fue el 100 %. El ranking incluye
las predicciones finitas no convergentes y permite identificarlas. Convergencia
no equivale a precisión. El detalle queda en `diagnostico_ajustes.csv`.

## Incertidumbre

Se guardaron 30 reajustes posteriores por caso Iterative, con media, SD y
cuantiles. Los intervalos entre cuantiles 2,5–97,5 % contienen 198/228 mediciones
Ishii (86,8 %) y 52/68 MMC2 (76,5 %). **No están calibrados para prometer una
cobertura del 95 %.** En particular, muestras fuera del rango aprendido y el
recorte a cero pueden producir intervalos poco representativos.

Los tres métodos tienen además rangos empíricos de residuos internos cuando
hay al menos cinco muestras disponibles. Tampoco garantizan cobertura nominal.
Las incertidumbres no incluyen volumen intracelular ni referencia Gerosa.

## Distribuciones

Hay histogramas para los diez metabolitos y tres métodos, junto con media,
mediana, SD, IQR, CV y distancia de Wasserstein. Antes y después comparten
observaciones: no se aplican pruebas para muestras independientes. Se entregan
pruebas pareadas exploratorias sobre medido–estimado, corrección de Holm y sus
supuestos. Un p-valor no significativo no demuestra equivalencia ni calidad.

## Verificación

- 22 pruebas de código pasaron: 16 previas y seis de la ampliación.
- Se comprobaron particiones sin grupos compartidos, valores reales de
  evaluación, conservación de observados y cuantiles/SD de los sorteos guardados.
- 126 celdas adicionales de fuentes contrastadas con Excel (antes de seleccionar
  la cohorte): 102 Ishii y 24 MMC2.
- Se comprobaron formato del notebook y sus celdas secuencialmente con Python
  para ambas bases. El kernel Jupyter no pudo iniciarse por una restricción de
  red local del entorno; no se afirma haber probado la interfaz de VS Code.

La metodología nueva usa bloques externos de grupos y una celda ocultada por
consulta. No sustituye la evaluación previa por series y patrones múltiples;
esta sigue disponible en la versión anterior conservada dentro de la carpeta.
Los cambios de resultados no deben atribuirse solo a los parámetros nuevos.
