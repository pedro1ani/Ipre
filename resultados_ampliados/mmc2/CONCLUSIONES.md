# Resultados ampliados: MMC2

El menor macro MAE_R externo corresponde a **Mediana**.

```text
              MAE_R    MAE_mM
method                       
Iterative  0.763574  0.829577
KNN        0.629370  0.693246
Mediana    0.547058  0.639915
```

Cada metabolito tiene el mismo peso. El ranking es exploratorio; no demuestra superioridad universal ni valida los huecos reales. La selección de parámetros se hizo dentro del entrenamiento de cada partición.

Iterative: 100.0 % de las predicciones externas proceden de ajustes que cumplieron el criterio de parada. El ranking incluye todas las predicciones finitas, también las no convergentes; este diagnóstico debe acompañar su interpretación.

## Incertidumbre

Los cuantiles posteriores contienen 52 de 68 mediciones ocultadas (76.5 %). La cobertura observada no alcanza el 95 % nominal de los cuantiles.

```text
   method metabolite  interval  n  coverage  mean_width_mM
Iterative     adp[c]  residual  7  0.571429       0.316370
Iterative     adp[c] posterior  7  1.000000       0.410047
Iterative     akg[c]  residual  7  0.285714       0.456183
Iterative     akg[c] posterior  7  0.428571       0.307733
Iterative     amp[c]  residual  7  0.857143       0.297163
Iterative     amp[c] posterior  7  1.000000       0.362859
Iterative     atp[c]  residual  7  0.857143       1.876424
Iterative     atp[c] posterior  7  1.000000       2.447578
Iterative     cit[c]  residual  7  0.857143       4.105478
Iterative     cit[c] posterior  7  0.714286       6.153550
Iterative     fum[c]  residual  7  0.857143       3.284833
Iterative     fum[c] posterior  7  0.428571       1.308119
Iterative     pep[c]  residual  7  0.714286       0.294879
Iterative     pep[c] posterior  7  1.000000       0.405255
Iterative     pyr[c]  residual  2  0.500000       5.358937
Iterative     pyr[c] posterior  6  0.166667       2.022282
Iterative     s7p[c]  residual  7  0.857143       6.164309
Iterative     s7p[c] posterior  7  0.857143       7.893120
Iterative    succ[c]  residual  2  0.500000       1.134226
Iterative    succ[c] posterior  6  1.000000       0.905036
      KNN     adp[c]  residual  7  0.714286       0.250576
      KNN     akg[c]  residual  7  0.714286       0.575210
      KNN     amp[c]  residual  7  0.714286       0.221118
      KNN     atp[c]  residual  7  0.714286       1.441353
      KNN     cit[c]  residual  7  0.714286       4.038667
      KNN     fum[c]  residual  7  0.714286       3.018229
      KNN     pep[c]  residual  7  0.857143       0.523600
      KNN     pyr[c]  residual  2  0.500000       3.030188
      KNN     s7p[c]  residual  7  0.714286       5.945215
      KNN    succ[c]  residual  2  0.500000       0.494662
  Mediana     adp[c]  residual  7  0.714286       0.359514
  Mediana     akg[c]  residual  7  0.714286       0.729586
  Mediana     amp[c]  residual  7  0.714286       0.311568
  Mediana     atp[c]  residual  7  0.857143       2.037150
  Mediana     cit[c]  residual  7  0.857143       3.314472
  Mediana     fum[c]  residual  7  0.714286       2.916283
  Mediana     pep[c]  residual  7  0.857143       0.576399
  Mediana     pyr[c]  residual  2  0.500000       3.230875
  Mediana     s7p[c]  residual  7  0.857143       6.000186
  Mediana    succ[c]  residual  2  0.500000       0.418128
```

Los intervalos residuales son rangos empíricos de errores internos, no intervalos con cobertura garantizada. Con menos de cinco muestras de calibración se omiten. Los intervalos posteriores de Iterative provienen de 30 reajustes con semillas distintas. Sus cuantiles son aproximados; se reporta su cobertura externa real. La SD es de imputación bajo el modelo, no la SD experimental. No incluye la incertidumbre del volumen celular ni de la referencia Gerosa.

## Distribuciones y pruebas

Los histogramas usan los mismos intervalos de clase antes/después. Se reportan media, mediana, SD, IQR, CV y distancia de Wasserstein. No se calcula t-test, Mann–Whitney ni KS independiente sobre antes/después porque las muestras se solapan. Para valores ocultados se reporta distancia entre ECDF y permutación pareada por grupo, junto con prueba de signos de sesgo y corrección de Holm. Estas pruebas son exploratorias: requieren intercambiabilidad o simetría de errores por grupo y los modelos de validación cruzada comparten entrenamiento. Un p alto no demuestra equivalencia de distribuciones ni imputaciones correctas.

## Límites del estudio

- Diez metabolitos, 30 muestras Ishii u ocho condiciones MMC2, por separado.
- Glucosa MMC2 se excluye de evaluación: es la referencia conocida.
- La evaluación nueva oculta una celda por consulta; las ausencias reales se imputan conjuntamente. No sustituye las pruebas de patrones reales y series analíticas conservadas en el flujo anterior.
- Cada candidato real se obtiene sin su bloque externo de grupos. Se usa la misma selección interna; no se reajusta usando la muestra consultada.
- ATP, ADP y AMP están completos en la cohorte Ishii: amplían la evaluación y los predictores, pero no generan huecos reales artificiales.
- La reducción de dispersión no se interpreta como una mejora.
- Los resultados y parámetros nuevos no tienen que coincidir con los antiguos.
