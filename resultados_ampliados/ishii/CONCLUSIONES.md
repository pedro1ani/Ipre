# Resultados ampliados: ISHII

El menor macro MAE_R externo corresponde a **KNN**.

```text
              MAE_R    MAE_mM
method                       
Iterative  0.133307  0.193695
KNN        0.102287  0.114105
Mediana    0.120449  0.136317
```

Cada metabolito tiene el mismo peso. El ranking es exploratorio; no demuestra superioridad universal ni valida los huecos reales. La selección de parámetros se hizo dentro del entrenamiento de cada partición.

Iterative: 65.4 % de las predicciones externas proceden de ajustes que cumplieron el criterio de parada. El ranking incluye todas las predicciones finitas, también las no convergentes; este diagnóstico debe acompañar su interpretación.

## Incertidumbre

Los cuantiles posteriores contienen 198 de 228 mediciones ocultadas (86.8 %). La cobertura observada no alcanza el 95 % nominal de los cuantiles.

```text
   method metabolite  interval  n  coverage  mean_width_mM
Iterative     adp[c]  residual 30  0.966667       0.831121
Iterative     adp[c] posterior 30  0.933333       0.519779
Iterative     akg[c]  residual 11  0.818182       0.113863
Iterative     akg[c] posterior 11  0.818182       0.087381
Iterative     amp[c]  residual 30  0.866667       0.901473
Iterative     amp[c] posterior 30  0.866667       0.592099
Iterative     atp[c]  residual 30  0.733333       1.963294
Iterative     atp[c] posterior 30  0.833333       1.485880
Iterative     cit[c]  residual 12  0.833333       0.054792
Iterative     cit[c] posterior 12  0.833333       0.051704
Iterative     fum[c]  residual 23  0.913043       0.104631
Iterative     fum[c] posterior 23  0.869565       0.076077
Iterative     pep[c]  residual 22  0.909091       0.219332
Iterative     pep[c] posterior 22  0.863636       0.183567
Iterative     pyr[c]  residual 14  0.928571       0.369064
Iterative     pyr[c] posterior 14  0.857143       0.297269
Iterative     s7p[c]  residual 29  0.931034       1.705796
Iterative     s7p[c] posterior 29  0.896552       1.146110
Iterative    succ[c]  residual 27  0.925926       0.145188
Iterative    succ[c] posterior 27  0.851852       0.084166
      KNN     adp[c]  residual 30  0.900000       0.609465
      KNN     akg[c]  residual 11  0.636364       0.072231
      KNN     amp[c]  residual 30  0.800000       0.624007
      KNN     atp[c]  residual 30  0.933333       1.756967
      KNN     cit[c]  residual 12  0.833333       0.078858
      KNN     fum[c]  residual 23  0.869565       0.076154
      KNN     pep[c]  residual 22  0.909091       0.235394
      KNN     pyr[c]  residual 14  0.857143       0.304660
      KNN     s7p[c]  residual 29  0.965517       1.776832
      KNN    succ[c]  residual 27  0.925926       0.111197
  Mediana     adp[c]  residual 30  0.933333       0.817905
  Mediana     akg[c]  residual 11  0.818182       0.066941
  Mediana     amp[c]  residual 30  0.866667       0.721510
  Mediana     atp[c]  residual 30  0.900000       2.132855
  Mediana     cit[c]  residual 12  0.833333       0.074024
  Mediana     fum[c]  residual 23  0.913043       0.074639
  Mediana     pep[c]  residual 22  0.954545       0.213393
  Mediana     pyr[c]  residual 14  0.928571       0.323095
  Mediana     s7p[c]  residual 29  0.896552       1.726309
  Mediana    succ[c]  residual 27  0.888889       0.109413
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
