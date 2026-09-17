
## Benchmarking interno — comparación de niveles de features

Se completó el benchmarking interno de Phase 1 utilizando el mismo subconjunto Horlbeck2016 CRISPRi, target `efficacy_score` y validación cruzada `GroupKFold(n_splits=5)` agrupada por `nearest_tss_gene`.

### Resultados completos

| Nivel de features | Modelo | Spearman ρ | R² | RMSE |
|---|---|---:|---:|---:|
| Position-only | Ridge | 0.05996 | 0.00057 | 0.43389 |
| Position-only | Random Forest | 0.09056 | -0.07854 | 0.45072 |
| Sequence-only (F1) | Ridge | 0.21021 | 0.04788 | 0.42348 |
| Sequence-only (F1) | Random Forest | 0.21601 | 0.01979 | 0.42969 |
| Multi-modal (F1+F2+F3+F5) | XGBoost | 0.32153 | 0.10042 | 0.41162 |

### Comparación metodológica

El benchmarking muestra una progresión clara al incorporar distintos niveles de información. El modelo **Position-only**, basado únicamente en la distancia al TSS y la condición de estar dentro del promotor de 2 kb, presenta una señal predictiva débil (Spearman ρ ≈ 0.06–0.09), indicando que la posición genómica por sí sola explica una fracción limitada de la variabilidad en la eficacia de los sgRNA. Al incorporar características **Sequence-only (F1)** —GC content, MFE de RNAfold, longitud de guía, máximo de G consecutivas y conteo de T— la correlación aumenta hasta aproximadamente ρ = 0.21, mostrando que las propiedades intrínsecas de la secuencia aportan información sustancialmente mayor. Finalmente, el modelo **Multi-modal (F1+F2+F3+F5)** alcanza ρ = 0.322 y el menor RMSE (0.4116), evidenciando que la combinación de características de secuencia, accesibilidad de cromatina, marcas epigenéticas y contexto/promotor proporciona información complementaria que mejora la predicción respecto a cualquier nivel individual evaluado.

### Conclusión

El benchmarking interno respalda la hipótesis de Phase 1 de que la eficacia de los sgRNA no está determinada únicamente por su posición genómica ni por características intrínsecas de la secuencia. La integración multimodal proporciona la mejor capacidad predictiva entre los modelos evaluados, con una mejora de Spearman de aproximadamente +0.105 frente a Sequence-only y +0.231 frente al mejor modelo Position-only.


## External benchmarking — RuleSet3

An external sequence-only benchmark was performed against RuleSet3 (`crisprScore`) using the canonical 20-nt Horlbeck2016 CRISPRi guides.

- Horlbeck2016 canonical guides: 5,072
- Exact hg38 sequence-context matches: 5,049
- Final 1:1 matched observations: 5,049
- Internal model: XGBoost F1+F2+F3+F5, evaluated using reconstructed 5-fold GroupKFold OOF predictions
- External baseline: RuleSet3 sequence model (Hsu2013 tracrRNA)

### Spearman correlation

| Model | N | Spearman ρ |
|---|---:|---:|
| XGBoost F1+F2+F3+F5 | 5,049 | 0.26731 |
| RuleSet3 | 5,049 | 0.27176 |

Observed difference (RuleSet3 − XGBoost):

- Δρ = +0.00445
- Paired bootstrap: 10,000 replicates, seed 42
- 95% CI: [-0.03055, 0.03891]
- p = 0.8128

The confidence interval includes zero and the paired bootstrap test provides no evidence of a statistically significant difference between RuleSet3 and the internal XGBoost model on this matched benchmark subset.

R² and RMSE were also calculated, but are not used as the primary external comparison because RuleSet3 predictions are on a different score scale and are not calibrated to the Horlbeck2016 CRISPRi activity-score scale.

**Conclusion:** On the 5,049-guide matched external benchmark, RuleSet3 and the internal multimodal XGBoost model show comparable rank-order performance. The observed +0.00445 Spearman advantage of RuleSet3 is small and not statistically significant.
