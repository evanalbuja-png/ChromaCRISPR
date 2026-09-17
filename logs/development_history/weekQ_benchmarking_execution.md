# Semana Q — Benchmarking (Sec. 12) — ejecución final

## Position-baseline (Sec. 12.3)
- Predictor: `log_dist_target_tss`
- LOCO (23 folds): **ρ = 0.0206 ± 0.0370**
- Proposal expected ~0.40–0.50: **no se reproduce** en K562 EEP (Horlbeck+Gasperini).
- Interpretación: distancia al TSS del gen target aporta casi cero ranking within/between el set actual; alineado con SHAP (F1/F2 dominan; `log_dist_nearest_expressed_tss` > target TSS dist).

## Rule Set 2 / Azimuth (Sec. 12.1)
- Azimuth: no importable en chromacrispr-phase1 (py3.10) tras intento acotado.
- CRISPOR HTTPS: connection reset; HTTP reachable but not used as primary.
- **Desviación documentada:** Rule Set **3** (`results/ruleset3_horlbeck_scored.csv`).
  - Mapped: 5,056 / 18,318 Horlbeck guides
  - LOCO ρ: **0.2587 ± 0.0654** (23 folds)
  - Global Spearman vs EEP: 0.255
  - canonical_ngg rate in scored file ~0.49

## DeepCRISPR (Sec. 12.2)
- Factibilidad only: requires **Python 3.6 + TensorFlow 1.3.0 + sonnet 1.9**; models in `trained_models/*.tar.gz`; Docker image available.
- **Bloqueante** en este entorno sin contenedor legacy. No install.

## Tabla comparativa final
| Method | LOCO ρ | vs XGB |
|--------|--------|--------|
| XGBoost 101 features | **0.377 ± 0.039** | — |
| Sequence-only (F1) | ~0.329 | −0.05 |
| Rule Set 3 (RS2 proxy) | **0.259 ± 0.065** | −0.12 |
| Position-baseline | **0.021 ± 0.037** | −0.36 |
| DeepCRISPR | N/A | bloqueante |

## Conclusión Sec. 12
ChromaCRISPR (XGB+EEP) supera baselines de secuencia (RS3) y posición en LOCO K562. El claim del proposal de que position-baseline ~0.45 no se sostiene aquí. DeepCRISPR y CRISPOR/Azimuth nativos quedan fuera por infraestructura, no por descarte científico a priori.

## Addendum — diagnóstico position-baseline (post-revisión)

### Qué se calculó
- LOCO con regresión lineal univariada y Spearman de `log_dist_target_tss` vs EEP.
- Replicado con `log_dist_nearest_expressed_tss`, Ridge y GBR (max_depth=2).

### Feature correcta
- Sec. 12.3 / F5: distancia guía→TSS del gen target ≈ `log_dist_target_tss`.
- `log_dist_nearest_expressed_tss` es muy similar (Spearman 0.90, 79% iguales).
- **No** se usó la feature “equivocada”.

### Signo
- Spearman(log_dist, EEP) ≈ −0.02 → mayor distancia, ligeramente peor EEP.
- Coincide con Horlbeck (promotor-proximal más efectivo). No hay inversión de signo.

### Resultados corregidos / ampliados
| Setup | LOCO ρ |
|-------|--------|
| Linear `log_dist_target_tss` | 0.021 ± 0.037 |
| Linear `log_dist_nearest_expressed_tss` | 0.013 ± 0.038 |
| GBR stumps `log_dist_target_tss` | **0.119 ± 0.042** |
| GBR stumps `log_dist_nearest_expressed_tss` | **0.122 ± 0.036** |

Aun el baseline no lineal univariado (~0.12) queda **muy por debajo** de seq-only (0.33), RS3 (0.26) y XGB (0.38), y del rango proposal 0.40–0.50.

### Resolución de la “contradicción” ablation/SHAP
- Ablation F5 = **modalidad completa**, no distancia sola.
- SHAP alto en distancia = efecto **multivariado/interacción**, no poder univariado.
- Target EEP within-gene + librería Horlbeck proximal → poca variación de distancia *útil para ranking univariado*.

**Conclusión Sec. 12.3 (confirmada):** position-only no es baseline competitivo en este diseño experimental; el número ρ≈0.02 (lineal) / ≈0.12 (no lineal) es real, no un artefacto de evaluación.
