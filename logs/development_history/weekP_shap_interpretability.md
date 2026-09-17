# Semana P — TreeSHAP interpretability (Sec. 9) — versión robusta

## Modelo
XGBoost + **EEP_percentile** (oficial Bloque 3). No usar control efficacy_score (Semana O).

## Scope SHAP
1. **Oficial:** LOCO OOF TreeSHAP (`pred_contribs`), n=19,446.
2. **Chequeo:** test chr21/22 n=490.
3. Ranking LOCO vs test490: Spearman **0.982** → test490 no distorsionaba el orden; se reporta LOCO como oficial.

## Importancia global (LOCO OOF, top 10)
| Feature | mean\|SHAP\| |
|---------|-------------|
| guide_3prime_pos4_C | 2.38 |
| ATAC_pm100_mean | 1.94 |
| guide_length | 1.77 |
| mfe_rnafold | 1.68 |
| gc_content | 1.61 |
| seed_pwm_score | 1.48 |
| g_run_max | 1.40 |
| log_dist_nearest_expressed_tss | 1.28 |
| ATAC_pm100_p90 | 0.96 |
| guide_3prime_pos3_C | 0.94 |

Σ|SHAP| por modalidad: **F1 19.5 | F2 5.4 | F5 3.1 | F3 1.5 | F4 0.01** — idéntico orden al ablation N.

## I1 — top 100 EEP (Sec. 9.2, no anécdotas)
Definición: sum SHAP de features cromatina F2 (ATAC_pm100/500/p90, H3K27ac_pm500, dnase); F5 (log_dist_*, ccre_PLS); F1 (todas seq).

| Patrón | Fracción |
|--------|----------|
| Cromatina SHAP sum > 0 | **57/100** |
| F5 SHAP sum > 0 | **60/100** |
| Sequence-dominated (seq > chromatin & seq > F5) | **61/100** |
| Cromatina+ y no seq-dominated | **26/100** |

**Interpretación vs Sec. 9.2:** el proposal esperaba que alta eficacia se explique *primariamente* por ATAC/marcas activas/TSS. Solo **~26%** cumple cromatina+ sin dominio de secuencia; **61%** son seq-dominated (EIF3D-like). Esto **no** se trata como bug de features: es coherente con ablation/RQ2 (F1≫F2) y debe reportarse como hallazgo de heterogeneidad de mecanismos, no como fallo I1 binario pass/fail.

## I2 H3K9me3
SHAP ≈ 0 en casos high-H3K9me3 (test); F3 grupal nulo en LOCO — sin señal de penalización heterocromatina aprendida.

## I3
No ejecutable (Condition B bloqueada).

## Coherencia ablation N
Confirmada a nivel global (modalidades) e individual (F3/F4 residuales). Ninguna feature F3/F4 con importancia alta que contradiga el ablation grupal.

## Artefactos
- `logs/weekP_shap_global_importance_loco.tsv`
- `logs/weekP_shap_importance_loco_vs_test.tsv`
- `logs/weekP_I1_top100.tsv`
- `logs/weekP_plots/`
