# Week S — Target redesign: EEP vs z-score continuo armonizado

**Fecha:** 2026-09-14  
**Estado:** CERRADO  
**Contexto:** Bloque 5 extendido. Hypothesis (Semana O): el colapso de Condition A se debe al target EEP (percentil within-gene), no solo a dominio de features.

## Diseño del target continuo

- Columna: `z_score_within_study` (ya en `K562_training_labels_harmonized_v1.csv` / matriz v1), documentada como **`z_score_harmonized_continuous`**.
- Origen: z-score within-study por dataset (Horlbeck, Gasperini) previo al percentil within-gene (Semana D). **No** se recalcula EEP; ambos targets coexisten.
- Overlap de secuencias Horlbeck∩Gasperini: 0 (Semana D); sin conflicto de duplicados cross-study.
- `EEP_percentile` **no modificado**.

## Protocolo

- Modelo: XGBoost, **mismos hiperparámetros Semana L / M4** (101 features):
  - max_depth=4, lr≈0.0297, n_estimators=289, min_child_weight=4,
  - subsample≈0.766, colsample_bytree≈0.567, reg_lambda≈0.318
- LOCO: 23 folds (chr1–22 + X), idéntico a Semana M.
- Condition A: train en las 19 446 guías K562; predict sobre `data/interim/conditionA_full/sanson_features_full.csv` (101/101 features); merge con LFC Sanson FIXED.
- Polaridad Condition A: se reporta **ρ(pred, −LFC)** (predicción alta = más knockdown = LFC más negativo).
- Permutation: 10 000 iteraciones sobre Δ = ρ_z − ρ_EEP (vs −LFC).

## Resultados principales

### 1. LOCO K562 (mismo feature set 101)

| Target | LOCO mean ρ | SD | Median |
|--------|-------------|-----|--------|
| EEP_percentile (Semana M) | 0.377 | 0.039 | 0.367 |
| **z_score continuo** | **0.429** | **0.044** | **0.412** |

→ Dentro de K562 el target continuo **mejora** (~+0.05), no hay trade-off de ranking interno.

Archivos: `results/xgb_oof_predictions_zscore_continuous.csv`, `results/weekS_loco_zscore_folds.csv`, `results/weekS_loco_zscore_summary.json`.

### 2. Condition A — transferencia Sanson2018 (n=100 942)

| Target | HT29 ρ(pred,−LFC) | A375 ρ(pred,−LFC) |
|--------|-------------------|-------------------|
| EEP (re-entrenado simétrico, params M4) | 0.008 | 0.012 |
| **z-score continuo** | **0.136** | **0.158** |

Permutation test (H1: z ≥ EEP):

| Cohort | Δ (z−EEP) | p_perm |
|--------|-----------|--------|
| HT29 | +0.127 | **≈0.0001** |
| A375 | +0.146 | **≈0.0001** |

Referencia previa (Semana O / control Horlbeck `efficacy_score` only): HT29≈0.15, A375≈0.17 — mismo orden de magnitud que z-score armonizado (que además incluye Gasperini).

Archivos: `results/weekS_conditionA_zscore_vs_eep_preds.csv`, `results/weekS_conditionA_zscore_summary.json`.

### 3. Ablation F1 vs FULL (target z-score, LOCO)

| Feature set | n_feat | LOCO mean | SD |
|-------------|--------|-----------|-----|
| F1 (secuencia + PAM one-hot proxy) | 35 | 0.371 | 0.041 |
| FULL | 101 | 0.429 | 0.044 |
| Δ FULL−F1 | — | **+0.058** | — |

Patrón cualitativo alineado con Semana N (F1 domina; cromatina aporta delta positivo estable bajo target continuo).

Archivo: `results/weekS_ablation_f1_vs_full_zscore.json`.

## Interpretación (preprint-ready)

1. **Causa dominante del fallo de transferencia Condition A:** diseño del target EEP (ranking within-gene), no ausencia total de señal de cromatina.
2. **Evidencia causal formal:** misma arquitectura, mismos features, mismos hiperparámetros; solo cambia el target → Δ Condition A altamente significativo (p≈10⁻⁴).
3. **No hay sacrificio en K562:** LOCO sube de 0.377 a 0.429.
4. **Limitación:** ρ de transferencia (~0.14–0.16) sigue siendo modesto en valor absoluto; el rediseño **rescata** señal respecto de EEP≈0, no resuelve del todo el gap cross-library/cell-line.
5. **RS2/DeepCRISPR** permanecen bloqueantes (Semana R); no reabiertos aquí.

## Tabla resumen para figura principal

| Modelo / target | LOCO K562 | Cond.A HT29 | Cond.A A375 |
|-----------------|-----------|-------------|-------------|
| XGB 101 + EEP | 0.377±0.039 | ~0.01 | ~0.01 |
| XGB 101 + z-score continuo | **0.429±0.044** | **0.136** | **0.158** |
| Control efficacy_score Horlbeck-only (prev.) | — | ~0.15 | ~0.17 |


## Sanity checks (pre-Tabla 1)

### 1. ¿El LOCO z está inflado por separabilidad Horlbeck vs Gasperini?

- `z_score_within_study` por origen: Horlbeck y Gasperini ambos mean≈0, sd≈1 (z within-study). No hay desalineación de escala entre estudios.
- Las 101 features **sí** predicen `dataset_origin` (RF depth=5, AUC≈1.0, ACC≈0.997 vs majority 0.942). Esto es esperable: Gasperini es mayormente distal/enhancer; Horlbeck promotor-proximal — F2/F3/F4/F5 codifican ese contexto genómico, no una fuga del label.
- **Test decisivo:** LOCO z **solo Horlbeck** (mismos params M4, 23 folds) = **0.448 ± 0.041**, *mayor* que LOCO H+G = 0.429 ± 0.044.
  - Si el LOCO full estuviera inflado por “detectar Gasperini”, el LOCO solo-Horlbeck sería menor. Ocurre lo contrario: añadir Gasperini diluye ligeramente el LOCO K562.
- **Conclusión:** se descarta inflación artificial del ρ=0.429 por heterogeneidad de estudio. El AUC de origen refleja separabilidad biológica de locus, no un artefacto del target continuo.

### 2. Comparación explícita vs control Semana O (Condition A)

| Setup | Train target | Train set | HT29 ρ(pred,−LFC) | A375 |
|-------|--------------|-----------|-------------------|------|
| Semana O control | `efficacy_score` Horlbeck (no armonizado) | Solo Horlbeck | ~0.15 | ~0.17 |
| Semana S | z-score continuo armonizado | Horlbeck + Gasperini | **0.136** | **0.158** |
| Semana S EEP (simétrico) | EEP_percentile | Horlbeck + Gasperini | 0.008 | 0.012 |

El resultado z-score armonizado es **ligeramente menor** que el control Semana O pese a más datos de entrenamiento. Hipótesis (no demostradas aquí): (i) bajo peso relativo de Gasperini (~6% de filas); (ii) heterogeneidad de ensayo (proportion_remaining transformado vs scores Horlbeck) que añade ruido al target combinado sin aportar mucha señal transferible a Sanson; (iii) el control O pudo beneficiarse de alineación más directa con la escala de efficacy Horlbeck. En cualquier caso, ambos superan de forma clara a EEP≈0; el rediseño de target (continuo vs percentil within-gene) sigue siendo el factor dominante frente a “añadir o no Gasperini”.

