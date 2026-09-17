# Week U — Bootstrap 95% CIs (final, preprint-ready)

**N_boot:** 1000  
**Methods:** (i) resample (y, pred) pairs for OOF Spearman and Condition A; (ii) resample fold-level ρ for LOCO-mean and ablation steps.

## 1. LOCO K562 (XGB 101 features)

| Target | Protocol | Ready-to-cite |
|--------|----------|---------------|
| EEP | OOF pair bootstrap | **0.374 [0.362, 0.386]** |
| EEP | mean of 23 fold ρ | **0.375 [0.360, 0.389]** |
| z-score continuo | OOF pair bootstrap | **0.430 [0.418, 0.442]** |
| z-score continuo | mean of 23 fold ρ | **0.429 [0.412, 0.447]** |

CIs de z y EEP **no se solapan** (EEP hi 0.386 < z lo 0.412 en fold-mean; OOF: 0.386 < 0.418).

## 2. Condition A (ρ(pred, −LFC), n=100,942)

| Target | HT29 | A375 |
|--------|------|------|
| EEP | **0.008 [0.002, 0.014]** | **0.012 [0.006, 0.018]** |
| z-score continuo | **0.136 [0.129, 0.142]** | **0.158 [0.152, 0.164]** |

CIs z vs EEP **no se solapan** en ninguna cohorte (evidencia formal de superioridad del target continuo en transferencia).

## 3. Ablation secuencial (target EEP, Semana N) — IC sobre media de 23 folds

| Step | n_features | Ready-to-cite |
|------|------------|---------------|
| F1 | 36 | 0.329 [0.318, 0.339] |
| F1+F2 | 57 | 0.360 [0.348, 0.373] |
| F1+F2+F3 | 90 | 0.360 [0.348, 0.374] |
| F1+F2+F3+F4 | 92 | 0.358 [0.347, 0.370] |
| FULL | 101 | 0.376 [0.364, 0.389] |

Notas: F1→F1+F2 y último tramo hacia FULL son los saltos con p_perm=0 en Semana N; CIs de pasos F2→F4 se solapan (delta≈0), coherente con p_perm no significativo.

## 4. Benchmarking

| Method | Ready-to-cite | Notes |
|--------|---------------|-------|
| XGB EEP (OOF) | **0.374 [0.362, 0.386]** | |
| XGB EEP (LOCO mean-fold) | **0.375 [0.360, 0.389]** | |
| Rule Set 3 (global vs EEP) | **0.255 [0.229, 0.279]** | n=5072 Horlbeck w/ PAM context |
| Rule Set 3 (LOCO mean-fold) | **0.259 [0.232, 0.286]** | 23 folds |
| Position-baseline linear (LOCO mean-fold) | **0.021 [0.005, 0.035]** | 23 folds; ρ≈0.02 |
| Position-baseline nonlinear | ~0.12 (Semana Q) | sin tabla fold-level de ρ GBR en results/; citar como punto estimado previo si no se recompute |

**Overlap:** XGB EEP CI vs RS3 CI → **NO overlap** (XGB lo 0.362 > RS3 hi 0.279). XGB vs position linear → **NO overlap**.

## Archivos fuente

- `results/xgb_oof_predictions_eep_m4.csv`
- `results/xgb_oof_predictions_zscore_continuous.csv`
- `results/weekS_conditionA_zscore_vs_eep_preds.csv`
- `logs/weekN_ablation_seq_folds.tsv`
- `results/weekQ_rs3_loco.csv`, `results/ruleset3_horlbeck_scored.csv`
- `results/weekQ_position_baseline_loco.csv`
