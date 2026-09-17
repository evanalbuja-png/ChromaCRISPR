# Semana M — LOCO 23 folds (Sec. 10.1)

**Modelo primario:** Ridge (hiperparámetros fijos Semana L; sin re-tune por fold).  
**Referencia:** XGBoost (params Optuna 100 trials).  
**Folds:** chr1–22 + X (Y ausente en matriz).  
**Features:** 101 accepted. Target: EEP_percentile.

## Métrica oficial del proyecto (de aquí en adelante)
| Model | LOCO mean ρ | SD |
|-------|-------------|-----|
| **Ridge** | **0.3025** | **0.0374** |
| XGBoost | 0.3770 | 0.0385 |

## Por cromosoma
Ver `logs/weekM_loco_folds.tsv`.

## Comparaciones
| Protocolo | Ridge ρ | XGB ρ |
|-----------|---------|-------|
| Val chr20 solo (Semana L) | 0.335 | 0.403 |
| **LOCO 23-fold mean** | **0.303** | **0.377** |
| Condition A GroupKFold ~5 (23 feats, Phase1) | ~0.28–0.30 | — |

LOCO Ridge ≈ val chr20 dentro de 1 SD. Ganancia XGB vs Condition A antigua atribuible en parte al feature set 101.

## Outliers
- Ridge: ninguno (|z|>2).
- XGB: chr21 (ρ=0.47, n=199) — N bajo; no estructural (Spearman ρ vs n_test NS).

## Decisión
- Primario sigue siendo **Ridge** (Semana L).
- Métrica de performance reportable: **ρ_LOCO = 0.303 ± 0.037** (Ridge).
- XGB LOCO 0.377 ± 0.039 disponible como techo de ranking para Condition A/B si se prioriza transferencia por performance.
