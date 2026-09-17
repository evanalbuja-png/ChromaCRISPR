# Week 9 — Position-only Baseline

## Dataset

- Dataset: Horlbeck2016
- Target: `efficacy_score`
- Predictor: `nearest_tss_distance`
- Samples: efficacy_score non-null subset
- Groups: `nearest_tss_gene`
- Cross-validation: 5-fold GroupKFold

## Models

- Ridge Regression (`alpha=1.0`)
- Random Forest (`n_estimators=500`, `random_state=42`, `n_jobs=-1`)

## Results

| Model | Spearman ρ | Pearson R | Pearson R² | R² | RMSE |
|---|---:|---:|---:|---:|---:|
| Random Forest | 0.090558 | 0.082000 | 0.006925 | -0.078538 | 0.450720 |
| Ridge | 0.059964 | 0.037290 | 0.001683 | 0.000573 | 0.433887 |

## Comparison with best complete model

Best complete model: XGBoost F1+F2+F3+F5, Spearman ρ ≈ 0.322.
