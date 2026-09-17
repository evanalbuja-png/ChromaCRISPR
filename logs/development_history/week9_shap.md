# Week 9 — SHAP Analysis

## Model

- Model: XGBoost
- Feature set: F1+F2+F3+F5
- n_estimators: 500
- max_depth: 3
- learning_rate: 0.05
- min_child_weight: 5
- random_state: 42
- Cross-validation: 5-fold GroupKFold
- Group variable: nearest_tss_gene

## TreeSHAP validation

TreeSHAP contributions were calculated using XGBoost native
`Booster.predict(..., pred_contribs=True)`.

The SHAP contributions plus the bias term reconstructed the model
predictions with a maximum numerical error below 1.4e-6 across all folds.

## Global SHAP importance

Global feature importance was calculated as mean absolute SHAP value
across the test predictions from all five folds.

| Rank | Feature | Mean absolute SHAP |
|---:|---|---:|
| 1 | nearest_tss_distance | 0.042178 |
| 2 | ATAC_mean | 0.034260 |
| 3 | guide_length | 0.031998 |
| 4 | gc_content | 0.031125 |
| 5 | g_run_max | 0.029731 |
| 6 | mfe_rnafold | 0.029447 |
| 7 | ATAC_max | 0.019390 |
| 8 | poly_t_count | 0.013710 |
| 9 | H3K27me3_mean | 0.012626 |
| 10 | H3K27ac_mean | 0.011861 |
| 11 | H3K4me3_mean | 0.011341 |
| 12 | H3K27ac_max | 0.010262 |
| 13 | H3K27me3_p90 | 0.009854 |
| 14 | ATAC_p90 | 0.009304 |
| 15 | H3K27ac_p90 | 0.008957 |

## Output

- `results/xgb_shap_feature_importance.csv`

## SHAP Dependence Plots

Dependence plots were generated for the following features:

- nearest_tss_distance
- ATAC_mean
- ATAC_max
- H3K27ac_mean
- mfe_rnafold

Outputs:

- `results/shap_dependence_nearest_tss_distance.png`
- `results/shap_dependence_ATAC_mean.png`
- `results/shap_dependence_ATAC_max.png`
- `results/shap_dependence_H3K27ac_mean.png`
- `results/shap_dependence_mfe_rnafold.png`
