# Week 8 — Baseline Models sobre Horlbeck2016

## Dataset

- Input: `data/processed/sgRNA_feature_matrix_phase1_with_efficacy.csv`
- Filas con `efficacy_score`: 18,386
- Features F1 + F2 + F3 + F5: 23
- Train/test: 80/20
- Random state: 42
- Train rows: 14,708
- Test rows: 3,678

## Features utilizadas

- `gc_content`
- `mfe_rnafold`
- `guide_length`
- `g_run_max`
- `poly_t_count`
- `ATAC_mean`
- `ATAC_max`
- `ATAC_p90`
- `ATAC_sum`
- `H3K27ac_mean`
- `H3K27ac_max`
- `H3K27ac_p90`
- `H3K27ac_sum`
- `H3K4me3_mean`
- `H3K4me3_max`
- `H3K4me3_p90`
- `H3K4me3_sum`
- `H3K27me3_mean`
- `H3K27me3_max`
- `H3K27me3_p90`
- `H3K27me3_sum`
- `nearest_tss_distance`
- `within_promoter_2kb`

## Modelos baseline

| Model | Spearman rho | Pearson R² | RMSE |
|---|---:|---:|---:|
| Ridge Regression | 0.261910 | 0.067347 | 0.412126 |
| LASSO | 0.256399 | 0.062422 | 0.413213 |
| Random Forest | 0.388259 | 0.124984 | 0.399188 |

## Configuración

- Ridge: `alpha=1.0`, StandardScaler + median imputation.
- LASSO: `alpha=0.01`, StandardScaler + median imputation.
- Random Forest: `n_estimators=200`, `random_state=42`, `n_jobs=-1`.
- Métrica de correlación: Spearman rho.
- Métrica de ajuste: R² de Pearson mediante `r2_score`.
- Error: RMSE.
