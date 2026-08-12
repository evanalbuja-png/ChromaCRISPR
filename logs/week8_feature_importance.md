# Week 8 — Random Forest Feature Importance

## Dataset

- Input: `data/processed/sgRNA_feature_matrix_phase1_with_efficacy.csv`
- Dataset: `Horlbeck2016`
- Target: `efficacy_score`
- Filas utilizadas: 18,386
- Features: 23
- Train/test: 80/20
- Random state: 42
- Random Forest: `n_estimators=200`, `n_jobs=-1`

## Top 10 features

| Rank | Feature | Importance |
|---:|---|---:|
| 1 | `nearest_tss_distance` | 0.10540204 |
| 2 | `mfe_rnafold` | 0.07737575 |
| 3 | `ATAC_max` | 0.06257259 |
| 4 | `ATAC_p90` | 0.05929537 |
| 5 | `ATAC_sum` | 0.05632187 |
| 6 | `ATAC_mean` | 0.05595513 |
| 7 | `H3K27ac_max` | 0.05535327 |
| 8 | `H3K27ac_p90` | 0.05144684 |
| 9 | `gc_content` | 0.05033285 |
| 10 | `H3K4me3_p90` | 0.04018842 |

## Ranking completo

| Rank | Feature | Importance |
|---:|---|---:|
| 1 | `nearest_tss_distance` | 0.10540204 |
| 2 | `mfe_rnafold` | 0.07737575 |
| 3 | `ATAC_max` | 0.06257259 |
| 4 | `ATAC_p90` | 0.05929537 |
| 5 | `ATAC_sum` | 0.05632187 |
| 6 | `ATAC_mean` | 0.05595513 |
| 7 | `H3K27ac_max` | 0.05535327 |
| 8 | `H3K27ac_p90` | 0.05144684 |
| 9 | `gc_content` | 0.05033285 |
| 10 | `H3K4me3_p90` | 0.04018842 |
| 11 | `H3K27me3_max` | 0.03946976 |
| 12 | `H3K4me3_sum` | 0.03748008 |
| 13 | `H3K4me3_mean` | 0.03732307 |
| 14 | `H3K27ac_sum` | 0.03727503 |
| 15 | `H3K27ac_mean` | 0.03639150 |
| 16 | `H3K27me3_p90` | 0.03583725 |
| 17 | `H3K27me3_sum` | 0.03350695 |
| 18 | `H3K27me3_mean` | 0.03344168 |
| 19 | `H3K4me3_max` | 0.03278532 |
| 20 | `guide_length` | 0.02411381 |
| 21 | `g_run_max` | 0.02357466 |
| 22 | `poly_t_count` | 0.01366354 |
| 23 | `within_promoter_2kb` | 0.00089321 |
