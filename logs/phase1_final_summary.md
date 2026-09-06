# ChromaCRISPR Phase 1 — Final Summary

## Phase 1 status

**Phase 1 formally completed.**

Phase 1 established and evaluated a multi-modal machine-learning framework for sgRNA efficacy prediction using sequence, chromatin accessibility, epigenetic and genomic-context features.

---

## 1. Phase 1 objective

The objective of Phase 1 was to establish a reproducible sgRNA efficacy prediction pipeline and determine whether integrating multiple biological feature modalities improves prediction compared with simpler feature sets.

The Phase 1 evaluation focused on the Horlbeck2016 CRISPRi dataset and used `efficacy_score` as the prediction target.

---

## 2. Final modeling dataset

Final dataset:

`data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv`

Dataset characteristics:

- Samples: **18,332 sgRNAs**
- Columns: **31**
- Non-null `efficacy_score`: **18,332**
- Unique `nearest_tss_gene` groups: **1,821**
- Cross-validation: **5-fold GroupKFold**
- Group variable: `nearest_tss_gene`

The final modeling dataset was restricted to Horlbeck2016 for the Phase 1 efficacy benchmark.

---

## 3. Implemented feature modalities

### F1 — Sequence features

Five sequence-derived features were implemented:

- `gc_content`
- `mfe_rnafold`
- `guide_length`
- `g_run_max`
- `poly_t_count`

### F2 — Chromatin accessibility

Four ATAC-seq features were implemented:

- `ATAC_mean`
- `ATAC_max`
- `ATAC_p90`
- `ATAC_sum`

### F3 — Epigenetic marks

Twelve features were implemented from three histone marks:

**H3K27ac**
- `H3K27ac_mean`
- `H3K27ac_max`
- `H3K27ac_p90`
- `H3K27ac_sum`

**H3K4me3**
- `H3K4me3_mean`
- `H3K4me3_max`
- `H3K4me3_p90`
- `H3K4me3_sum`

**H3K27me3**
- `H3K27me3_mean`
- `H3K27me3_max`
- `H3K27me3_p90`
- `H3K27me3_sum`

### F5 — Genomic context

Three genomic-context features were implemented:

- `nearest_tss_gene`
- `nearest_tss_distance`
- `within_promoter_2kb`

F4 (Hi-C) was not included in the final Phase 1 model.

---

## 4. Main ablation results

A complete feature-ablation benchmark was performed with Random Forest using the same 5-fold GroupKFold framework.

| Feature set | Spearman ρ | R² | RMSE |
|---|---:|---:|---:|
| F1 | 0.216015 | 0.019793 | 0.429689 |
| F1+F2 | 0.261788 | 0.063891 | 0.419901 |
| F1+F2+F3 | 0.267895 | 0.072602 | 0.417939 |
| F1+F2+F3+F5 | **0.310199** | **0.094731** | **0.412926** |

The ablation shows a progressive improvement as additional feature modalities are incorporated. The largest improvement occurs after adding F5 to the F1+F2+F3 configuration.

A separate XGBoost ablation table was not generated. Therefore, no XGBoost-specific ablation values are reported here.

---

## 5. Final multi-modal model

The best Phase 1 model was:

**XGBoost — F1+F2+F3+F5**

Configuration:

- `n_estimators = 500`
- `max_depth = 3`
- `learning_rate = 0.05`
- `min_child_weight = 5`
- `random_state = 42`
- 5-fold GroupKFold
- Groups: `nearest_tss_gene`

Performance:

- **Spearman ρ = 0.32153**
- **R² = 0.10042**
- **RMSE = 0.41162**

Five-fold Spearman values:

- 0.331641
- 0.336326
- 0.301200
- 0.327819
- 0.310646

Mean Spearman:

**ρ ≈ 0.322**

---

## 6. Internal benchmarking

The internal benchmark compared three information levels:

| Feature level | Model | Spearman ρ | R² | RMSE |
|---|---|---:|---:|---:|
| Position-only | Ridge | 0.05996 | 0.00057 | 0.43389 |
| Position-only | Random Forest | 0.09056 | -0.07854 | 0.45072 |
| Sequence-only (F1) | Ridge | 0.21021 | 0.04788 | 0.42348 |
| Sequence-only (F1) | Random Forest | 0.21601 | 0.01979 | 0.42969 |
| Multi-modal (F1+F2+F3+F5) | XGBoost | **0.32153** | **0.10042** | **0.41162** |

The benchmark demonstrates that sequence features provide substantially more predictive signal than genomic position alone, while the multi-modal feature combination provides the strongest performance among the evaluated configurations.

---

## 7. SHAP analysis

SHAP analysis was performed on the final XGBoost F1+F2+F3+F5 model using native XGBoost TreeSHAP.

TreeSHAP validation confirmed that SHAP contributions plus the bias term reconstructed model predictions with a maximum numerical error below **1.4e-6** across all folds.

### Top global features

| Rank | Feature | Mean absolute SHAP |
|---:|---|---:|
| 1 | `nearest_tss_distance` | 0.042178 |
| 2 | `ATAC_mean` | 0.034260 |
| 3 | `guide_length` | 0.031998 |
| 4 | `gc_content` | 0.031125 |
| 5 | `g_run_max` | 0.029731 |
| 6 | `mfe_rnafold` | 0.029447 |
| 7 | `ATAC_max` | 0.019390 |
| 8 | `poly_t_count` | 0.013710 |
| 9 | `H3K27me3_mean` | 0.012626 |
| 10 | `H3K27ac_mean` | 0.011861 |

Additional SHAP features and dependence plots are documented in:

- `logs/week9_shap.md`
- `results/xgb_shap_feature_importance.csv`

Dependence plots were generated for `nearest_tss_distance`, `ATAC_mean`, `ATAC_max`, `H3K27ac_mean`, and `mfe_rnafold`.

---

## 8. Additional efficacy dataset investigation

Additional datasets were investigated during Phase 1 follow-up.

### Sanson2018

No usable quantitative efficacy target was recovered. The inspected score fields contained zero non-null quantitative values.

### Replogle2022

Supplementary resources provided guide/library and coordinate information but no quantitative efficacy measurement suitable for integration as an equivalent target.

### Gasperini2019

A processed dataset contained:

- 1,676 rows
- 1,530 unique guides
- 1,676 non-null values
- `efficacy_score` mean: 0.705822
- standard deviation: 0.182020
- minimum: 0.024877
- median: 0.762175
- maximum: 0.985834

However, validation showed that the original variable was `Diff_expression_test_fold_change`, propagated through `Spacer → Target_Site → Diff_expression_test_fold_change`.

Because multiple sgRNAs targeting the same `Target_Site` inherit the same value, and because this measurement represents target-gene expression effect rather than an efficacy score directly equivalent to the Horlbeck2016 CRISPRi activity target, Gasperini2019 was **not incorporated as an additional efficacy target in Phase 1**.

No Phase 1 dataset or model was modified as a consequence.

Detailed investigation:

- `logs/week9_followup.md`

---

## 9. Limitations

The following limitations are recognized at the end of Phase 1:

1. The final efficacy modeling benchmark relies on Horlbeck2016 as the quantitative efficacy dataset.
2. Additional datasets investigated during Phase 1 did not provide directly equivalent quantitative efficacy measurements suitable for straightforward target integration.
3. F4 Hi-C features were not incorporated into the final model.
4. The best observed Spearman correlation (~0.322) indicates meaningful but limited predictive performance; the model does not provide highly accurate individual sgRNA efficacy prediction.
5. Internal benchmarking does not constitute external validation on an independent efficacy dataset.
6. The SHAP results describe feature contributions to the trained model and should not be interpreted as direct evidence of biological causality.
7. The feature space remains limited to the modalities implemented in Phase 1.

---

## 10. Phase 2 backlog

The following items are explicitly deferred to Phase 2:

### Hi-C integration (F4)

Hi-C was prepared during the project but was not included in the final Phase 1 model. Its integration requires additional feature engineering and validation before it can be evaluated fairly against the established baseline.

### Additional efficacy datasets

Sanson2018, Replogle2022 and Gasperini2019 should not be forced into the Phase 1 target simply to increase sample size. Future work can investigate whether alternative quantitative formulations or task-specific modeling make these datasets compatible.

### External validation

Future evaluation should test the model on an independent dataset with a quantitatively compatible efficacy endpoint.

### Additional guide-scoring resources

CRISPOR and DeepCRISPR scores were not available in a form suitable for the Phase 1 benchmark and remain candidates for future feature expansion.

### Dataset expansion

Sanson and other datasets may be revisited once a defensible harmonization strategy for efficacy targets is established.

---

## 11. Final Phase 1 conclusion

Phase 1 successfully established a reproducible multi-modal sgRNA efficacy prediction framework and completed the planned internal evaluation, including:

- Feature engineering for F1, F2, F3 and F5
- Final Horlbeck2016 modeling dataset
- Position-only baseline
- Sequence-only benchmark
- Random Forest feature ablation
- XGBoost multi-modal modeling
- Internal benchmarking
- TreeSHAP validation and global feature importance
- Investigation of additional efficacy datasets
- Documentation of limitations and Phase 2 items

The best Phase 1 model is the **XGBoost F1+F2+F3+F5 multi-modal model**, achieving **Spearman ρ ≈ 0.322** under 5-fold GroupKFold grouped by `nearest_tss_gene`.

**Phase 1 is complete and ready for formal closure.**
