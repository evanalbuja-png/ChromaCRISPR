# Week 9 Follow-up — External Efficacy Dataset Investigation

## Objective

Investigate whether Sanson2018, Replogle2022, and Gasperini2019 contain a usable quantitative efficacy measurement that could be incorporated into ChromaCRISPR Phase 1.

## Results

### Sanson2018

Files inspected:

- `data/interim/sanson2018_recovered_official.csv`
- `data/interim/sanson2018_brunello_reference.csv`
- `data/processed/sanson2018_reference_FINAL.csv`

`data/interim/sanson2018_recovered_official.csv` contains the columns:

- `score_raw`
- `score_norm`
- `score_type`

However:

- `score_raw`: 0 non-null values
- `score_norm`: 0 non-null values
- `score_type`: entirely null

Therefore, no usable efficacy score was recovered from the available Sanson2018 files.

### Replogle2022

Files inspected:

- `data/raw/crispr_datasets/Replogle2022/NIHMS1812939-supplement-11.xlsx`
- `data/processed/replogle2022_reference_FINAL.csv`

The original supplementary workbook contains:

- `TabA_K562_day8_library`
- `TabB_K562_day6_library`
- `TabC_RPE1_day7_library`

These tables contain sgRNA pair/library information but no quantitative efficacy measurement.

The processed reference file also contains only guide/coordinate annotation fields and no efficacy score.

No additional Replogle2022 files were found under `data/`.

Therefore, no usable efficacy score was identified in the currently available Replogle2022 data.

### Gasperini2019

File inspected:

- `data/interim/gasperini2019_efficacy.csv`

This file contains a quantitative `efficacy_score`.

Statistics:

- Rows: 1,676
- Unique guides: 1,530
- Non-null `efficacy_score`: 1,676
- Experiments:
  - S1: 258
  - S2: 1,418
- Mean efficacy score: 0.705822
- Standard deviation: 0.182020
- Minimum: 0.024877
- Median: 0.762175
- Maximum: 0.985834

Therefore, Gasperini2019 provides a potentially usable efficacy dataset.

## Conclusion

| Dataset | Usable efficacy identified? |
|---|---|
| Sanson2018 | No |
| Replogle2022 | No |
| Gasperini2019 | Yes |

No integration or model training was performed in this step.

## Validation of Gasperini2019 efficacy score

The provenance of `data/interim/gasperini2019_efficacy.csv` was reviewed using:
- `scripts/extract_horlbeck_gasperini_efficacy_week8.py`
- `scripts/parsers/gasperini.py`
- Original Gasperini2019 supplementary tables.

The original Gasperini tables do not contain a variable named `efficacy_score`.

The extracted `efficacy_score` was created by renaming:
`Diff_expression_test_fold_change` → `efficacy_score`

The association was performed as:
`Spacer → Target_Site → Diff_expression_test_fold_change`

Therefore, multiple sgRNAs targeting the same `Target_Site` inherit the same score. This explains the repeated `guide_sequence` observations in the extracted file.

This variable represents the effect of the perturbation on target-gene expression and is not equivalent to the Horlbeck2016 `CRISPRi activity score`.

### Decision

`Diff_expression_test_fold_change` from Gasperini2019 will NOT be incorporated as an additional `efficacy_score` target in the current ChromaCRISPR Phase 1 training dataset.

Gasperini2019 remains potentially useful as an independent source of enhancer-gene regulatory information, but it will not be treated as an additional guide-efficacy dataset.

No dataset or model was modified as part of this validation.

## Step 2 — Inspección de efficacy en otros datasets
- Inspección ejecutada para Sanson2018, Replogle2022 y Gasperini2019.
- Se revisaron archivos, hojas Excel, columnas, filas, columnas numéricas, candidatos de efficacy/activity/phenotype/log2FC/depletion y posibles columnas de unión.
- Resultados detallados: salida del comando de inspección ejecutado en esta sesión.
