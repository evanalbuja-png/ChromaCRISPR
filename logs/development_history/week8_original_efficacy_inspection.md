# ChromaCRISPR Phase 1 — Week 8

## Paso 5A — Inspección de datasets originales de eficacia

Se realizó una búsqueda dirigida de los datasets Horlbeck2016, Gasperini2019 y Replogle2022.

**No se realizó ningún merge ni modificación de la matriz Phase 1.**

## Horlbeck2016

### `NONE`

- Candidatos de score: `NONE`
- Candidatos de llave: `NONE`

## Gasperini2019

### `data/interim/gasperini_coordinate_validation_sample.csv`

- Candidatos de score: `NONE`
- Candidatos de llave: `guide_sequence`

### `data/interim/gasperini_validation_ensembl_intersect.txt`

- Candidatos de score: `READ_ERROR`
- Candidatos de llave: `NONE`

### `data/interim/gasperini_validation_intersect.txt`

- Candidatos de score: `READ_ERROR`
- Candidatos de llave: `NONE`

### `data/interim/gasperini_validation_window10kb_ensembl_intersect.txt`

- Candidatos de score: `NONE`
- Candidatos de llave: `gene_id "ENSG00000188460"; gene_version "4"; gene_name "ACTBP11"; gene_source "havana"; gene_biotype "pseudogene";`

### `data/interim/gasperini_validation_window10kb_intersect.txt`

- Candidatos de score: `READ_ERROR`
- Candidatos de llave: `NONE`

## Replogle2022

### `NONE`

- Candidatos de score: `NONE`
- Candidatos de llave: `NONE`

## Estado

- Sanson2018 ya posee `score_norm` en `data/interim/sanson2018_recovered_official.csv`.
- Horlbeck2016, Gasperini2019 y Replogle2022 fueron inspeccionados en busca de sus variables originales.
- Todavía no se generó `efficacy_score`.
- Todavía no se generó `sgRNA_feature_matrix_phase1_with_efficacy.csv`.
- No se realizaron imputaciones, normalizaciones ni merges.
