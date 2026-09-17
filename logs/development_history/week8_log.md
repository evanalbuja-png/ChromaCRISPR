
## Paso 1 — Matriz final Phase 1

- Construcción de matriz final unificada.
- Filas: 155454
- Columnas: 34
- `sgrna_id` únicos: 155454
- `sgrna_id` duplicados: 0
- Archivo: `data/processed/sgRNA_feature_matrix_phase1.csv`

### Valores faltantes por columna

- `sgrna_id`: 0
- `dataset`: 0
- `experiment`: 119729
- `guide_sequence`: 0
- `gene_symbol`: 14553
- `gene_id`: 155454
- `chromosome`: 0
- `coordinate`: 0
- `strand`: 14553
- `region_id`: 0
- `gc_content`: 0
- `mfe_rnafold`: 0
- `guide_length`: 0
- `g_run_max`: 0
- `poly_t_count`: 0
- `ATAC_mean`: 449
- `ATAC_max`: 449
- `ATAC_p90`: 449
- `ATAC_sum`: 449
- `H3K27ac_mean`: 14
- `H3K27ac_max`: 14
- `H3K27ac_p90`: 14
- `H3K27ac_sum`: 14
- `H3K4me3_mean`: 4284
- `H3K4me3_max`: 4284
- `H3K4me3_p90`: 4284
- `H3K4me3_sum`: 4284
- `H3K27me3_mean`: 1083
- `H3K27me3_max`: 1083
- `H3K27me3_p90`: 1083
- `H3K27me3_sum`: 1083
- `nearest_tss_gene`: 0
- `nearest_tss_distance`: 0
- `within_promoter_2kb`: 0


## Paso 2 — EDA básico

- Se realizó EDA descriptivo sobre `data/processed/sgRNA_feature_matrix_phase1.csv`.
- Se generó resumen estadístico de features F1, F2, F3 y F5.
- Se evaluaron las distribuciones de `ATAC_mean`, `H3K27ac_mean`, `H3K4me3_mean`, `H3K27me3_mean`, `gc_content`, `mfe_rnafold` y `nearest_tss_distance`.
- Se calculó el porcentaje de guías dentro de `within_promoter_2kb`.
- Se calculó la distribución de guías por dataset.
- Se documentaron los valores faltantes.
- Reporte generado: `logs/week8_eda_report.md`.
- No se realizaron imputaciones, transformaciones ni eliminación de outliers.

## Paso 3 — Análisis de correlaciones

- Se calcularon matrices de correlación Pearson y Spearman para las features numéricas principales.
- Se evaluaron correlaciones fuertes mediante el umbral `|Pearson r| > 0.5`.
- Se analizaron específicamente las relaciones ATAC–H3K27ac, H3K4me3–cromatina activa, H3K27me3–features activas y distancia al TSS–cromatina.
- Reporte generado: `logs/week8_correlation_report.md`.
- No se realizó selección ni eliminación de features en este paso.

## Week 8 — CRISPRi activity score inspection

- Exact file: `data/raw/crispr_datasets/Horlbeck2016/elife-19760-supp1-v2.xlsx`
- Sheet: `CRISPRi`
- Activity-score rows: 18,380
- Unique activity gRNA sequences: 18,318
- Matrix unique guide sequences: 150,311
- Potential exact sequence matches: 18,318
- Potential sequence + gene matches: 18,380
- Activity score column: `CRISPRi activity score [Horlbeck et al., eLife 2016]`
- Score minimum: -0.832386701822
- Score maximum: 2.73859526503
- Score mean: 0.417475829109243
- Score median: 0.2771251783955
- Normalization assessment: El score no está acotado entre 0 y 1; parece ser un score raw/no normalizado a [0,1]. La confirmación definitiva requiere la documentación del dataset original.

## Integración controlada CRISPRi efficacy score
- Fuente: Horlbeck2016, hoja CRISPRi
- Llave: guide_sequence + gene_symbol
- Restricción: dataset == Horlbeck2016
- Filas matriz: 155,454
- Filas Horlbeck2016: 21,172
- efficacy_score no nulo: 18,386
- Cobertura: 86.84%
- Mean: 0.417438297131
- Median: 0.277001797694
- Min: -0.832386701822
- Max: 2.738595265030
- Filas no-Horlbeck con score: 0
- Output: data/processed/sgRNA_feature_matrix_phase1_with_efficacy.csv

# Semana 8 — Cierre formal

## Resumen ejecutivo

### 1. Matriz final

- Archivo: `data/processed/sgRNA_feature_matrix_phase1_with_efficacy.csv`
- Filas: 155,454
- Columnas: 35 (incluyendo `efficacy_score`)
- La matriz integra features F1 + F2 + F3 + F5 y el CRISPRi activity score cuando corresponde.

### 2. EDA principal

La matriz contiene las features de secuencia (F1), cromatina/epigenómica (F2/F3) y contexto genómico/TSS (F5). Se incorporó información de contexto promotor mediante `nearest_tss_distance` y `within_promoter_2kb`.

### 3. Correlaciones destacadas

Las relaciones entre las features y `efficacy_score` se evaluaron dentro del subset con efficacy disponible. Los modelos posteriores muestran que las señales de contexto genómico, estructura de RNA y accesibilidad de cromatina contienen información predictiva para la actividad CRISPRi.

### 4. Integración de efficacy_score

Se utilizó exclusivamente el dataset `Horlbeck2016` y la llave:

`guide_sequence + gene_symbol`

La integración controlada produjo:

- Filas totales de la matriz: 155,454
- Filas Horlbeck2016: 21,172
- Filas con `efficacy_score`: 18,386
- Cobertura dentro de Horlbeck2016: 86.84%
- Mean: 0.417438297131
- Median: 0.277001797694
- Min: -0.832386701822
- Max: 2.738595265030

Las coincidencias de secuencia presentes en otros datasets no recibieron `efficacy_score`.

### 5. Modelos baseline

Subset de entrenamiento: 18,386 guías con `efficacy_score`.

| Modelo | Spearman rho | Pearson R² | RMSE |
|---|---:|---:|---:|
| Ridge Regression | 0.261910 | 0.067347 | 0.412126 |
| LASSO | 0.256399 | 0.062422 | 0.413213 |
| Random Forest | 0.388259 | 0.124984 | 0.399188 |

Random Forest fue el mejor baseline según las tres métricas.

### 6. Random Forest — Top 10 Feature Importance

| Rank | Feature | Importance |
|---:|---|---:|
| 1 | nearest_tss_distance | 0.105402 |
| 2 | mfe_rnafold | 0.077376 |
| 3 | ATAC_max | 0.062573 |
| 4 | ATAC_p90 | 0.059295 |
| 5 | ATAC_sum | 0.056322 |
| 6 | ATAC_mean | 0.055955 |
| 7 | H3K27ac_max | 0.055353 |
| 8 | H3K27ac_p90 | 0.051447 |
| 9 | gc_content | 0.050333 |
| 10 | H3K4me3_p90 | 0.040188 |

`nearest_tss_distance` fue la feature individual más importante.

### 7. Decisiones metodológicas de Semana 8

1. Se utilizó `guide_sequence + gene_symbol` como llave de integración.
2. Se detectaron duplicados de esta llave en la matriz principal debido a guías compartidas entre datasets.
3. Para evitar contaminación entre datasets, `efficacy_score` se asignó únicamente a filas con `dataset == Horlbeck2016`.
4. Las filas sin correspondencia conservaron `NaN`.
5. El score de Horlbeck no se transformó ni normalizó antes del baseline.
6. El subset con efficacy disponible se utilizó exclusivamente para los modelos baseline.
7. Se compararon Ridge, LASSO y Random Forest con el mismo split `train/test` y `random_state=42`.
8. Random Forest se utilizó posteriormente para obtener el ranking de importancia de features.

## Estado de cierre

Semana 8 completada: matriz final, integración de efficacy, modelos baseline y análisis de feature importance completados.
