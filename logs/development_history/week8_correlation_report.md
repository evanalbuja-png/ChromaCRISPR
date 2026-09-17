# ChromaCRISPR Phase 1 — Week 8

## Paso 3 — Análisis de correlaciones

### Dataset

- Archivo: `data/processed/sgRNA_feature_matrix_phase1.csv`
- Filas: 155,454
- Features analizadas: 11
- Correlación Pearson: sí
- Correlación Spearman: sí
- Método: correlación por pares completos (`pairwise complete observations`)


## 1. Matriz de correlación Pearson

|                      |   ATAC_mean |   ATAC_max |   ATAC_p90 |   H3K27ac_mean |   H3K4me3_mean |   H3K27me3_mean |   gc_content |   mfe_rnafold |   guide_length |   nearest_tss_distance |   within_promoter_2kb |
|:---------------------|------------:|-----------:|-----------:|---------------:|---------------:|----------------:|-------------:|--------------:|---------------:|-----------------------:|----------------------:|
| ATAC_mean            |      1.0000 |     0.9301 |     0.9472 |         0.3522 |        -0.0019 |          0.0069 |       0.2325 |       -0.1149 |         0.0564 |                -0.0995 |                0.1656 |
| ATAC_max             |      0.9301 |     1.0000 |     0.9869 |         0.2649 |         0.0010 |          0.0072 |       0.1770 |       -0.0864 |         0.0515 |                -0.0662 |                0.1109 |
| ATAC_p90             |      0.9472 |     0.9869 |     1.0000 |         0.2772 |         0.0005 |          0.0067 |       0.1848 |       -0.0907 |         0.0543 |                -0.0737 |                0.1226 |
| H3K27ac_mean         |      0.3522 |     0.2649 |     0.2772 |         1.0000 |         0.0042 |          0.0225 |       0.0290 |       -0.0097 |         0.1038 |                -0.0935 |                0.1309 |
| H3K4me3_mean         |     -0.0019 |     0.0010 |     0.0005 |         0.0042 |         1.0000 |         -0.0023 |       0.0055 |       -0.0016 |        -0.0029 |                -0.0053 |               -0.0106 |
| H3K27me3_mean        |      0.0069 |     0.0072 |     0.0067 |         0.0225 |        -0.0023 |          1.0000 |       0.0048 |       -0.0034 |         0.0059 |                -0.0100 |               -0.0038 |
| gc_content           |      0.2325 |     0.1770 |     0.1848 |         0.0290 |         0.0055 |          0.0048 |       1.0000 |       -0.5727 |        -0.0039 |                -0.1542 |                0.2405 |
| mfe_rnafold          |     -0.1149 |    -0.0864 |    -0.0907 |        -0.0097 |        -0.0016 |         -0.0034 |      -0.5727 |        1.0000 |        -0.0642 |                 0.0775 |               -0.1235 |
| guide_length         |      0.0564 |     0.0515 |     0.0543 |         0.1038 |        -0.0029 |          0.0059 |      -0.0039 |       -0.0642 |         1.0000 |                -0.0205 |                0.0167 |
| nearest_tss_distance |     -0.0995 |    -0.0662 |    -0.0737 |        -0.0935 |        -0.0053 |         -0.0100 |      -0.1542 |        0.0775 |        -0.0205 |                 1.0000 |               -0.5538 |
| within_promoter_2kb  |      0.1656 |     0.1109 |     0.1226 |         0.1309 |        -0.0106 |         -0.0038 |       0.2405 |       -0.1235 |         0.0167 |                -0.5538 |                1.0000 |



## 2. Matriz de correlación Spearman

|                      |   ATAC_mean |   ATAC_max |   ATAC_p90 |   H3K27ac_mean |   H3K4me3_mean |   H3K27me3_mean |   gc_content |   mfe_rnafold |   guide_length |   nearest_tss_distance |   within_promoter_2kb |
|:---------------------|------------:|-----------:|-----------:|---------------:|---------------:|----------------:|-------------:|--------------:|---------------:|-----------------------:|----------------------:|
| ATAC_mean            |      1.0000 |     0.9402 |     0.9567 |         0.5605 |         0.0042 |          0.0045 |       0.2221 |       -0.1311 |         0.0951 |                -0.1289 |                0.1651 |
| ATAC_max             |      0.9402 |     1.0000 |     0.9885 |         0.4855 |         0.0091 |          0.0113 |       0.1763 |       -0.1047 |         0.0890 |                -0.0982 |                0.1135 |
| ATAC_p90             |      0.9567 |     0.9885 |     1.0000 |         0.5041 |         0.0077 |          0.0110 |       0.1852 |       -0.1099 |         0.0917 |                -0.1062 |                0.1237 |
| H3K27ac_mean         |      0.5605 |     0.4855 |     0.5041 |         1.0000 |         0.0058 |          0.0162 |       0.1193 |       -0.0730 |         0.1515 |                -0.0563 |                0.1490 |
| H3K4me3_mean         |      0.0042 |     0.0091 |     0.0077 |         0.0058 |         1.0000 |         -0.0112 |      -0.0042 |        0.0033 |        -0.0026 |                -0.0071 |               -0.0105 |
| H3K27me3_mean        |      0.0045 |     0.0113 |     0.0110 |         0.0162 |        -0.0112 |          1.0000 |       0.0262 |       -0.0178 |         0.0051 |                -0.0107 |                0.0010 |
| gc_content           |      0.2221 |     0.1763 |     0.1852 |         0.1193 |        -0.0042 |          0.0262 |       1.0000 |       -0.5988 |         0.0430 |                -0.1460 |                0.2328 |
| mfe_rnafold          |     -0.1311 |    -0.1047 |    -0.1099 |        -0.0730 |         0.0033 |         -0.0178 |      -0.5988 |        1.0000 |        -0.0744 |                 0.0894 |               -0.1374 |
| guide_length         |      0.0951 |     0.0890 |     0.0917 |         0.1515 |        -0.0026 |          0.0051 |       0.0430 |       -0.0744 |         1.0000 |                 0.0105 |                0.0382 |
| nearest_tss_distance |     -0.1289 |    -0.0982 |    -0.1062 |        -0.0563 |        -0.0071 |         -0.0107 |      -0.1460 |        0.0894 |         0.0105 |                 1.0000 |               -0.6499 |
| within_promoter_2kb  |      0.1651 |     0.1135 |     0.1237 |         0.1490 |        -0.0105 |          0.0010 |       0.2328 |       -0.1374 |         0.0382 |                -0.6499 |                1.0000 |



## 3. Correlaciones fuertes (|Pearson r| > 0.5)

| Feature 1 | Feature 2 | Pearson r | Spearman rho |
|---|---|---:|---:|

| `ATAC_max` | `ATAC_p90` | 0.9869 | 0.9885 |

| `ATAC_mean` | `ATAC_p90` | 0.9472 | 0.9567 |

| `ATAC_mean` | `ATAC_max` | 0.9301 | 0.9402 |

| `gc_content` | `mfe_rnafold` | -0.5727 | -0.5988 |

| `nearest_tss_distance` | `within_promoter_2kb` | -0.5538 | -0.6499 |


## 4. ATAC vs H3K27ac

| Feature 1 | Feature 2 | Pearson r | Spearman rho |
|---|---|---:|---:|

| `ATAC_mean` | `H3K27ac_mean` | 0.3522 | 0.5605 |

| `ATAC_max` | `H3K27ac_mean` | 0.2649 | 0.4855 |

| `ATAC_p90` | `H3K27ac_mean` | 0.2772 | 0.5041 |


## 5. H3K4me3 vs ATAC / H3K27ac

| Feature 1 | Feature 2 | Pearson r | Spearman rho |
|---|---|---:|---:|

| `H3K4me3_mean` | `ATAC_mean` | -0.0019 | 0.0042 |

| `H3K4me3_mean` | `ATAC_max` | 0.0010 | 0.0091 |

| `H3K4me3_mean` | `ATAC_p90` | 0.0005 | 0.0077 |

| `H3K4me3_mean` | `H3K27ac_mean` | 0.0042 | 0.0058 |


## 6. H3K27me3 vs features activas

| Feature 1 | Feature 2 | Pearson r | Spearman rho |
|---|---|---:|---:|

| `H3K27me3_mean` | `ATAC_mean` | 0.0069 | 0.0045 |

| `H3K27me3_mean` | `ATAC_max` | 0.0072 | 0.0113 |

| `H3K27me3_mean` | `ATAC_p90` | 0.0067 | 0.0110 |

| `H3K27me3_mean` | `H3K27ac_mean` | 0.0225 | 0.0162 |

| `H3K27me3_mean` | `H3K4me3_mean` | -0.0023 | -0.0112 |


## 7. Distancia al TSS vs cromatina

| Feature 1 | Feature 2 | Pearson r | Spearman rho |
|---|---|---:|---:|

| `nearest_tss_distance` | `ATAC_mean` | -0.0995 | -0.1289 |

| `nearest_tss_distance` | `ATAC_max` | -0.0662 | -0.0982 |

| `nearest_tss_distance` | `ATAC_p90` | -0.0737 | -0.1062 |

| `nearest_tss_distance` | `H3K27ac_mean` | -0.0935 | -0.0563 |

| `nearest_tss_distance` | `H3K4me3_mean` | -0.0053 | -0.0071 |

| `nearest_tss_distance` | `H3K27me3_mean` | -0.0100 | -0.0107 |


## 8. Interpretación automática básica

Las correlaciones con `|Pearson r| > 0.5` se consideran fuertes para este análisis descriptivo y se listan arriba.

La comparación Pearson/Spearman permite detectar si las relaciones observadas son principalmente lineales o si también existe asociación monotónica.


> Este análisis es exploratorio. Las correlaciones no implican causalidad y no se realizan eliminación de features, selección de variables ni transformación de datos en este paso.
