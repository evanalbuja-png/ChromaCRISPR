# Week 12 — Diagnóstico overlap Horlbeck↔Sanson (526 guías)

## Objetivo

Explicar por qué el modelo full (ρ≈−0.13) rinde peor que efficacy_score crudo (ρ≈−0.35) vs LFC en el overlap de secuencia idéntica.

## 1. Correlaciones en el overlap

- **lfc_HT29**: modelo vs LFC ρ=-0.1303 | efficacy_score vs LFC ρ=-0.3626
- **lfc_A375**: modelo vs LFC ρ=-0.1572 | efficacy_score vs LFC ρ=-0.3522
- **Modelo vs efficacy_score (target Horlbeck)**: ρ=0.2984

Si ρ(modelo, efficacy_score) es alto pero ρ(modelo, LFC) es bajo → el modelo predice bien el target de training pero ese target no se transfiere igual de bien que el score crudo (poco probable si son el mismo número).

Si ρ(modelo, efficacy_score) ya es bajo en el overlap → degradación del modelo *dentro* del subset (generalización / features).

## 2. Missingness / imputación

| Set | N | % con ≥1 missing | mean n_missing | % con ≥1 missing cromatina |
|-----|---|------------------|----------------|----------------------------|
| Sanson full | 100,942 | 3.22% | 0.167 | 3.22% |
| Overlap 526 | 526 | 2.85% | 0.114 | 2.85% |

### Missing rate por feature (overlap vs full)

| Feature | % missing full | % missing overlap |
|---------|----------------|-------------------|
| ATAC_mean | 0.43% | 0.00% |
| ATAC_max | 0.43% | 0.00% |
| ATAC_p90 | 0.43% | 0.00% |
| ATAC_sum | 0.43% | 0.00% |
| H3K27ac_mean | 0.01% | 0.00% |
| H3K27ac_max | 0.01% | 0.00% |
| H3K27ac_p90 | 0.01% | 0.00% |
| H3K27ac_sum | 0.01% | 0.00% |
| H3K4me3_mean | 3.02% | 2.85% |
| H3K4me3_max | 3.02% | 2.85% |
| H3K4me3_p90 | 3.02% | 2.85% |
| H3K4me3_sum | 3.02% | 2.85% |
| H3K27me3_mean | 0.70% | 0.00% |
| H3K27me3_max | 0.70% | 0.00% |
| H3K27me3_p90 | 0.70% | 0.00% |
| H3K27me3_sum | 0.70% | 0.00% |

## 3. Modelo vs efficacy_score en el overlap

- Spearman(pred, efficacy_score) = **0.2984**
- Residual (pred − score): mean=-0.2642, std=0.4321, MAE=0.4073

## 4. Top 10 |predicción − efficacy_score|

```
      guide_sequence  efficacy_score  pred_efficacy  abs_err  lfc_HT29  lfc_A375  ATAC_mean  H3K27ac_mean  H3K4me3_mean  nearest_tss_distance  _n_missing
GCAGCCACCCGGTCCCCCTC        2.490990       0.689658 1.801332 -2.619250 -2.085241   4.668720     16.809328      1.701818                  46.0           0
GTGTGTCGCATACCGCCCTC        2.354901       0.754407 1.600494 -2.678830 -3.681404   6.309431     27.657454      1.451222                  91.0           0
GGGCGGCGGATGGAGGTCAG        2.139713       0.603575 1.536138 -2.875585 -2.600597   5.894368     40.303056      1.930303                 389.0           0
GGCCCTGCGGTGTGACTCGC        1.889178       0.447541 1.441637 -3.044292 -3.202831   6.952569      0.553580      4.111520                9948.0           0
GCGGTACCCGGGCCCCGATG        2.134728       0.709462 1.425266 -2.222569 -3.044854   4.451493      2.931335      5.081200                  78.0           0
GCGCACCTCACTAGTCACGA        2.119035       0.698879 1.420156 -2.170683 -2.525328   7.452223     15.284645      1.640920                  35.0           0
GCGAGACCCCCTAGTAACAG        1.835133       0.451210 1.383923 -2.294019 -1.945733   3.188374     25.771054      1.941667                  69.0           0
GTCGTACTGACCGAGCGGGG        2.104114       0.737872 1.366241 -1.880549 -2.258967   7.417250      5.014606      2.442927                  60.0           0
GCGGCGCGCAAGGAAAGATC        2.026759       0.688297 1.338462 -3.245265 -2.089289   4.868210     41.781012      1.552716                  33.0           0
GCCTCGTGAACGAAGCGGTG        1.800767       0.462566 1.338201 -2.907981 -2.269081   4.626499     27.233278      0.865667                  29.0           0
```

## 5. Clasificación de causa (con evidencia)

**Causa más probable: generalización del ranking / optimización global**

Evidencia: Spearman(pred, efficacy_score) en overlap = 0.298 (razonable vs OOF ~0.30), pero pred vs LFC (-0.130) << score crudo vs LFC (-0.363). El score crudo conserva señal fenotípica punto-a-punto que el modelo, optimizado para el ranking global de ~18k guías, no preserva igual de bien en este subset de transferencia.

Nota: el score crudo *es* el target de esas 526 guías; cualquier compresión/regresión a la media del modelo reduce la correlación con un segundo fenotipo (LFC Sanson) si la relación efficacy→LFC es aproximadamente lineal y monotónica.
