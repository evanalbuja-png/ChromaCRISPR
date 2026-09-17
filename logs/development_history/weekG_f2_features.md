# Semana G — F2 Accessibility Features (cerrada, completa)

**Fecha:** 2026-09-10 / 2026-09-11  
**Base:** 19.446 guías de K562_training_labels_harmonized_v1.csv

## Interpretación de n=18
- 18 oficiales = 6 ventanas (±100, ±200, ±500, ±1kb, ±2kb, ±5kb) × (mean, max, p90).
- 3 extras (mismo estatus que guide_length): dnase_200bp_mean, atac_peak_summit_overlap, distance_to_nearest_atac_peak.

## Fuentes ENCODE
| Dato | Accession | Archivo |
|---|---|---|
| ATAC signal | ENCSR868FGK → ENCFF019IPA | ATAC-seq_K562.bigWig |
| ATAC peaks (IDR) | ENCSR868FGK → ENCFF993BAP | ATAC_K562_IDR_peaks.bed.gz |
| DNase signal | ENCSR000EKS → ENCFF452XDU | DNase-seq_K562.bigWig |

## Reproducibilidad
Spearman(ATAC_pm500_{mean,max,p90} vs FIXED) = **1.000000** (n=18.318).

## Resultado
`data/interim/F2_accessibility_features_v1.csv` — 19.446 filas, 18 oficiales + 3 extras, 0 NaNs.
- Summit overlap rate: 29.2 %
- Distance to nearest peak: mediana 92 bp

## Criterio de éxito
CUMPLIDO — F2 completo (18 + 3 extras).
