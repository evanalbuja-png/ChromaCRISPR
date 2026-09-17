# Semana H — F3 Histone + ChromHMM (cerrada)

**Fecha:** 2026-09-11  
**Interpretación:** Opción A — 31 oficiales (lectura literal Sec. 6.3).  
**Discrepancia:** proposal declara n=42; texto suma 31. Documentada, no forzada.

## Oficiales implementadas (31)
| Grupo | Features | N |
|---|---|---|
| 5 marks × 3 windows × mean | H3K27ac, H3K4me3, H3K4me1, H3K27me3, H3K9me3 × (±500, ±1kb, ±2kb) | 15 |
| CTCF ±1 kb mean | CTCF_pm1000_mean | 1 |
| ChromHMM 15-state one-hot | chromhmm_{Quies,TxWk,Enh1,...} | 15 |

## Accessions ENCODE / Roadmap
| Track | Experimento | Archivo |
|---|---|---|
| H3K27ac | ENCSR000AKP → ENCFF502GAA | H3K27ac_K562.bigWig |
| H3K4me3 | ENCSR000AKQ → ENCFF000BXD | H3K4me3_K562.bigWig |
| H3K27me3 | ENCSR000AKS → ENCFF900IMQ | H3K27me3_K562.bigWig |
| H3K9me3 | ENCSR000DWD → ENCFF949BXN | H3K9me3_K562.bigWig |
| H3K4me1 | **ENCSR000EWC** → ENCFF761XBZ | H3K4me1_K562.bigWig (proposal ENCSR000DWA era Control, incorrecto) |
| CTCF | **ENCSR000AKO** → ENCFF979PWH | CTCF_K562.bigWig (proposal ENCSR000AKR era H3K36me3, incorrecto) |
| ChromHMM 15-state | ENCSR365YNI → ENCFF106BGJ | ChromHMM_K562_15state.bed.gz |

## Reproducibilidad
- H3K27ac_pm500 vs FIXED: Spearman **1.000000**
- H3K4me3_pm500 vs FIXED: **0.999648**
- H3K27me3_pm500 vs FIXED: **0.999762**

## NaNs
74 guías (chr17/18/20, coords altas) sin señal en H3K4me3/H3K27me3 bigWigs → imputadas a 0.

## Archivo
`data/interim/F3_histone_chromhmm_features_v1.csv` — 19.446 × 31 oficiales (+ ids + chromhmm_state).

## Criterio de éxito
CUMPLIDO.

## Imputación (corregida, Sec. 7.3)
- 74 guías con NaN en H3K4me3 y H3K27me3 (3 ventanas cada uno; mismo set de guías; chr17/18/20 coords altas).
- Patrón: las 3 ventanas faltan siempre juntas por mark → 1 flag binario por mark (no 6).
- Método: **mediana del training set** (no 0) + columnas `H3K4me3_missing`, `H3K27me3_missing`.
- Flags son extras de pipeline (Sec. 7.3), no cuentan dentro de las 31 oficiales de F3.

## Verificación de accessions corregidos (fuente: metadata ENCODE JSON)
Consultado vía `https://www.encodeproject.org/experiments/{ACC}/?format=json`:
| Accession del proposal | Campo description / target en ENCODE | Conclusión |
|---|---|---|
| ENCSR000DWA | "Control ChIP-seq on human K562" | **No es H3K4me1** |
| ENCSR000AKR | "H3K36me3 ChIP-seq on human K562" | **No es CTCF** |
| ENCSR000EWC | "H3K4me1 ChIP-seq on human K562" | Usado (ENCFF761XBZ fold-change) |
| ENCSR000AKO | "CTCF ChIP-seq on human K562" | Usado (ENCFF979PWH fold-change) |

Verificación no se basó en nombres de archivo locales, sino en el campo de metadata del experimento en el portal ENCODE.
