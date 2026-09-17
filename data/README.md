# Data sources (do not redistribute third-party raw files)

This repository **does not** ship raw CRISPR screen tables, ENCODE bigWigs, or Hi-C contact maps.
Those remain under their original licenses/terms. Below are the exact accessions and locations used in this study.

## CRISPR screens

| Dataset | Citation | Role | Public source |
|---------|----------|------|----------------|
| Horlbeck 2016 | eLife CRISPRi growth phenotype (γ) | Training | Paper supplementary (eLife 19760); local QC → 18,318 guides |
| Gasperini 2019 | Cell enhancer Perturb-seq | Training (high-confidence subset) | Paper Table S2; GEO-related materials as published |
| Sanson 2018 | Nat Commun Dolcetto CRISPRi | External validation only (HT29, A375) | Paper supplementary / Addgene library docs |
| Replogle 2022 | Cell genome-scale Perturb-seq | **Excluded** (no compact official guide-level efficacy table; full objects ~97–160 GB) | — |

## ENCODE / epigenomic tracks (K562, GRCh38)

| Signal | Experiment | File role (as used) |
|--------|------------|---------------------|
| ATAC-seq | ENCSR868FGK | fold-change / signal bigWig |
| DNase-seq | ENCSR000EKS | signal (e.g. ENCFF452XDU) |
| H3K27ac | ENCSR000AKP | fold enrichment |
| H3K4me3 | ENCSR000AKQ | fold enrichment |
| H3K27me3 | ENCSR000AKS | fold enrichment |
| H3K9me3 | ENCSR000DWD | fold enrichment |
| H3K4me1 | ENCSR000EWC | fold enrichment (corrected vs original proposal accession) |
| CTCF | ENCSR000AKO | fold enrichment (corrected vs original proposal accession) |
| ChromHMM 15-state | ENCSR365YNI | bed9 segmentation (K562) |

Portal: https://www.encodeproject.org/

## Hi-C

| Resource | Accession | Notes |
|----------|-----------|-------|
| K562 in situ Hi-C | 4DN 4DNFITUOMFUQ | 10 kb contacts; A/B exploratory only |

## Derived files shipped in this repo

See `data/processed/README.md`:

- `K562_training_matrix_v1.csv` — 19,446 guides × 101 features + labels  
- `K562_training_labels_harmonized_v1.csv` — EEP + within-study z-score  

These are **author-derived compilations**, released under CC-BY-4.0 (`LICENSE-DATA`).
