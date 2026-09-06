# ChromaCRISPR Phase 1

## Estado del proyecto
Semana 1 - Configuración de infraestructura computacional

## Entorno
- WSL2 (Windows 11)
- Conda/Mamba
- Entorno dedicado: chromacrispr-phase1 (Python 3.10)

## Objetivo
Preparación de pipeline reproducible para:
- ATAC-seq / ChIP-seq (bigWig processing)
- Hi-C (cooler, hic2cool)
- Feature engineering genómico
- Machine Learning (XGBoost, SHAP, Optuna)

## Estructura
Pipeline reproducible basado en Snakemake.

Referencia: ChromaCRISPR Phase 1 Proposal

## Semana 3 - QC ENCODE + Hi-C inicializado (SO1 avance)

Durante la Semana 3 se realizó:

### QC de datos ENCODE K562

- Verificación de integridad de archivos BigWig.
- Evaluación de headers y cromosomas.
- Cálculo de estadísticas básicas de señal.
- Generación de reporte:


data/interim/qc_encode_report.md


### Preparación inicial Hi-C

Dataset:


4DNFITUOMFUQ.hic


Procesamiento:


hic2cool conversion
.cool generation
KR balancing


Resultado:


data/processed/hic_k562/4DNFITUOMFUQ.cool


Configuración:

- Resolución: 10 kb
- Formato: Cooler v3
- Cromosomas: 24

Esto completa el avance de SO1:
Dataset curation and harmonization.

## Week 4 - Unified sgRNA dataset QC (SO1)

Completed dataset curation for CRISPR guide datasets DS1-DS4:

- Horlbeck et al. 2016
- Sanson et al. 2018
- Gasperini et al. 2019
- Replogle et al. 2022

Pipeline implemented:

- Dataset loading
- Schema normalization
- Sequence QC
- Guide length validation
- Duplicate handling
- Dataset harmonization

Final dataset:
data/interim/sgRNA_unified.csv


Statistics:

- Total guides: 183,723
- Experiments: 9
- Genome builds: hg19/hg38
- QC status: passed

SO1 milestone:
Dataset curation completed.


# Semana 6 Status

## Epigenomic feature integration completed

Semana 6 completada.

Features integradas:

- F2:
  - ATAC-seq K562

- F3:
  - H3K27ac
  - H3K4me3
  - H3K27me3

Dataset final generado:

data/interim/features/sgRNA_with_f2_f3_repressive.csv
Coberturas:
- ATAC: 99.711%
- H3K27ac: 99.991%
- H3K4me3: 97.244%
- H3K27me3: 99.303%

La integración fue realizada mediante `region_id` con corrección de coordenadas 1-based → 0-based y merges many-to-one.

# Phase 1 — CLOSED

Phase 1 is formally completed.

### Final achievements

- Final multi-modal sgRNA feature matrix integrating:
  - F1: sequence features
  - F2: ATAC-seq chromatin accessibility
  - F3: histone-mark features
  - F5: genomic/TSS context
- Random Forest feature ablation completed.
- Final XGBoost multi-modal model evaluated with 5-fold GroupKFold.
- Internal benchmarking completed across position-only, sequence-only and multi-modal feature levels.
- SHAP analysis and TreeSHAP numerical validation completed.

### Best Phase 1 result

**XGBoost — F1+F2+F3+F5**

**Spearman ρ ≈ 0.322**

Final model performance:
- Spearman ρ: 0.32153
- R²: 0.10042
- RMSE: 0.41162

### Phase 2 backlog

The following remain deferred to Phase 2:

- F4 — Hi-C integration.
- Expansion and harmonization of additional efficacy datasets.
- External validation on an independent compatible efficacy dataset.
- Additional guide-scoring features/resources.

For the complete closure record, see:

`logs/phase1_final_summary.md`
