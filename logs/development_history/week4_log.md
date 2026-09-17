# ChromaCRISPR Phase 1
# Week 4 Log - sgRNA Dataset Curation

## Objective

Advance SO1:
Dataset curation and quality control of CRISPR screening datasets.

## Completed tasks

### Dataset integration

Integrated four CRISPR datasets:

- Horlbeck2016
- Sanson2018
- Gasperini2019
- Replogle2022


### Pipeline

Implemented:

- Dataset-specific parsers
- Unified schema
- Sequence validation
- Guide length QC
- Duplicate inspection
- Metadata preservation


## Final dataset

File:

data/interim/sgRNA_unified.csv


Statistics:

Initial loaded guides:
184086


Invalid guides removed:

363


Final guides:

183723


Retention:

99.8%


## Dataset composition

Gasperini2019:
16307

Horlbeck2016:
20809

Replogle2022:
32546

Sanson2018:
114061


## Problems encountered

1. Replogle dataset contained paired sgRNA libraries.
   Solution:
   Converted paired guides into unified guide representation.

2. Some repeated sequences mapped to different genes.
   Solution:
   Preserved biological duplicates across targets.

3. Genome build differences:
   hg19 and hg38 annotations retained.


## Pipeline changes

Added:

- parsers/replogle.py
- dataset manifest
- QC validation workflow


## Status

SO1 dataset curation milestone completed.
Ready for downstream feature engineering.
