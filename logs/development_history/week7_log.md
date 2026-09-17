## Week 7 - F1 Sequence Features Extraction

Date:
2026-08-02

## Objective

Extract and integrate sequence-level features (F1) from the main sgRNA dataset.

Input dataset:

data/interim/features/sgRNA_with_f2_f3_repressive.csv


Rows processed:

155,454 sgRNAs


## Generated features

The following sequence features were extracted:

- GC content
- guide length
- RNAfold minimum free energy (MFE)
- maximum G homopolymer length
- poly-T motif count


## RNAfold calculation

ViennaRNA RNAfold 2.7.2 was used.

Calculation mode:

Batch processing of unique guide sequences.


Unique sequences processed:

150,311


Results obtained:

150,311


Final coverage:

155,454 / 155,454 guides


## Guide length validation

The original dataset_length column was evaluated.

A total of 21,172 mismatches were detected.

All mismatches corresponded to Horlbeck2016.

Observed pattern:

dataset_length = sequence_length + 3


This was interpreted as guide + PAM length.

Decision:

guide_length was calculated exclusively from guide_sequence length.


## Final output

Generated file:

data/interim/features/f1_sequence_features.csv


Columns:

- guide_sequence
- region_id
- gc_content
- mfe_rnafold
- guide_length
- g_run_max
- poly_t_count


Validation:

Rows preserved:
155,454 / 155,454

Missing values:
0


Status:

F1 Sequence Features completed.

## F5.2 - TSS Genomic Context Features Completed

Date: 2026-08-02

### Objective
Integration of genomic context features (F5) using GENCODE v49 TSS annotations.

### Reference annotation

- Genome build: hg38 / GRCh38
- Annotation source: GENCODE v49
- Input file:
  - data/reference/hg38/gencode_v49_TSS.tsv
- Total TSS annotations:
  - 78,691

### Methodology

TSS proximity features were calculated using:

- bedtools v2.31.1
- bedtools closest -d

sgRNA coordinates were converted to BED format using:

- chromosome
- start_0based
- end_0based

A unique identifier (`sgrna_id`) was introduced to preserve the original experimental records because `region_id` is not unique across the integrated CRISPR datasets.

### Generated features

Output:

data/interim/features/f5_tss_features.csv

Features:

- sgrna_id
- region_id
- nearest_tss_gene
- nearest_tss_strand
- nearest_tss_distance
- within_promoter_2kb

Promoter definition:

- ±2 kb around nearest TSS

### Validation

Input sgRNA records:

155,454

Output records:

155,454

Unique sgRNA IDs:

155,454

Missing records:

0

TSS ties detected:

316

Tie resolution:

Multiple equally distant TSS assignments were resolved deterministically by keeping the first bedtools closest assignment.

### Status

F5.2 TSS context extraction completed successfully.
Ready for F1 + F5 feature integration.

## F1 + F5 Feature Integration Completed

Date: 2026-08-02

### Objective

Integration of sequence features (F1) and genomic context features (F5) into the main CRISPR guide feature matrix.

### Input datasets

Base dataset:

data/interim/features/sgRNA_with_f2_f3_repressive.csv

Features integrated:

F1:
- data/interim/features/f1_sequence_features.csv

F5 genomic context:
- data/interim/features/f5_genomic_context_features.csv

F5 TSS context:
- data/interim/features/f5_tss_features.csv

### Output

Generated:

data/interim/features/sgRNA_with_f1_f5.csv

### Integration method

All feature tables were integrated using a stable unique identifier:

- sgrna_id

The original row ordering was preserved for F1 and F5.1 features after validation.

### Validation

Expected rows:
155454

Final rows:
155454

Duplicate sgrna_id:
0

Missing values in critical integrated features:
0

Integrated feature groups:

- F1 sequence features
- F2 ATAC accessibility
- F3 histone marks
- F5 genomic context

### Status

F1 + F5 integration completed successfully.

SO2 feature matrix expansion in progress.
