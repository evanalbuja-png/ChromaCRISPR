# ChromaCRISPR Phase 1 — Week 8

## Inspección de datasets originales de eficacia

## Horlbeck2016

Archivo: `data/raw/crispr_datasets/Horlbeck2016/elife-19760-supp1-v2.xlsx`

### Hojas

- `CRISPRi`
- `CRISPRa`

### Hoja `CRISPRi`

**Columnas:**

- `0`: `gene symbol`
- `1`: `chromosome`
- `2`: `PAM genomic coordinate [hg19]`
- `3`: `strand targeted`
- `4`: `sgRNA length (including PAM)`
- `5`: `sgRNA sequence`
- `6`: `CRISPRi activity score [Horlbeck et al., eLife 2016]`

**Posibles columnas de guía:**

- `sgRNA length (including PAM)`
- `sgRNA sequence`

**Posibles columnas de eficacia/score:**

- `CRISPRi activity score [Horlbeck et al., eLife 2016]`

**Primeras filas:**

```text
gene symbol chromosome  PAM genomic coordinate [hg19] strand targeted  sgRNA length (including PAM)        sgRNA sequence  CRISPRi activity score [Horlbeck et al., eLife 2016]
       AARS      chr16                       70323441               +                            24 GCGCTCTGATTGGACGGAGCG                                              0.019320
       AARS      chr16                       70323216               +                            24 GCCCCAGGATCAGGCCCCGCG                                              0.348892
       AARS      chr16                       70323296               +                            24 GGCCGCCCTCGGAGAGCTCTG                                              0.912409
```

### Hoja `CRISPRa`

**Columnas:**

- `0`: `gene symbol`
- `1`: `chromosome`
- `2`: `PAM genomic coordinate [hg19]`
- `3`: `strand targeted`
- `4`: `sgRNA length (including PAM)`
- `5`: `sgRNA sequence`
- `6`: `CRISPRa activity score`

**Posibles columnas de guía:**

- `sgRNA length (including PAM)`
- `sgRNA sequence`

**Posibles columnas de eficacia/score:**

- `CRISPRa activity score`

**Primeras filas:**

```text
gene symbol chromosome  PAM genomic coordinate [hg19] strand targeted  sgRNA length (including PAM)        sgRNA sequence  CRISPRa activity score
        AHR       chr7                       17337908               +                            24 GacaactggtagacaacCAAT               -0.005671
        AHR       chr7                       17337955               +                            24 GACGGTGGGTCAGCTAACTTG                0.072547
        AHR       chr7                       17338075               +                            23  GTAAGGACGCCCCCCCCCGC                0.181175
```

## Sanson2018

Archivo: `data/raw/crispr_datasets/Sanson2018/41467_2018_7901_MOESM6_ESM.xlsx`

### Hojas

- `SetA raw reads`
- `SetB raw reads`
- `SetA sgRNA annotations`
- `SetB sgRNA annotations`

### Hoja `SetA raw reads`

**Columnas:**

- `0`: `Supplementary Data 3. Dolcetto CRISPRi screening data`
- `1`: `Unnamed: 1`
- `2`: `Unnamed: 2`
- `3`: `Unnamed: 3`
- `4`: `Unnamed: 4`
- `5`: `Unnamed: 5`
- `6`: `Unnamed: 6`
- `7`: `Unnamed: 7`

**Posibles columnas de guía:**

- Ninguna detectada automáticamente.

**Posibles columnas de eficacia/score:**

- Ninguna detectada automáticamente.

**Primeras filas:**

```text
Supplementary Data 3. Dolcetto CRISPRi screening data Unnamed: 1 Unnamed: 2 Unnamed: 3 Unnamed: 4 Unnamed: 5 Unnamed: 6 Unnamed: 7
                                                  NaN       pDNA       HT29       HT29       HT29       A375       A375       A375
                                       sgRNA Sequence       pDNA       RepA       RepB       RepC       RepA       RepB       RepC
                                 AAAAAAAAAATACTGAGAGA        503        264        231        275        483        608        337
```

### Hoja `SetB raw reads`

**Columnas:**

- `0`: `Unnamed: 0`
- `1`: `pDNA`
- `2`: `HT29`
- `3`: `HT29.1`
- `4`: `HT29.2`
- `5`: `A375`
- `6`: `A375.1`
- `7`: `A375.2`

**Posibles columnas de guía:**

- Ninguna detectada automáticamente.

**Posibles columnas de eficacia/score:**

- Ninguna detectada automáticamente.

**Primeras filas:**

```text
          Unnamed: 0 pDNA HT29 HT29.1 HT29.2 A375 A375.1 A375.2
      sgRNA Sequence pDNA RepA   RepB   RepC RepA   RepB   RepC
AAAAAAAAAATTTCCTAGCG  257  275    201    187  394    137    155
AAAAAAAAAGCTGTGCGCAG  417  433    474    340  445    270    498
```

### Hoja `SetA sgRNA annotations`

**Columnas:**

- `0`: `sgRNA Sequence`
- `1`: `Annotated Gene Symbol`
- `2`: `Annotated Gene ID`

**Posibles columnas de guía:**

- `sgRNA Sequence`

**Posibles columnas de eficacia/score:**

- Ninguna detectada automáticamente.

**Primeras filas:**

```text
      sgRNA Sequence Annotated Gene Symbol  Annotated Gene ID
AAAAAAAAAATACTGAGAGA                 GATA3               2625
AAAAAAAAGAGGAGGGACGG                  ANKH              56172
AAAAAAAATTTCCTAGCGTG              SMARCAD1              56916
```

### Hoja `SetB sgRNA annotations`

**Columnas:**

- `0`: `sgRNA Sequence`
- `1`: `Annotated Gene Symbol`
- `2`: `Annotated Gene ID`

**Posibles columnas de guía:**

- `sgRNA Sequence`

**Posibles columnas de eficacia/score:**

- Ninguna detectada automáticamente.

**Primeras filas:**

```text
      sgRNA Sequence Annotated Gene Symbol  Annotated Gene ID
AAAAAAAAAATTTCCTAGCG              SMARCAD1              56916
AAAAAAAAAGCTGTGCGCAG          LOC100129112          100129112
AAAAAAAACTGTCCCGCAAC              MAP1LC3C             440738
```

## Gasperini2019_S1

Archivo: `data/raw/crispr_datasets/Gasperini2019/NIHMS1038673-supplement-TableS1.xlsx`

### Hojas

- `A_Pilot_gRNA_library`
- `B_Pilot_145_enhancergene_pairs`

### Hoja `A_Pilot_gRNA_library`

**Columnas:**

- `0`: `Spacer`
- `1`: `Target_Site`
- `2`: `chr.candidate_enhancer`
- `3`: `start.candidate_enhancer`
- `4`: `stop.candidate_enhancer`
- `5`: `category`

**Posibles columnas de guía:**

- `Spacer`

**Posibles columnas de eficacia/score:**

- Ninguna detectada automáticamente.

**Primeras filas:**

```text
              Spacer Target_Site chr.candidate_enhancer  start.candidate_enhancer  stop.candidate_enhancer           category
CAAGCTGGTTAAAAACCCCG    chr1.576                   chr1                   2472004                  2474362 candidate_enhancer
GCAAAGTAGGTCTTTCCCAG    chr1.576                   chr1                   2472004                  2474362 candidate_enhancer
GGCCCGCCAGGAAAGCTGCA   chr1.1432                   chr1                   8257629                  8258384 candidate_enhancer
```

### Hoja `B_Pilot_145_enhancergene_pairs`

**Columnas:**

- `0`: `Target_Site`
- `1`: `ENSG`
- `2`: `target_gene_short`
- `3`: `Diff_expression_test_raw_pval`
- `4`: `Diff_expression_test_fold_change`
- `5`: `chr.candidate_enhancer`
- `6`: `start.candidate_enhancer`
- `7`: `stop.candidate_enhancer`

**Posibles columnas de guía:**

- Ninguna detectada automáticamente.

**Posibles columnas de eficacia/score:**

- `Diff_expression_test_fold_change`

**Primeras filas:**

```text
Target_Site            ENSG target_gene_short  Diff_expression_test_raw_pval  Diff_expression_test_fold_change chr.candidate_enhancer  start.candidate_enhancer  stop.candidate_enhancer
  chr4.1626 ENSG00000109255               NMU                            0.0                          0.206331                   chr4                  56595855                 56596581
  chr4.1627 ENSG00000109255               NMU                            0.0                          0.217810                   chr4                  56596632                 56597095
   chrX.232 ENSG00000205542            TMSB4X                            0.0                          0.369964                   chrX                  12973885                 12974748
```

## Gasperini2019_S2

Archivo: `data/raw/crispr_datasets/Gasperini2019/NIHMS1038673-supplement-TableS2.xlsx`

### Hojas

- `Key`
- `S2A_AtScale_library_gRNA.cs`
- `B_AtScale_664_enhancergenepairs`

### Hoja `Key`

**Columnas:**

- `0`: `Tab S2A | Annotated table of the gRNAs and candidate enhancers for the at-scale experiment, as well as all controls`

**Posibles columnas de guía:**

- `Tab S2A | Annotated table of the gRNAs and candidate enhancers for the at-scale experiment, as well as all controls`

**Posibles columnas de eficacia/score:**

- Ninguna detectada automáticamente.

**Primeras filas:**

```text
                                      Tab S2A | Annotated table of the gRNAs and candidate enhancers for the at-scale experiment, as well as all controls
                                                                            Spacer | sequence of the 19 (for TSS controls) or 20 bp (for all else) spacer
Target_Site | ID for targeted locus. Candidate enhancers’ IDs are derived from the original DHS peaks used to define theenhancersite (ENCODE ENCFF001UWQ)
                                                                                    chr.candidate_enhancer  | chromosomal location of candidate enhancers
```

### Hoja `S2A_AtScale_library_gRNA.cs`

**Columnas:**

- `0`: `Spacer`
- `1`: `Target_Site`
- `2`: `chr.candidate_enhancer`
- `3`: `start.candidate_enhancer`
- `4`: `stop.candidate_enhancer`
- `5`: `Category`

**Posibles columnas de guía:**

- `Spacer`

**Posibles columnas de eficacia/score:**

- Ninguna detectada automáticamente.

**Primeras filas:**

```text
              Spacer Target_Site  chr.candidate_enhancer  start.candidate_enhancer  stop.candidate_enhancer Category
AATGAGGAGCAAACGAAAAT     control                     NaN                       NaN                      NaN      NTC
ACGAAATGTTTCATGACCAA     control                     NaN                       NaN                      NaN      NTC
ATAGATTTACGTTACTCTCT     control                     NaN                       NaN                      NaN      NTC
```

### Hoja `B_AtScale_664_enhancergenepairs`

**Columnas:**

- `0`: `Target_Site`
- `1`: `ENSG`
- `2`: `target_gene_short`
- `3`: `Diff_expression_test_raw_pval`
- `4`: `Diff_expression_test_fold_change`
- `5`: `Diff_expression_test_Empirical_pval`
- `6`: `Diff_expression_test_Empirical_adjusted_pval`
- `7`: `high_confidence_subset`
- `8`: `chr.candidate_enhancer`
- `9`: `start.candidate_enhancer`
- `10`: `stop.candidate_enhancer`

**Posibles columnas de guía:**

- Ninguna detectada automáticamente.

**Posibles columnas de eficacia/score:**

- `Diff_expression_test_fold_change`

**Primeras filas:**

```text
Target_Site            ENSG target_gene_short  Diff_expression_test_raw_pval  Diff_expression_test_fold_change  Diff_expression_test_Empirical_pval  Diff_expression_test_Empirical_adjusted_pval  high_confidence_subset chr.candidate_enhancer  start.candidate_enhancer  stop.candidate_enhancer
  chr2.2482 ENSG00000115977              AAK1                   1.451572e-03                          0.756542                             0.002719                                      0.098652                    True                   chr2                  69056234                 69056865
  chrX.2695 ENSG00000101986             ABCD1                   7.351840e-04                          0.669369                             0.001825                                      0.073014                    True                   chrX                 153250743                153251468
 chr10.2252 ENSG00000138316          ADAMTS14                   2.140000e-09                          0.447355                             0.000449                                      0.032498                   False                  chr10                  72426863                 72427518
```

## Replogle2022

Archivo: `data/raw/crispr_datasets/Replogle2022/NIHMS1812939-supplement-11.xlsx`

### Hojas

- `TabA_K562_day8_library`
- `TabB_K562_day6_library`
- `TabC_RPE1_day7_library`

### Hoja `TabA_K562_day8_library`

**Columnas:**

- `0`: `unique sgRNA pair ID`
- `1`: `gene`
- `2`: `transcript`
- `3`: `ensembl gene id`
- `4`: `sgID_A`
- `5`: `targeting sequence A`
- `6`: `sgID_B`
- `7`: `targeting sequence B`
- `8`: `duplicated guide pair?`
- `9`: `either guide duplicated?`

**Posibles columnas de guía:**

- `unique sgRNA pair ID`
- `targeting sequence A`
- `targeting sequence B`
- `duplicated guide pair?`
- `either guide duplicated?`

**Posibles columnas de eficacia/score:**

- Ninguna detectada automáticamente.

**Primeras filas:**

```text
       unique sgRNA pair ID gene transcript ensembl gene id                  sgID_A targeting sequence A                  sgID_B targeting sequence B  duplicated guide pair?  either guide duplicated?
  0_A1BG_P1_ENSG00000121410 A1BG         P1 ENSG00000121410   A1BG_+_58858964.23-P1 GCTCCGGGCGACGTGGAGTG   A1BG_-_58858788.23-P1 GGGGCACCCAGGAGCGGTAG                   False                     False
  1_A1BG_P2_ENSG00000121410 A1BG         P2 ENSG00000121410   A1BG_-_58864840.23-P2 GCCGGTGCAGTGAGTGTCTG   A1BG_-_58864822.23-P2 GATGATGGTCGCGCTCACTC                   False                     False
2_AAAS_P1P2_ENSG00000094914 AAAS       P1P2 ENSG00000094914 AAAS_-_53715438.23-P1P2 GAGGACGAGTACGCGGTCCC AAAS_+_53715355.23-P1P2 GCCTCGCCGTTTGTCCCTTG                   False                     False
```

### Hoja `TabB_K562_day6_library`

**Columnas:**

- `0`: `unique sgRNA pair ID`
- `1`: `gene`
- `2`: `transcript`
- `3`: `ensembl gene id`
- `4`: `sgID_A`
- `5`: `targeting sequence A`
- `6`: `sgID_B`
- `7`: `targeting sequence B`
- `8`: `duplicated guide pair?`
- `9`: `either guide duplicated?`

**Posibles columnas de guía:**

- `unique sgRNA pair ID`
- `targeting sequence A`
- `targeting sequence B`
- `duplicated guide pair?`
- `either guide duplicated?`

**Posibles columnas de eficacia/score:**

- Ninguna detectada automáticamente.

**Primeras filas:**

```text
         unique sgRNA pair ID  gene transcript ensembl gene id                   sgID_A targeting sequence A                   sgID_B targeting sequence B  duplicated guide pair?  either guide duplicated?
  2_AAAS_P1P2_ENSG00000094914  AAAS       P1P2 ENSG00000094914  AAAS_-_53715438.23-P1P2 GAGGACGAGTACGCGGTCCC  AAAS_+_53715355.23-P1P2 GCCTCGCCGTTTGTCCCTTG                   False                     False
  8_AAMP_P1P2_ENSG00000127837  AAMP       P1P2 ENSG00000127837 AAMP_+_219134851.23-P1P2 GGTCGCGCAGAGCTGACTCT AAMP_+_219134841.23-P1P2 GGCTGACTCTGGGAGGCGTT                   False                     False
10_AARS2_P1P2_ENSG00000124608 AARS2       P1P2 ENSG00000124608 AARS2_+_44281027.23-P1P2 GAGTGGCAGCTGCAGCCCGG AARS2_+_44281044.23-P1P2 GGCTACGATGGCAGCGTCAG                   False                     False
```

### Hoja `TabC_RPE1_day7_library`

**Columnas:**

- `0`: `unique sgRNA pair ID`
- `1`: `gene`
- `2`: `transcript`
- `3`: `ensembl gene id`
- `4`: `sgID_A`
- `5`: `targeting sequence A`
- `6`: `sgID_B`
- `7`: `targeting sequence B`
- `8`: `duplicated guide pair?`
- `9`: `either guide duplicated?`

**Posibles columnas de guía:**

- `unique sgRNA pair ID`
- `targeting sequence A`
- `targeting sequence B`
- `duplicated guide pair?`
- `either guide duplicated?`

**Posibles columnas de eficacia/score:**

- Ninguna detectada automáticamente.

**Primeras filas:**

```text
       unique sgRNA pair ID gene transcript ensembl gene id                   sgID_A targeting sequence A                   sgID_B targeting sequence B  duplicated guide pair?  either guide duplicated?
2_AAAS_P1P2_ENSG00000094914 AAAS       P1P2 ENSG00000094914  AAAS_-_53715438.23-P1P2 GAGGACGAGTACGCGGTCCC  AAAS_+_53715355.23-P1P2 GCCTCGCCGTTTGTCCCTTG                   False                     False
8_AAMP_P1P2_ENSG00000127837 AAMP       P1P2 ENSG00000127837 AAMP_+_219134851.23-P1P2 GGTCGCGCAGAGCTGACTCT AAMP_+_219134841.23-P1P2 GGCTGACTCTGGGAGGCGTT                   False                     False
9_AAR2_P1P2_ENSG00000131043 AAR2       P1P2 ENSG00000131043  AAR2_-_34824434.23-P1P2 GTGGGGCGAGGCGGTGAGTG  AAR2_+_34824488.23-P1P2 GGACTCTGAGCCGAGAAGAG                   False                     False
```
