# Week 11 — ENCODE epigenomic coverage: HT-29 and A375

## Scope

Verification of ENCODE coverage for HT-29 and A375 before deciding whether the
K562-derived F2/F3 chromatin feature extraction can be extended to external
cell lines.

Target chromatin features:

- Accessibility: ATAC-seq or DNase-seq
- H3K27ac
- H3K4me3
- H3K4me1
- H3K27me3
- H3K9me3

No data were downloaded during this inventory. Only ENCODE REST API metadata
were queried.

---

## HT-29

ENCODE canonical biosample name: `HT-29`

Direct experiment inventory:

| Accession | Status | Assay | Relevant for F2/F3 |
|---|---|---|---|
| ENCSR814KRX | released | DNase-seq | Yes — accessibility substitute |
| ENCSR807KQN | released | small RNA-seq | No |
| ENCSR295GUQ | released | RAMPAGE | No |
| ENCSR971GPJ | released | total RNA-seq | No |

### Accessibility

HT-29 has DNase-seq experiment `ENCSR814KRX`.

The experiment contains four biological replicates. Released GRCh38 bigWig
files with `read-depth normalized signal` were identified for the biological
replicates.

Representative released GRCh38 bigWig files:

| Biological replicate | BigWig |
|---|---|
| 1 | ENCFF528SGI |
| 2 | ENCFF441MDS |
| 3 | ENCFF846SFA |
| 4 | ENCFF074JTF |

These files are directly associated with `ENCSR814KRX`.

Archived hg19 bigWig files were not considered for the GRCh38 feature
extraction route.

### Histone marks

No direct HT-29 experiment was identified for:

- H3K27ac
- H3K4me3
- H3K4me1
- H3K27me3
- H3K9me3

Therefore HT-29 has accessibility coverage but no verified coverage for the
five histone marks required to reproduce the K562 F2/F3 feature set.

---

## A375

ENCODE canonical biosample name: `A375`

Direct experiment inventory:

| Accession | Status | Assay | Relevant for F2/F3 |
|---|---|---|---|
| ENCSR376XXO | released | small RNA-seq | No |
| ENCSR080XAB | released | RAMPAGE | No |
| ENCSR504VXC | released | total RNA-seq | No |

### Accessibility

No ATAC-seq or DNase-seq experiment was identified in the direct A375
experiment inventory.

### Histone marks

No direct A375 experiment was identified for:

- H3K27ac
- H3K4me3
- H3K4me1
- H3K27me3
- H3K9me3

Therefore A375 currently has no verified ENCODE chromatin coverage suitable
for reproducing the K562 F2/F3 feature set.

---

## K562 reference

The original K562 feature set uses:

| Feature | ENCODE accession |
|---|---|
| ATAC-seq | ENCSR868FGK |
| H3K27ac | ENCSR000AKP |
| H3K4me3 | ENCSR000AKQ |
| H3K4me1 | ENCSR000DWA |
| H3K27me3 | ENCSR000AKS |
| H3K9me3 | ENCSR000DWD |
| CTCF | ENCSR000AKR |
| DNase-seq | ENCSR000EKS |

The external validation requirement is therefore not met by a simple
one-to-one chromatin feature mapping in either HT-29 or A375.

---

## Coverage summary

| Cell line | Accessibility | H3K27ac | H3K4me3 | H3K4me1 | H3K27me3 | H3K9me3 |
|---|---|---|---|---|---|---|
| K562 | Yes | Yes | Yes | Yes | Yes | Yes |
| HT-29 | Yes (DNase) | No | No | No | No | No |
| A375 | No | No | No | No | No | No |

---

## Interpretation

HT-29 provides a potentially usable accessibility resource through DNase-seq,
including released GRCh38 normalized bigWig signal across four biological
replicates. However, the five histone marks required for the K562-derived F2/F3
feature set are not available in the direct HT-29 ENCODE inventory.

A375 has neither ATAC-seq nor DNase-seq nor the five target histone marks in
its direct ENCODE experiment inventory.

Consequently, neither HT-29 nor A375 currently provides sufficient matched
epigenomic coverage to reproduce the complete K562 F2/F3 feature extraction
under a uniform pipeline.

## Decision status

**F2/F3 external validation with HT-29 and A375: NOT VIABLE as a complete
matched-feature replication.**

HT-29 may remain useful for a future accessibility-only analysis, but that
would constitute a different validation question and should not be presented
as replication of the complete K562 F2/F3 feature set.

No pipeline changes were made.

No ENCODE data files were downloaded.

Final Week 11 route decision remains subject to the project protocol:
document the cell-type/epigenomic mismatch and evaluate whether to return to
Option A or close the Sanson2018 external-validation path.
