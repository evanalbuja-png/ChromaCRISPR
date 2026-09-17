## Week 6 - Step 1: F2 ATAC feature extraction - FINAL STATUS

Date: 2026-07-26

Status: COMPLETED AND APPROVED

Summary:
- Extracted numerical ATAC-seq features from deepTools computeMatrix output.
- Parsed deepTools metadata JSON after removing '@' prefix.
- Extracted:
  - ATAC_mean
  - ATAC_max
  - ATAC_p90
  - ATAC_sum

Integration strategy:
- Created region_id:
  chr:start-end
- Association performed using inner join.
- Deduplication of region_id applied before merging.

Deduplication results:
- ATAC matrix:
  - Original regions: 155005
  - Duplicated region_id: 15274
  - Unique regions: 139731

- BED:
  - Original regions: 155454
  - Duplicated region_id: 15290
  - Unique regions: 140164

Final integration:
- Matched regions: 139731
- Coverage over unique BED: 99.691%

Decision:
The duplicated region_id entries are accepted as expected behavior for the current high-density sgRNA dataset.
Potential causes:
- multiple guides mapping to identical coordinates
- overlapping sgRNA regions

Further investigation of duplicate origins is deferred to optional QC.

Output:
data/interim/features/f2_atac_features.csv

## Week 6 - Step 2: Coordinate offset correction

Date: 2026-07-26

Issue:
Initial F2 integration showed low coverage (12.212%).

Diagnosis:
Coordinate comparison between sgRNA_unified_FINAL.csv and F2 features identified a systematic -1 bp offset.

Cause:
CSV coordinates were interpreted as 1-based while BED/deepTools coordinates are 0-based.

Correction:
Created region_id using:
start = coordinate - 1
end = start + 1

Result:
F2 integration recalculated after coordinate normalization.

## Week 6 - Step 3: F2 ATAC Feature QC

Date: 2026-07-26

Generated QC report:
logs/f2_features_qc_report.md

Evaluated:
- ATAC_mean distribution
- ATAC_max distribution
- ATAC_p90 distribution
- ATAC_sum distribution
- Zero and low signal percentages
- Dataset-stratified distributions

Input:
data/interim/features/sgRNA_with_f2.csv

## Step 4 - F3 H3K27ac feature extraction

Status: Completed

Input:
- data/interim/features/histones/H3K27ac_matrix.gz

Parameters:
- upstream/downstream: 500 bp
- bin size: 10 bp
- reference point: center

Features extracted:
- H3K27ac_mean
- H3K27ac_max
- H3K27ac_p90
- H3K27ac_sum

Output:
- data/interim/features/f3_h3k27ac_features.csv

QC:
- Matrix regions processed: <reported by script>
- Coverage: <reported by script>

## H3K27ac integration (F3)

Input:
- sgRNA_with_f2.csv
- f3_h3k27ac_features.csv

Method:
- region_id chr:start-end
- 1-based to 0-based correction
- pandas many_to_one merge

Output:
- data/interim/features/sgRNA_with_f2_f3.csv

Coverage:
- H3K27ac matched guides: <reported>
- H3K27ac coverage: <reported>


## F3 H3K4me3 extraction
- Matrix rows: 151170
- Unique regions: 136211
- Coverage: 97.180%
- Duplicate region_id removed: 14959


## H3K4me3 F3 integration
- Input dataset: data/interim/features/sgRNA_with_f2_f3.csv
- H3K4me3 regions: 136211
- Final guides: 155454
- Guides with H3K4me3: 151170
- Coverage: 97.244%


## F3 H3K27me3 extraction
- Matrix rows: 154371
- Unique regions: 139150
- Coverage: 99.277%
- Duplicate region_id removed: 15221


## F3 H3K27me3 integration
- Input: data/interim/features/sgRNA_with_f2_f3_full.csv
- Output: data/interim/features/sgRNA_with_f2_f3_repressive.csv
- Final guides: 155454
- H3K27me3 matched guides: 154371
- H3K27me3 coverage: 99.303%


# Semana 6 - Resumen ejecutivo final

## Objetivo
Integración de features epigenómicas F2 (ATAC-seq) y F3 (histone marks) para el dataset principal de sgRNA.

## Features generadas

### F2 - ATAC-seq K562
Archivo:
- `data/interim/features/f2_atac_features.csv`

Características:
- Ventana: ±500 bp
- Bin size: 10 bp
- Reference point: center
- Features:
  - ATAC_mean
  - ATAC_max
  - ATAC_p90
  - ATAC_sum

Cobertura:
- 139731 regiones únicas
- 99.69% respecto al BED único

Decisión:
- Se aceptó deduplicación por `region_id`.
- Los duplicados (~10%) fueron documentados como comportamiento esperado de datasets sgRNA de alta densidad.

---

## Integración F2

Archivo:
- `data/interim/features/sgRNA_with_f2.csv`

Cobertura final:
- 155005 / 155454 guías
- 99.711%

Decisión:
- Se corrigió discrepancia sistemática 1-based → 0-based.
- Se utilizó relación many-to-one entre guías y regiones epigenómicas.

---

# F3 - Histonas

## H3K27ac (marca activa)

Archivo:
- `data/interim/features/f3_h3k27ac_features.csv`

Cobertura:
- 140152 regiones únicas
- 99.991%

Integrado en:
- `data/interim/features/sgRNA_with_f2_f3.csv`

Cobertura:
- 155440 / 155454
- 99.991%

---

## H3K4me3 (promotores)

Archivo:
- `data/interim/features/f3_h3k4me3_features.csv`

Cobertura:
- 136211 regiones únicas

Integrado en:
- `data/interim/features/sgRNA_with_f2_f3_full.csv`

Cobertura:
- 151170 / 155454
- 97.244%

---

## H3K27me3 (represiva)

Archivo:
- `data/interim/features/f3_h3k27me3_features.csv`

Cobertura:
- 139150 regiones únicas

Integrado en:
- `data/interim/features/sgRNA_with_f2_f3_repressive.csv`

Cobertura:
- 154371 / 155454
- 99.303%

---

# Decisiones metodológicas

- Se mantuvo `region_id = chr:start-end` como identificador principal.
- Se aplicó corrección 1-based → 0-based para coordenadas provenientes del CSV de guías.
- Los merges fueron realizados con validación `many_to_one`.
- Las regiones sin match fueron conservadas con valores faltantes.
- Los duplicados de regiones fueron eliminados manteniendo la primera ocurrencia.
- Se evitó expandir el scope investigando duplicados de sgRNA.

---

# Commits principales Semana 6

- `4c79b6a`
  - Use region_id matching for F2 ATAC feature integration

- `c7aba15`
  - Fix 1-based to 0-based coordinate offset for F2 integration

- `d446a7c`
  - Integrate H3K27ac F3 features with F2 dataset

- `bc01dd8`
  - Integrate H3K4me3 F3 features with F2 F3 dataset

- `88141e5`
  - Extract F3 H3K27me3 histone repression features

- `cd33bf7`
  - Integrate H3K27me3 F3 repressive histone features

## Estado:
Semana 6 completada.
