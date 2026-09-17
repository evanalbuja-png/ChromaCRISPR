# Semana K — Consolidación matriz entrenamiento

**Fecha:** 2026-09-12  
**Matriz:** `data/processed/K562_training_matrix_v1.csv` — **19.446 × 120**  
**Manifest:** `data/processed/K562_feature_manifest_v1.tsv`

## Join
- Clave labels/F4/F5: (guide_sequence, chromosome, coordinate, gene_target) — 0 dups.
- F1/F2/F3 deduplicados por secuencia/locus antes del merge (evita fan-out multi-gen).
- Fan-out: **No**. Filas finales = 19.446.

## Target
- Primario: `EEP_percentile`
- Auxiliares (no features): raw_score, z_score_within_study

## NaN residuales (esperados)
| Columna | N | Nota |
|---|---|---|
| strand | 1128 | Gasperini |
| pam_* (3nt, is_NGG, one-hot, match_type) | 1013 | sin match genómico |
| E1 / compartment_A | 786 | flagged; no modelado base |
| hic_contact_* | 70 | sin contacto |
| genomic_dist / tss_* | 50–60 | TSS no mapeado |

## Flagged (no modelado base)
- F4 compartment E1 / compartment_A (sanity ATAC falló)

## Criterio de éxito
CUMPLIDO.

## Verificación final pre-cierre Bloque 2 (2026-09-12)

1. **PAM NaN en matriz consolidada:** pam_3nt / pam_is_NGG / 12 one-hot = 1013 NaN cada una; rate NGG skipna=0.999891; máscaras idénticas.

2. **Join F4/F5:** clave (guide_sequence, chromosome, coordinate, gene_target).  
   64 guide_sequence multi-gen (144 filas).  
   DISTINCT hic_contact_log: 56/64; DISTINCT log_dist_target_tss: 60/64.  
   Ejemplos ELK1/UXT, CYSTM1/UBE2D2 con valores distintos → join correcto, no colapso.

3. **Inventario 120 cols:** metadata 13 + accepted 101 + flagged 3 (E1, compartment_A/AB) + target EEP + label_aux 2.  
   Anti-leakage: ningún accepted con |ρ|>0.99 vs EEP (máx ~0.11). EEP no duplicado.

**Bloque 2 CERRADO.**
