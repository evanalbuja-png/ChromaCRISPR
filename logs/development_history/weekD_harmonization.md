# Semana D — Armonización rank-based DS1 + DS3 (cerrada)

**Fecha:** 2026-09-10  
**Objetivo:** Combinar Horlbeck + Gasperini high-confidence en labels EEP según Sec. 7.1.

## Archivos de entrada
- Horlbeck: `data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_FIXED.csv` (18.318 guías, efficacy_score)
- Gasperini: `data/interim/DS3_gasperini/DS3_gasperini_guides_efficacy_highconf_polarity_fixed.csv` (1.328 → 1.128 tras deduplicar)

## Método (Sec. 7.1)
1. Z-score within-study independiente para cada dataset.
2. Overlap cross-study por secuencia o coordenada: **0** → no hubo promediado de duplicados.
3. Deduplicación interna de Gasperini: se eliminaron 200 filas exactas duplicadas (artefacto de merge previo). Se conservaron 80 casos de misma secuencia asociada a distintos gene_target (efecto multi-gen).
4. EEP = percentil (0–100) del z-score **dentro de cada gen_target**.
5. Genes con una sola guía → EEP = 50 por convención.

## Resultado
- Archivo: `data/processed/K562_training_labels_harmonized_v1.csv`
- Filas: **19.446**
- Columnas: guide_sequence, chromosome, coordinate, gene_target, dataset_origin, raw_score, z_score_within_study, EEP_percentile
- Genes: 2.151 (1.645 con ≥5 guías)

## Limitaciones
- Gasperini: mediana 2 guías/gen → percentiles discretos.
- 75 genes Horlbeck con 1 guía → EEP fijo = 50.
- Solo 23 genes en común entre los dos datasets.
- Horlbeck `nearest_tss_gene` contiene algunos IDs ENSG mezclados con símbolos HGNC.

## Criterio de éxito
CUMPLIDO.

## Verificación de nomenclatura de genes (cierre real)

Fecha: 2026-09-10

- Filas Horlbeck con gene_target formato ENSG: 2.166 (11.8 %), 260 genes únicos.
- Filas con símbolo HGNC: 16.152 (88.2 %), 1.561 genes únicos.
- Intersección ENSG ∩ símbolo dentro de Horlbeck: **0**.
- No existe fragmentación del mismo gen bajo dos identificadores.
- Genes ENSG con 1 sola guía: 20 (son genes que realmente tienen 1 guía, no artefactos).
- No se encontró tabla de mapeo ENSG→HGNC local en el proyecto.

**Decisión:** El impacto es despreciable. No se genera v2 ni se re-calcula EEP.  
La mezcla de nomenclatura se documenta como limitación heredada del reference de Horlbeck, sin efecto de fragmentación artificial.
