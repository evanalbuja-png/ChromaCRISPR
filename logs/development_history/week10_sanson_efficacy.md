# Week 10 — Sanson2018 Efficacy Score Derivation (SO1 / Phase 2)

## 1. Archivos de raw reads utilizados
- `data/raw/crispr_datasets/Sanson2018/41467_2018_7901_MOESM6_ESM.xlsx`
- Hojas: `SetA raw reads`, `SetB raw reads`, `SetA sgRNA annotations`, `SetB sgRNA annotations`
- Diseño: Dolcetto CRISPRi (dCas9-KRAB), líneas celulares **HT29** y **A375** (no K562), pDNA (T0) + 3 réplicas biológicas por línea celular por set.

## 2. Método de conteo / normalización
- CPM normalization por columna: `(count + pseudocount) / (total_reads + pseudocount*n_guides) * 1e6`, pseudocount = 1.
- Score de eficacia por guía: `log2(mean(CPM réplicas post-selección) / CPM_pDNA)`, calculado independientemente para HT29 y A375.
- SetA y SetB combinados y deduplicados por secuencia exacta de sgRNA (57,050 + 57,011 → 114061 guías únicas).
- Guides filtradas: solo secuencias ACGT válidas con counts numéricos completos en las 7 columnas (pDNA + 3 réplicas x 2 líneas).

## 3. Cobertura
- **Total guías con score derivado: 114,061** (supera el criterio de éxito de ≥5,000).
- Overlap por secuencia exacta con matriz Phase 1 (Horlbeck, K562): **551 guías** de 18,318 en Phase 1 y 114,061 en Sanson2018 (~3.0% de Horlbeck, ~0.48% de Sanson2018).

## 4. Distribución del score
- lfc_HT29: mean=-0.205, std=0.930, rango=[-8.51, 3.47]
- lfc_A375: mean=-0.192, std=0.863, rango=[-5.48, 4.22]

## 5. Correlación con efficacy_score de Horlbeck (n=551 guías en overlap)
| Comparación | Spearman ρ | p-value | Pearson r | p-value |
|---|---|---|---|---|
| lfc_HT29 vs efficacy_score (K562) | -0.357 | 5.3e-18 | -0.355 | 8.1e-18 |
| lfc_A375 vs efficacy_score (K562) | -0.348 | 4.2e-17 | -0.373 | 1.3e-19 |

Signo negativo esperado y consistente: `lfc` mide depleción (más negativo = mayor represión = guía más efectiva), mientras que `efficacy_score` de Horlbeck está en la orientación opuesta.

## 6. Respuestas a las preguntas abiertas

**P1 — ¿Se puede calcular un LFC/depletion score limpio a nivel de sgRNA individual?**
Sí. Los raw reads tienen la estructura estándar de screen de dropout (pDNA + réplicas post-selección), permitiendo CPM + log2FC sin ambigüedad. Cobertura: 114,061 guías.

**P2 — ¿Cuántas hacen match con la matriz de features ya existente?**
551 de 114,061 por secuencia exacta (0.48%). El overlap bajo se debe a que Horlbeck (K562) y Sanson2018 (HT29/A375, librería Dolcetto) son librerías de guías diseñadas independientemente — coinciden solo donde la heurística de diseño produjo la misma secuencia de 20nt para el mismo gen.

**P3 — ¿Es usable para entrenamiento conjunto o solo validación externa?**
**Solo validación externa / expansión exploratoria, no entrenamiento conjunto directo con las features de cromatina actuales.** Razones:
1. Cell type mismatch: el score refleja contexto epigenómico de HT29/A375, pero las features F2/F3 de Phase 1 (ATAC-seq, ChIP-seq) están extraídas de K562. Usar el score de Sanson2018 con features de cromatina K562 introduciría una discordancia biológica (efficacy medido en una línea, contexto epigenómico de otra) — viola directamente el principio de trazabilidad y el diseño cell-type-conditioned del proyecto.
2. La correlación moderada (ρ≈0.35-0.37) con Horlbeck, aunque estadísticamente robusta, es consistente con lo esperado *entre* líneas celulares distintas — no sorprendentemente alta, no sorprendentemente baja. No hay evidencia de que sea "casi Horlbeck" transferible 1:1.
3. Uso correcto: (a) como dataset de **validación cruzada externa** una vez que se extraigan features de cromatina propias de HT29/A375 (si ENCODE tiene cobertura suficiente — pendiente de verificar, no evaluado esta semana); (b) como insumo directo para el experimento de transferencia cross-cell-type ya contemplado en SO5 del proposal original.

## 7. Decisión de cierre
El score de Sanson2018 **se acepta como dataset de validación externa potencial**, no se integra al training set de Phase 1/2 en este momento. Ítem de backlog: evaluar disponibilidad de ATAC-seq/ChIP-seq ENCODE para HT29 y A375 antes de intentar cualquier uso downstream de este score.

## 8. Artefactos generados
- `data/processed/sanson2018_lfc_scores.csv` (114,061 filas)
- `data/processed/sanson2018_horlbeck_overlap.csv` (551 filas)
- `scripts/week10_step3_compute_lfc.py`, `scripts/week10_step4_overlap_correlation.py`
