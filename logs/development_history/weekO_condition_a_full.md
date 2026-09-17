# Semana O — Condition A (101 features) + control de target

## Hallazgo principal (control)
El colapso ALL-set con EEP (ρ≈0.01) se debe **sobre todo al target within-gene**, no al domain shift de cromatina como causa dominante.

| Train target | Features | ALL \|ρ\| HT29 | ALL A375 | Overlap HT29 |
|--------------|----------|----------------|----------|--------------|
| EEP_percentile | 101 | 0.012 | 0.013 | 0.148 |
| **efficacy_score Horlbeck** | **101** | **0.148** | **0.171** | 0.106 |
| EEP seq-only | 36 | 0.101 | 0.105 | 0.067 |
| Phase1 ref | ~23 | ~0.22 | ~0.24 | — |

## Evidencia adicional
- Genes Sanson con ≥2 guías en training: **9.4 %** de filas (1.576 genes).
- EEP nunca entrenó comparación entre genes; Sanson ALL es casi todo between-gene.
- Feature shift ATAC/H3K27ac existe pero no explica el salto 0.01→0.15 al cambiar solo el target.

## Conclusión corregida
1. **Causa dominante del colapso EEP→Sanson ALL:** target relativo por gen (EEP).
2. **Domain shift cromatina:** factor secundario (0.17 aún < 0.22 Phase1).
3. **Dimensionalidad 101:** con target continuo, Condition A vuelve a rango ~0.15–0.17, no al techo 0.35–0.45 del proposal.
4. Para transferencia cross-library/cross-cell, el label debe ser **comparable between-gene** (efficacy continuo o LFC armonizado), no solo EEP.

## Condition B
Sigue bloqueada (ENCODE HT29/A375).

## RQ transfer (ajustada)
Más features no bastan; el **diseño del target** condiciona Condition A. EEP es correcto para LOCO K562 within-gene; para Condition A reportar también el control con efficacy_score continuo.
