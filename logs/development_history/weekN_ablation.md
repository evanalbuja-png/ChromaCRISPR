# Semana N — Ablation study (Sec. 10.3)

**Modelo:** XGBoost (params Semana L, fijos). **Eval:** LOCO 23 folds.  
**Permutación:** 10,000 over paired fold Δρ.

## Grupos reales (manifest)
F1=36, F2=21, F3=33, F4=2, F5=9 (total 101).

## Secuencial
| Step | ρ_LOCO | Δ | p |
|------|--------|---|---|
| F1 | 0.329±0.026 | — | — |
| +F2 | 0.360±0.031 | +0.032 | ≈0 |
| +F3 | 0.360±0.033 | −0.001 | 0.70 |
| +F4 | 0.358±0.031 | −0.002 | 0.18 |
| +F5 FULL | 0.376±0.033 | +0.018 | ≈0 |

## Inverso (contribución marginal)
| Removed | ρ without | Δ vs FULL | p |
|---------|-----------|-----------|---|
| F1 | 0.223 | 0.153 | ≈0 |
| F2 | 0.355 | 0.021 | ≈0 |
| F3 | 0.376 | ~0 | 0.89 |
| F4 | 0.377 | ~0 | 0.60 |
| F5 | 0.358 | 0.018 | 0.0001 |

## RQ1 (¿cromatina mejora?)
**Sí.** Accesibilidad (F2) incrementa ρ de forma significativa sobre secuencia sola; F5 añade otra ganancia significativa. FULL 0.376 > F1 0.329.

## RQ2 (¿qué modalidad más?)
Orden de contribución marginal: **F1 >> F2 ≈ F5 >> F3 ≈ F4 (~0)**.  
F3/F4 no aportan en presencia de F1+F2+F5 (redundancia ATAC–histona y Hi-C poco informativo en set mayormente <5 kb TSS).

## Nota
FULL LOCO 0.376 ≈ Semana M XGB 0.377 (consistente).

## Nota F4
Δ≈0 no implica que Hi-C sea irrelevante en general. En este training set ~90% de guías tienen hic_within_5kb_tss (zero-pad Sec. 6.3) y solo 2 features F4 aceptadas; la señal distal (~Gasperini) se diluye en el LOCO global. Coherente con Semana I.
