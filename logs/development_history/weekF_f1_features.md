# Semana F — F1 Sequence Features (cerrada)

**Fecha:** 2026-09-10  
**Base:** 19.446 guías de K562_training_labels_harmonized_v1.csv

## Contabilidad vs Sec. 6.3

| Feature Sec. 6.3 | Estado | Nombre de columna | Notas |
|---|---|---|---|
| GC content (fraction) | Implementada | gc_content | Verificada |
| G runs (max homopolymer) | Implementada | g_run_max | Verificada |
| poly-T runs (TTTT count) | Implementada | poly_t_count | Verificada |
| MFE (RNAfold ViennaRNA) | Implementada | mfe_rnafold | Re-calculada con binario RNAfold; coincide con FIXED |
| Melting temperature (SantaLucia 1998) | Implementada | tm_santalucia | Biopython Tm_NN (DNA_NN4) |
| Position-weight matrix seed region | Implementada | seed_pwm_score | PWM no supervisado (frecuencias posicionales sobre las 19.446 guías, sin condicionar en EEP). Seed = pos 1–8. |
| PAM context (4 nt surrounding NGG, one-hot) | **Pendiente** | — | Requiere secuencia genómica flanqueante. No disponible en esta etapa. |
| guide_length | Extra (no en Sec. 6.3) | guide_length | Conservada como extra no oficial |
| Proxy extremo 3′ (16 cols one-hot) | Proxy documentado | guide_3prime_pos{1-4}_{A,C,G,T} | **No es PAM context**. Nombre honesto. No cuenta para las 24. |

**Oficiales completas: 6 / 24.**  
PAM context real queda pendiente (se completará cuando existan flancos genómicos).

## Verificación de control (3 guías)
Todas las features existentes coincidieron exactamente entre la matriz FIXED y el recálculo independiente (gc, g_run, poly_t, length, mfe).

## Sanity check PWM (anti-leakage)
Spearman(seed_pwm_score, EEP_percentile) = **0.0746** (p = 1.25e-25).  
Rho bajo → sin evidencia de filtración de la label. Aceptado.

## Archivo generado
`data/interim/F1_sequence_features_v1.csv` — 19.446 filas × 24 columnas (1 id + 6 oficiales + 1 extra + 16 proxy).

## Criterio de éxito
CUMPLIDO — features honestamente contabilizadas, sin leakage, auditoría transparente.
