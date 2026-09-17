# Semana I — F4 Hi-C (cerrada, parcial)

**Fecha:** 2026-09-11  
**Fuente Hi-C:** 4DNFITUOMFUQ (Rao et al.), cooler 10 kb + mcool 100 kb, KR→weight multiplicativo (1/KR).

## Conteo vs n=8
Texto Sec. 6.3 ≈ 5–6 features; proposal declara n=8. Discrepancia documentada (como F3).

## Implementadas
| Feature | Estado | Evidencia |
|---|---|---|
| hic_contact_log (KR 10kb, log1p, guía–TSS) | **OK** | Spearman(dist, contact)=**−0.90** |
| hic_within_5kb_tss (zero-pad rule) | **OK** | 17.496 guías (<5kb); informativo sobre todo en Gasperini |
| E1 (eigenvector 100kb) | Calculado | Spearman(E1, GC)=**+0.59** (orientación OK) |
| compartment_A / compartment_AB | Calculado pero **NO validado biológicamente** | ATAC A (4.73) < B (5.14); MW A>B p=1.0; Spearman(E1,ATAC)=−0.10 |

## No implementadas (pendiente real)
- same_TAD (requiere Arrowhead / DI)
- loop_anchor overlap Peakachu
- distance_to_nearest_loop_anchor

## Decisión sobre compartment
Se **conservan** E1 y compartment_A en el CSV para análisis exploratorio posterior, pero se marcan como **no confiables para modelado** hasta re-validar (posible re-orientación, otra resolución, o track publicado 4DN). No se cuenta compartment dentro de las features F4 “aceptadas” para Bloque 3.

## Features F4 aceptadas para entrenamiento
1. `hic_contact_log`
2. `hic_within_5kb_tss`  
(+ E1/compartment opcionales, flagged)

## Archivo
`data/interim/F4_hic_features_v1.csv`

## Criterio de éxito
Parcialmente cumplido: contacto Hi-C con sanity de distancia fuerte; compartment intentado y **refutado** por cruzado con ATAC — documentado sin forzar.

## Corrección post-dedup
`compartment_A` mostraba 0 NaN pese a que `E1` (su fuente) tenía 786 — causado por 
comparación booleana silenciosa sobre NaN (NaN > 0 → False). Corregido propagando 
NaN de E1 a compartment_A explícitamente. Evita interpretar falsamente esas 786 
filas como "compartimento B confirmado" cuando en realidad falta el dato. 
Shape final: 19.446 filas, NaN consistentes entre E1 y compartment_A (786 cada uno).
