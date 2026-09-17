# Semana Q — Inventario benchmarking (Sec. 12)

## Tabla de viabilidad

| Tool | Sec. | Disponibilidad | Notas |
|------|------|----------------|-------|
| Position-baseline | 12.3 | **Sin fricción** | Feature `log_dist_target_tss` (y nearest expressed) en `K562_training_matrix_v1.csv`. Artefactos Phase1: `results/week9_position_baseline_*.csv`. |
| CRISPOR Rule Set 2 | 12.1 | **Esfuerzo moderado** | REST `https://crispor.tefor.net` → Connection reset. `http://crispor.tefor.net` → 200 (HTTP only). **Alternativa local:** PyPI `azimuth==2.0` (Doench RS2 / Microsoft) — no instalado. **Existente:** `results/ruleset3_horlbeck_scored.csv` es Rule Set **3**, no RS2; no sustituye Sec. 12.1 sin documentar desviación. |
| DeepCRISPR | 12.2 | **Alto esfuerzo / posible bloqueante** | GitHub `bm2-lab/DeepCRISPR` alcanzable. Stack típico TensorFlow 1.x + pesos pre-entrenados (nuclease, no CRISPRi). No hay paquete en el env ni artefactos locales. Requiere inventario de deps/tamaño de pesos **antes** de instalar. |

## Recomendación para esta semana
1. **Implementar ya:** position-baseline (Spearman LOCO o test vs EEP; comparar |ρ| con XGB 0.377).
2. **Implementar con instalación controlada:** Azimuth (RS2) sobre secuencias K562 — reportar tamaño de install; si falla en py3.10, documentar y evaluar scores RS3 ya calculados solo como **desviación documentada**, no como RS2.
3. **DeepCRISPR:** un solo paso de factibilidad (clone shallow + README deps + tamaño pesos). Si TF1 o pesos inaccesibles → **bloqueante documentado** (mismo estándar Replogle/Condition B), no forzar.

## No hecho aún
Instalaciones, llamadas API masivas, ni tablas finales Sec. 12.
