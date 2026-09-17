#!/usr/bin/env python
# scripts/patch_phase2_report_paso125.py
"""Añade sección Paso 12.5 al reporte final Condition A."""

from pathlib import Path

REPORT = Path("logs/phase2_condition_a_final_report.md")
text = REPORT.read_text(encoding="utf-8")

section = """

---

## 9. Paso 12.5 — Alineación explícita con el proposal (desviaciones documentadas)

### 9.1 Líneas celulares de validación externa

El proposal (Sec. 7.1 DS2; Sec. 10.2) especifica validación externa en **Jurkat y MCF7** (Sanson et al. 2018).

**Datos efectivamente disponibles en el proyecto:** suplementario Dolcetto (`41467_2018_7901_MOESM6_ESM.xlsx`) con LFC guía-nivel en **HT29 y A375**. No se encontró archivo procesable con Jurkat/MCF7 en `data/raw/crispr_datasets/Sanson2018`.

**Decisión:** se sustituye Jurkat/MCF7 por HT29/A375 como líneas de Condition A, por disponibilidad de fenotipos. El diseño experimental (modelo K562 + features K562 → otra línea) se mantiene; solo cambia la identidad de la línea de destino. Esta sustitución es explícita y no implícita.

### 9.2 Dimensionalidad de features vs proposal

| Grupo | Proposal (Sec. 6.3) | Phase 1 FIXED usado |
|-------|---------------------|---------------------|
| F1 secuencia | 24 | 5 |
| F2 accesibilidad | 18 | 4 |
| F3 histonas | 42 | 12 |
| F4 Hi-C | 8 | 0 (no incluido) |
| F5 contexto genómico | 12 | 2 |
| **Total** | **104** | **23** |

Las expectativas de ρ para Condition A en Sec. 11.3 (≈0.35–0.45) se formularon asumiendo el modelo completo de 104 features. El ρ observado (|ρ|≈0.22–0.24) debe interpretarse en ese contexto: parte de la brecha respecto al proposal es atribuible a la reducción de features de Phase 1, no solo a la transferencia cross-cell-type.

### 9.3 Test de permutación (Sec. 10.2)

Comparación apareada full vs sequence-only, set completo Sanson (N=100,942), 10,000 permutaciones, contraste one-sided sobre Δ|ρ|:

| Cell | ρ seq-only | ρ full | Δ\|ρ\| | p |
|------|------------|--------|-------|---|
| HT29 | −0.0333 | −0.2163 | 0.1829 | 0.0001 |
| A375 | −0.0622 | −0.2377 | 0.1755 | 0.0001 |

La superioridad del modelo con cromatina K562 sobre el baseline de secuencia es estadísticamente significativa en ambas líneas.

Artefacto: `results/conditionA_permutation_test.json`.

### 9.4 Split del OOF interno K562

Los OOF reportados (sequence-only ρ=0.282; full ρ=0.297) provienen de **GroupKFold (5 folds) con `groups=chromosome`**, no de split aleatorio por guía.

No se ejecutó el LOCO exhaustivo de 23 folds (chr1–22+X) descrito en Sec. 7.3 / 10.1. GroupKFold por cromosoma es la aproximación usada en Phase 1 para Condition A; se declara como limitación heredada respecto al protocolo LOCO completo del proposal, no como validación inflada por autocorrelación espacial de un split aleatorio.

*(Los baselines de week 8 que usaban `train_test_split` aleatorio no son la fuente de estos OOF.)*

### 9.5 Bootstrap CIs (Sec. 10.5)

No calculados en este cierre. Quedan en backlog para el manuscrito si se requieren IC de 1,000 resamples sobre ρ.

---
"""

marker = "*Fin del reporte."
if section.strip() in text:
    print("Sección 9 ya presente; no se duplica.")
elif marker in text:
    text = text.replace(marker, section + "\n" + marker)
    REPORT.write_text(text, encoding="utf-8")
    print(f"Parche aplicado: {REPORT}")
else:
    text = text.rstrip() + "\n" + section
    REPORT.write_text(text, encoding="utf-8")
    print(f"Sección 9 añadida al final: {REPORT}")

print("=== PASO 12.5 DOCUMENTADO EN REPORTE ===")