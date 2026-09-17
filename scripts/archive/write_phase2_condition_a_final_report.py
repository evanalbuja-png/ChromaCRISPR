#!/usr/bin/env python
# scripts/write_phase2_condition_a_final_report.py
"""
Consolida Condition A + diagnóstico + contexto week10/11 en
logs/phase2_condition_a_final_report.md
"""

from pathlib import Path
from datetime import date

OUT = Path("logs/phase2_condition_a_final_report.md")
WEEK10 = Path("logs/week10_sanson_efficacy.md")
WEEK11 = Path("logs/week11_encode_coverage_ht29_a375.md")
DIAG = Path("logs/week12_overlap_diagnostic.md")

def read_if_exists(path: Path) -> str:
    if path.exists():
        return path.read_text(encoding="utf-8").strip()
    return f"*(Archivo {path} no encontrado en el momento de consolidar.)*"

def main():
    week10 = read_if_exists(WEEK10)
    week11 = read_if_exists(WEEK11)
    diag = read_if_exists(DIAG)

    report = f"""# Phase 2 — Condition A: Final Report

**Proyecto:** ChromaCRISPR Phase 1  
**Documento:** Fuente única de verdad para Condition A (Sección 10.2 del proposal)  
**Fecha de consolidación:** {date.today().isoformat()}  
**Estado:** Cerrado — no reabrir decisiones ya tomadas

---

## 1. Diseño experimental (Condition A)

Según la Sección 10.2 del research proposal:

> Aplicar el modelo entrenado en K562 (Horlbeck) **con features de cromatina K562** a guías evaluadas en otra línea celular, y correlacionar predicciones contra el LFC experimental observado.

**Implementación ejecutada:**

| Elemento | Detalle |
|----------|---------|
| Training | Horlbeck 2016 CRISPRi K562 (`efficacy_score`), matriz FIXED |
| Modelo | XGBoost (23 features: F1 + F2 + F3 + F5) |
| Guías de test | Sanson 2018, coordenadas hg38 válidas |
| N test | 100,942 guías |
| Features de cromatina | ATAC, H3K27ac, H3K4me3, H3K27me3 de **K562** (ENCODE), ventana ±500 bp, bin 10 bp |
| Outcomes | `lfc_HT29`, `lfc_A375` |
| Control | Modelo sequence-only (F1 + F5, sin cromatina), mismo diseño |
| Overlap limpio | 526 guías con secuencia idéntica Horlbeck ∩ Sanson |

El mismatch de línea celular (features K562 → fenotipo HT29/A375) es **parte del diseño**, no un defecto de ejecución.

---

## 2. Bug de pipeline encontrado y corregido (H3K4me3)

### Hallazgo

En `scripts/extract_f3_h3k4me3_features.py`, tras `computeMatrix --skipZeros`, los valores de señal se asignaban al BED por **índice posicional**:

```python
atac["region_id"] = bed["region_id"].iloc[:len(atac)].values  # incorrecto
```

en lugar de usar el campo `name` de la matriz deepTools. Con filas omitidas por `skipZeros`, H3K4me3 quedó desalineado respecto a las coordenadas reales.

**Evidencia:** Spearman entre valores almacenados y re-extracción en las mismas coordenadas ≈ **0.02–0.03** (prácticamente ruido). ATAC y H3K27ac eran reproducibles (ρ = 1.0).

### Fix

1. Re-extracción de H3K4me3 y H3K27me3 sobre `sgRNA_hg38.bed` (±500 bp, bin 10).
2. Unión por `chrom:start-end` (name de deepTools).
3. Reconstrucción de `sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_FIXED.csv`.
4. Re-entrenamiento del XGBoost → `models/xgb_phase1_horlbeck_FIXED.joblib`.

OOF GroupKFold post-fix: Spearman ≈ **0.297** (casi idéntico al pre-fix ≈ 0.30). H3K4me3 aportaba poco a la predicción interna en K562, pero la feature debía ser correcta antes de interpretar biología o transferencia.

---

## 3. Métricas finales

### 3.1 Dentro de K562 (validación interna, GroupKFold)

| Modelo | OOF Spearman ρ |
|--------|----------------|
| Sequence-only (F1 + F5) | 0.282 |
| Full (F1 + F2 + F3 + F5) | 0.297 |

Ganancia marginal de cromatina *dentro* de la misma línea celular.

### 3.2 Condition A — set completo Sanson (N = 100,942)

| Modelo | HT29 Spearman | A375 Spearman | HT29 Pearson | A375 Pearson |
|--------|---------------|---------------|--------------|--------------|
| Sequence-only | −0.033 | −0.062 | −0.082 | −0.096 |
| **Full (cromatina K562)** | **−0.216** | **−0.238** | **−0.256** | **−0.269** |

**Polaridad:** `efficacy_score` alto (Horlbeck) correlaciona negativamente con LFC Sanson en el overlap (ρ ≈ −0.35). El signo negativo de las predicciones es el esperado para CRISPRi (LFC más negativo = mayor eficacia).

### 3.3 Overlap de secuencia idéntica (N = 526)

| Comparación | HT29 ρ | A375 ρ |
|-------------|--------|--------|
| efficacy_score crudo vs LFC | −0.363 | −0.352 |
| Modelo full vs LFC | −0.130 | −0.157 |
| Sequence-only vs LFC | −0.149 | −0.159 |
| Modelo full vs efficacy_score | 0.298 | — |

---

## 4. Interpretación principal

1. **La cromatina K562 aporta señal en transferencia.**  
   Full \|ρ\| ≈ 0.22–0.24 vs sequence-only \|ρ\| ≈ 0.03–0.06 en el set completo. Esa diferencia es la evidencia central de Condition A a favor del argumento de ChromaCRISPR: el contexto epigenómico mejora la predicción de eficacia CRISPRi más allá de secuencia y distancia al TSS, incluso cuando el epigenoma no es de la línea de destino.

2. **Discrepancia vs proposal (H3 / Sección 11.3).**  
   El proposal anticipaba ρ ≈ 0.35–0.45 para Condition A. El valor observado (\|ρ\| ≈ 0.22–0.24) está por debajo. Se reporta como resultado empírico real, no como fallo de ejecución ni como artefacto del bug (el bug fue corregido antes de estas métricas finales).

3. **Anomalía del overlap (modelo ≪ score crudo).**  
   Diagnosticada en `logs/week12_overlap_diagnostic.md`. Resumen:
   - Missingness en overlap ≈ full (2.85% vs 3.22%) → **no** es imputación diferencial.
   - Spearman(pred, efficacy_score) en overlap = 0.298 ≈ OOF global → el modelo **sí** reproduce el target Horlbeck en el subset.
   - El modelo comprime extremos (regresión a la media; MAE residual ≈ 0.41 en overlap).
   - **Causa clasificada:** optimización del ranking global / regresión a la media. El score experimental crudo conserva correlación punto-a-punto con LFC Sanson que el regresor, entrenado sobre ~18k guías, no preserva igual de bien en este subset de transferencia.

---

## 5. Decisión sobre Condition B

**Condition B no se ejecuta en esta fase.**

Motivo (week 11): cobertura ENCODE insuficiente para HT29 y A375 respecto al set completo de marcas usado en el modelo K562 (ATAC + H3K27ac + H3K4me3 + H3K27me3). Forzar una Condition B parcial (p. ej. solo DNase-seq de HT29) compararía un epigenoma incompleto contra el modelo full-feature de K562 y **no** sería Condition B según el diseño del proposal (Sección 10.2 exige epigenoma matched completo).

Cualquier exploración futura con cobertura parcial debe etiquetarse como análisis exploratorio distinto, no como Condition B.

---

## 6. Artefactos de referencia

| Artefacto | Ruta |
|-----------|------|
| Matriz training FIXED | `data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_FIXED.csv` |
| Modelo full FIXED | `models/xgb_phase1_horlbeck_FIXED.joblib` |
| Modelo sequence-only | `models/xgb_sequence_only_horlbeck.joblib` |
| Matriz features Sanson (Condition A) | `data/interim/conditionA/sanson_feature_matrix_conditionA.csv` |
| Predicciones Condition A | `results/conditionA_sanson_predictions_FIXED.csv` |
| Resumen seq vs full | `results/conditionA_sequence_vs_full_summary.json` |
| Diagnóstico overlap | `logs/week12_overlap_diagnostic.md` |
| Tabla overlap | `results/week12_overlap_diagnostic_table.csv` |

---

## 7. Contexto previo (week 10 y week 11)

### 7.1 Week 10 — Sanson efficacy / overlap

{week10}

### 7.2 Week 11 — Cobertura ENCODE HT29 / A375

{week11}

### 7.3 Week 12 — Diagnóstico overlap (extracto)

{diag}

---

## 8. Mensaje para manuscrito / research statement (1 párrafo)

Condition A del diseño Phase 1 evaluó un modelo XGBoost entrenado en eficacia CRISPRi K562 (Horlbeck) usando features de secuencia, accesibilidad y marcas de histona de K562, aplicado a 100,942 guías Sanson 2018 con fenotipos en HT29 y A375. Tras corregir un error de alineación en la feature H3K4me3 del corpus de entrenamiento, el modelo multimodal alcanzó Spearman ρ = −0.22 (HT29) y −0.24 (A375) frente a LFC experimental, superando de forma clara un baseline solo de secuencia (ρ ≈ −0.03 a −0.06). El signo negativo es consistente con la polaridad del outcome CRISPRi. La magnitud está por debajo del rango anticipado en el proposal (0.35–0.45); se interpreta como límite real de transferencia con epigenoma no matched. Condition B (epigenoma matched a la línea de destino) no se ejecutó por cobertura ENCODE incompleta de las marcas requeridas en HT29/A375.

---

*Fin del reporte. Decisiones de diseño, fix de H3K4me3, métricas Condition A, clasificación de la anomalía del overlap y exclusión de Condition B quedan fijadas en este documento.*
"""

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(report, encoding="utf-8")
    print(f"Escrito: {OUT}")
    print(f"Tamaño: {OUT.stat().st_size:,} bytes")
    print("=== REPORTE FINAL CONDITION A CONSOLIDADO ===")

if __name__ == "__main__":
    main()
