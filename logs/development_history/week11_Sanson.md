## Condition A — Modelo K562 + features K562 aplicados a Sanson2018 (HT29 / A375)

### Diseño
Según Sección 10.2 del proposal: entrenar en K562 (Horlbeck), extraer features de cromatina K562 en coordenadas de las guías Sanson, predecir y correlacionar contra LFC observado en HT29 y A375. El mismatch de línea celular es parte del diseño.

### Corrección de pipeline (H3K4me3)
Se detectó un bug en `extract_f3_h3k4me3_features.py`: tras `computeMatrix --skipZeros`, los valores se asignaban al BED por **índice posicional** en lugar de usar el campo `name` de la matriz. Eso desalineó H3K4me3 respecto a las coordenadas reales (Spearman stored vs re-extract ≈ 0.02–0.03).

**Fix:** re-extracción de H3K4me3 y H3K27me3 sobre `sgRNA_hg38.bed` (±500 bp, bin 10), unión por `chrom:start-end` (name de deepTools), reconstrucción de la matriz consolidada y re-entrenamiento del XGBoost. ATAC y H3K27ac ya eran reproducibles (ρ = 1.0).

### Métricas Condition A (modelo FIXED)

| Modelo | N | HT29 Spearman | A375 Spearman |
|--------|---|---------------|---------------|
| Sequence-only (F1+F5) | 100,942 | −0.033 | −0.062 |
| Full (F1+F2+F3+F5, cromatina K562) | 100,942 | −0.216 | −0.238 |

Polaridad: `efficacy_score` alto en Horlbeck correlaciona negativamente con LFC Sanson (ρ ≈ −0.35 en el overlap de 526 guías idénticas). El signo negativo de las predicciones es el esperado para CRISPRi.

### Interpretación
1. **Aporte de cromatina K562 en transferencia:** el modelo full supera de forma clara al sequence-only en el set completo de Sanson (|ρ| ≈ 0.22 vs ≈ 0.05). Las features de accesibilidad/histonas de K562 aportan señal transferible más allá de secuencia y distancia al TSS.
2. **Magnitud vs proposal (H3):** el proposal anticipaba ρ ≈ 0.35–0.45 para Condition A. El valor observado (|ρ| ≈ 0.22–0.24) está por debajo; se reporta como discrepancia real, no como fallo de ejecución.
3. **Overlap 526:** el score experimental crudo de Horlbeck (ρ ≈ −0.35) sigue siendo superior al modelo (ρ ≈ −0.13). El modelo no captura toda la información fenotípica disponible en el training set para esas guías.
4. **OOF interno K562:** sequence-only 0.282 vs full 0.297 → ganancia marginal de cromatina *dentro* de la misma línea; la ganancia relativa es mayor en el setting de transferencia (Condition A).

### Artefactos
- Matriz FIXED: `data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_FIXED.csv`
- Modelo full: `models/xgb_phase1_horlbeck_FIXED.joblib`
- Modelo sequence-only: `models/xgb_sequence_only_horlbeck.joblib`
- Predicciones Condition A: `results/conditionA_sanson_predictions_FIXED.csv`
- Resumen seq vs full: `results/conditionA_sequence_vs_full_summary.json`

### No realizado en este bloque
- Condition B (epigenoma matched HT29/A375): cobertura incompleta de las marcas de histona en ENCODE para esas líneas; no se fuerza una Condition B parcial.