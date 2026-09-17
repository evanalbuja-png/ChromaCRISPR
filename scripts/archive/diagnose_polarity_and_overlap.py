#!/usr/bin/env python
# scripts/diagnose_polarity_and_overlap.py
"""
1. Polaridad del efficacy_score vs LFC
2. Correlación directa en el overlap de 551 guías (reproducir Semana 10)
3. Comparar rendimiento del modelo SOLO en esas 551 vs en las 100k
"""

from pathlib import Path
import json
import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

# Paths
HORLBECK = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv")
SANSON_MAT = Path("data/interim/conditionA/sanson_feature_matrix_conditionA.csv")
OVERLAP = Path("data/processed/sanson2018_horlbeck_overlap.csv")  # si existe
MODEL = Path("models/xgb_phase1_horlbeck_conditionA.joblib")
FEATURES = Path("models/xgb_phase1_feature_list.json")

def main():
    print("=" * 60)
    print("1. POLARIDAD DEL TARGET (Horlbeck efficacy_score)")
    print("=" * 60)
    h = pd.read_csv(HORLBECK)
    print(f"Horlbeck shape: {h.shape}")
    print(f"efficacy_score: mean={h['efficacy_score'].mean():.3f}, "
          f"std={h['efficacy_score'].std():.3f}, "
          f"min={h['efficacy_score'].min():.3f}, max={h['efficacy_score'].max():.3f}")
    print(f"Cuartiles: {h['efficacy_score'].quantile([0.25, 0.5, 0.75]).to_dict()}")

    # Si hay LFC crudo en alguna columna de Horlbeck, correlacionar
    for col in h.columns:
        if "lfc" in col.lower() or "gamma" in col.lower() or "log2" in col.lower():
            rho, _ = spearmanr(h["efficacy_score"], h[col], nan_policy="omit")
            print(f"  Spearman efficacy_score vs {col}: {rho:.4f}")

    print("\n" + "=" * 60)
    print("2. OVERLAP DIRECTO (reproducir Semana 10)")
    print("=" * 60)

    s = pd.read_csv(SANSON_MAT)
    # Unir por secuencia exacta
    merged = h.merge(s, left_on="guide_sequence", right_on="guide_sequence",
                     how="inner", suffixes=("_h", "_s"))
    print(f"Guías con secuencia idéntica Horlbeck ∩ Sanson: {len(merged):,}")

    if len(merged) > 0:
        for cell in ["lfc_HT29", "lfc_A375"]:
            rho, p = spearmanr(merged["efficacy_score"], merged[cell], nan_policy="omit")
            print(f"  Spearman efficacy_score (Horlbeck) vs {cell}: ρ = {rho:.4f}  (n={merged[cell].notna().sum()})")

        # También correlación entre features de cromatina del mismo guide
        print("\n  Correlación de features de cromatina (mismo guide, K562 vs K562):")
        for f in ["ATAC_mean", "H3K27ac_mean", "H3K4me3_mean", "nearest_tss_distance"]:
            if f in merged.columns or f + "_h" in merged.columns:
                # después del merge con suffixes
                col_h = f if f in merged.columns else f + "_h"
                col_s = f + "_s" if f + "_s" in merged.columns else f
                if col_h in merged.columns and col_s in merged.columns:
                    rho, _ = spearmanr(merged[col_h], merged[col_s], nan_policy="omit")
                    print(f"    {f}: ρ = {rho:.4f}")

    print("\n" + "=" * 60)
    print("3. RENDIMIENTO DEL MODELO SOLO EN EL OVERLAP vs FULL SET")
    print("=" * 60)

    model = joblib.load(MODEL)
    with open(FEATURES) as f:
        feats = json.load(f)

    # Predicciones ya las tenemos en results/conditionA_sanson_predictions.csv
    pred = pd.read_csv("results/conditionA_sanson_predictions.csv")

    # Subset overlap
    if len(merged) > 0:
        overlap_seqs = set(merged["guide_sequence"])
        pred_overlap = pred[pred["guide_sequence"].isin(overlap_seqs)].copy()
        print(f"Predicciones en overlap: {len(pred_overlap):,}")

        for cell in ["lfc_HT29", "lfc_A375"]:
            rho_full, _ = spearmanr(pred[cell], pred["pred_efficacy"], nan_policy="omit")
            rho_ov, _ = spearmanr(pred_overlap[cell], pred_overlap["pred_efficacy"], nan_policy="omit")
            print(f"  {cell}:")
            print(f"    Full set (100k)     ρ = {rho_full:.4f}")
            print(f"    Solo overlap (~551) ρ = {rho_ov:.4f}")

    print("\n=== FIN DIAGNÓSTICO POLARIDAD + OVERLAP ===")

if __name__ == "__main__":
    main()