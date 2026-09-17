#!/usr/bin/env python
# scripts/fix5_conditionA_with_fixed_model.py
"""
Condition A: modelo FIXED (K562) + features K562 sobre guías Sanson
→ correlacionar vs LFC HT29 / A375.
"""

from pathlib import Path
import json
import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import r2_score, mean_squared_error

MODEL = Path("models/xgb_phase1_horlbeck_FIXED.joblib")
FEATURES = Path("models/xgb_phase1_feature_list_FIXED.json")
MATRIX = Path("data/interim/conditionA/sanson_feature_matrix_conditionA.csv")
OUT_PRED = Path("results/conditionA_sanson_predictions_FIXED.csv")
OUT_SUMMARY = Path("results/conditionA_sanson_summary_FIXED.json")

def main():
    model = joblib.load(MODEL)
    with open(FEATURES) as f:
        feats = json.load(f)

    df = pd.read_csv(MATRIX)
    print(f"Sanson matrix: {df.shape}")

    X = df[feats].copy()
    X = X.fillna(X.median())
    df["pred_efficacy"] = model.predict(X)

    cols = ["region_id", "guide_sequence", "gene_symbol", "chromosome",
            "lfc_HT29", "lfc_A375", "pred_efficacy"] + feats
    df[[c for c in cols if c in df.columns]].to_csv(OUT_PRED, index=False)

    results = {}
    print(f"\n{'Cell':<8} {'N':>8} {'Spearman':>10} {'Pearson':>10} {'R²':>10}")
    print("-" * 50)
    for cell, col in [("HT29", "lfc_HT29"), ("A375", "lfc_A375")]:
        mask = df[col].notna() & np.isfinite(df[col])
        y = df.loc[mask, col].values
        p = df.loc[mask, "pred_efficacy"].values
        rho, p_rho = spearmanr(y, p)
        r, p_r = pearsonr(y, p)
        r2 = r2_score(y, p)
        rmse = np.sqrt(mean_squared_error(y, p))
        results[cell] = {
            "n": int(mask.sum()),
            "spearman_rho": float(rho),
            "spearman_pval": float(p_rho),
            "pearson_r": float(r),
            "r2": float(r2),
            "rmse": float(rmse),
        }
        print(f"{cell:<8} {mask.sum():>8,} {rho:>10.4f} {r:>10.4f} {r2:>10.4f}")

    # Overlap 526: modelo vs score crudo
    horlbeck = pd.read_csv(
        "data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_FIXED.csv",
        usecols=["guide_sequence", "efficacy_score"]
    )
    ov = df.merge(horlbeck, on="guide_sequence", how="inner")
    print(f"\n--- Overlap n={len(ov)} ---")
    for cell in ["lfc_HT29", "lfc_A375"]:
        rho_pred, _ = spearmanr(ov["pred_efficacy"], ov[cell], nan_policy="omit")
        rho_raw, _ = spearmanr(ov["efficacy_score"], ov[cell], nan_policy="omit")
        print(f"  {cell}: modelo ρ={rho_pred:.4f}  |  efficacy_score crudo ρ={rho_raw:.4f}")

    with open(OUT_SUMMARY, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nGuardado: {OUT_PRED}")
    print(f"Resumen:  {OUT_SUMMARY}")
    print("=== CONDITION A (FIXED) COMPLETADA ===")

if __name__ == "__main__":
    main()