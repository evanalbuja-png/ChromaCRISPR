#!/usr/bin/env python
# scripts/predict_conditionA_sanson.py
"""
Condition A: modelo K562 (Horlbeck) + features K562 aplicados a guías Sanson
→ correlacionar predicciones vs LFC real en HT29 y A375.
"""

from pathlib import Path
import json
import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import r2_score, mean_squared_error

MODEL_PATH = Path("models/xgb_phase1_horlbeck_conditionA.joblib")
FEATURES_PATH = Path("models/xgb_phase1_feature_list.json")
MATRIX_PATH = Path("data/interim/conditionA/sanson_feature_matrix_conditionA.csv")
OUT_PRED = Path("results/conditionA_sanson_predictions.csv")
OUT_SUMMARY = Path("results/conditionA_sanson_summary.json")

def main():
    print("=== 1. Cargando modelo y features ===")
    model = joblib.load(MODEL_PATH)
    with open(FEATURES_PATH) as f:
        features = json.load(f)
    print(f"Modelo: {MODEL_PATH}")
    print(f"Features ({len(features)}): {features}")

    print("\n=== 2. Cargando matriz Sanson ===")
    df = pd.read_csv(MATRIX_PATH)
    print(f"Shape: {df.shape}")

    # Imputar missing con la mediana del propio set (robusto y simple)
    X = df[features].copy()
    medians = X.median()
    X = X.fillna(medians)
    print(f"Missing residuales tras imputación: {X.isna().sum().sum()}")

    print("\n=== 3. Predicción ===")
    preds = model.predict(X)
    df["pred_efficacy"] = preds

    # Guardar predicciones
    cols_out = ["region_id", "guide_sequence", "gene_symbol", "chromosome",
                "lfc_HT29", "lfc_A375", "pred_efficacy"] + features
    df[cols_out].to_csv(OUT_PRED, index=False)
    print(f"Predicciones guardadas: {OUT_PRED}")

    # --------------------------------------------------
    print("\n=== 4. Correlaciones Condition A ===")
    results = {}

    for cell, col in [("HT29", "lfc_HT29"), ("A375", "lfc_A375")]:
        mask = df[col].notna() & np.isfinite(df[col])
        y = df.loc[mask, col].values
        p = df.loc[mask, "pred_efficacy"].values

        rho, pval_rho = spearmanr(y, p)
        r, pval_r = pearsonr(y, p)
        r2 = r2_score(y, p)
        rmse = np.sqrt(mean_squared_error(y, p))

        results[cell] = {
            "n": int(mask.sum()),
            "spearman_rho": float(rho),
            "spearman_pval": float(pval_rho),
            "pearson_r": float(r),
            "pearson_pval": float(pval_r),
            "r2": float(r2),
            "rmse": float(rmse),
        }
        print(f"\n--- {cell} ---")
        print(f"  N          : {mask.sum():,}")
        print(f"  Spearman ρ : {rho:.4f}  (p = {pval_rho:.2e})")
        print(f"  Pearson r  : {r:.4f}  (p = {pval_r:.2e})")
        print(f"  R²         : {r2:.4f}")
        print(f"  RMSE       : {rmse:.4f}")

    # También correlación entre los dos LFC (control de calidad)
    mask2 = df["lfc_HT29"].notna() & df["lfc_A375"].notna()
    rho_cells, _ = spearmanr(df.loc[mask2, "lfc_HT29"], df.loc[mask2, "lfc_A375"])
    print(f"\nSpearman LFC_HT29 vs LFC_A375 (mismo set): {rho_cells:.4f}")

    with open(OUT_SUMMARY, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResumen guardado: {OUT_SUMMARY}")
    print("\n=== CONDITION A COMPLETADA ===")

if __name__ == "__main__":
    main()