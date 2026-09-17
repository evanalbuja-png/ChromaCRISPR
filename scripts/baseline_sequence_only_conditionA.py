#!/usr/bin/env python
# scripts/baseline_sequence_only_conditionA.py
"""
Baseline sequence-only (F1 + F5) entrenado en Horlbeck FIXED,
evaluado en Sanson (Condition A design).
Compara contra el modelo full (cromatina + secuencia).
"""

from pathlib import Path
import json
import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.model_selection import GroupKFold
from xgboost import XGBRegressor

HORLBECK = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_FIXED.csv")
SANSON = Path("data/interim/conditionA/sanson_feature_matrix_conditionA.csv")
OUT_DIR = Path("models")
RESULTS = Path("results")
OUT_DIR.mkdir(exist_ok=True)
RESULTS.mkdir(exist_ok=True)

SEQ_FEATURES = [
    "gc_content", "mfe_rnafold", "guide_length", "g_run_max", "poly_t_count",
    "nearest_tss_distance", "within_promoter_2kb",
]
TARGET = "efficacy_score"
GROUP = "chromosome"
RANDOM_STATE = 42

XGB_PARAMS = {
    "n_estimators": 400,
    "max_depth": 5,
    "learning_rate": 0.05,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "min_child_weight": 5,
    "reg_alpha": 0.1,
    "reg_lambda": 1.0,
    "objective": "reg:squarederror",
    "n_jobs": -1,
    "random_state": RANDOM_STATE,
    "tree_method": "hist",
}

def main():
    print("=== 1. Entrenar sequence-only en Horlbeck FIXED ===")
    h = pd.read_csv(HORLBECK)
    h = h.dropna(subset=[TARGET] + SEQ_FEATURES).copy()
    print(f"N train: {len(h):,}")

    X = h[SEQ_FEATURES]
    y = h[TARGET]
    groups = h[GROUP]

    gkf = GroupKFold(n_splits=5)
    oof = np.zeros(len(X))
    for fold, (tr, va) in enumerate(gkf.split(X, y, groups), 1):
        m = XGBRegressor(**XGB_PARAMS)
        m.fit(X.iloc[tr], y.iloc[tr], verbose=False)
        oof[va] = m.predict(X.iloc[va])
        rho = spearmanr(y.iloc[va], oof[va]).correlation
        print(f"  Fold {fold}: Spearman={rho:.4f}")

    rho_oof = spearmanr(y, oof).correlation
    print(f"OOF sequence-only (Horlbeck): Spearman={rho_oof:.4f}")

    model = XGBRegressor(**XGB_PARAMS)
    model.fit(X, y, verbose=False)
    joblib.dump(model, OUT_DIR / "xgb_sequence_only_horlbeck.joblib")
    with open(OUT_DIR / "xgb_sequence_only_features.json", "w") as f:
        json.dump(SEQ_FEATURES, f, indent=2)

    print("\n=== 2. Predecir en Sanson ===")
    s = pd.read_csv(SANSON)
    Xs = s[SEQ_FEATURES].copy().fillna(s[SEQ_FEATURES].median())
    s["pred_seq_only"] = model.predict(Xs)

    # Cargar predicciones full del modelo FIXED
    full_pred = pd.read_csv("results/conditionA_sanson_predictions_FIXED.csv")
    s = s.merge(
        full_pred[["guide_sequence", "pred_efficacy"]],
        on="guide_sequence", how="left"
    )
    s = s.rename(columns={"pred_efficacy": "pred_full"})

    print("\n=== 3. Correlaciones Condition A ===")
    print(f"{'Modelo':<20} {'Cell':<8} {'Spearman':>10} {'Pearson':>10}")
    print("-" * 52)

    summary = {"sequence_only_oof_horlbeck": float(rho_oof), "conditionA": {}}

    for label, col_pred in [("sequence-only", "pred_seq_only"), ("full (chromatin)", "pred_full")]:
        summary["conditionA"][label] = {}
        for cell, col_y in [("HT29", "lfc_HT29"), ("A375", "lfc_A375")]:
            mask = s[col_y].notna() & s[col_pred].notna()
            rho, _ = spearmanr(s.loc[mask, col_y], s.loc[mask, col_pred])
            r, _ = pearsonr(s.loc[mask, col_y], s.loc[mask, col_pred])
            print(f"{label:<20} {cell:<8} {rho:>10.4f} {r:>10.4f}")
            summary["conditionA"][label][cell] = {
                "spearman": float(rho), "pearson": float(r), "n": int(mask.sum())
            }

    # Overlap
    horl = h[["guide_sequence", "efficacy_score"]]
    ov = s.merge(horl, on="guide_sequence", how="inner")
    print(f"\n--- Overlap n={len(ov)} ---")
    for cell in ["lfc_HT29", "lfc_A375"]:
        rho_seq, _ = spearmanr(ov["pred_seq_only"], ov[cell], nan_policy="omit")
        rho_full, _ = spearmanr(ov["pred_full"], ov[cell], nan_policy="omit")
        rho_raw, _ = spearmanr(ov["efficacy_score"], ov[cell], nan_policy="omit")
        print(f"  {cell}: seq-only={rho_seq:.4f}  full={rho_full:.4f}  raw_score={rho_raw:.4f}")

    with open(RESULTS / "conditionA_sequence_vs_full_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    s[["guide_sequence", "lfc_HT29", "lfc_A375", "pred_seq_only", "pred_full"]].to_csv(
        RESULTS / "conditionA_sequence_vs_full_predictions.csv", index=False
    )
    print("\n=== BASELINE SEQUENCE-ONLY COMPLETADO ===")

if __name__ == "__main__":
    main()