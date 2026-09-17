#!/usr/bin/env python
# scripts/fix4_retrain_xgb_on_fixed.py
"""
Re-entrena XGBoost Phase 1 sobre la matriz Horlbeck FIXED
(H3K4me3 y H3K27me3 corregidos).
"""

from pathlib import Path
import json
import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.model_selection import GroupKFold
from xgboost import XGBRegressor

INPUT = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_FIXED.csv")
OUT_DIR = Path("models")
OUT_DIR.mkdir(exist_ok=True)

MODEL_PATH = OUT_DIR / "xgb_phase1_horlbeck_FIXED.joblib"
FEATURES_PATH = OUT_DIR / "xgb_phase1_feature_list_FIXED.json"
METRICS_PATH = OUT_DIR / "xgb_phase1_train_metrics_FIXED.json"

FEATURES = [
    "gc_content", "mfe_rnafold", "guide_length", "g_run_max", "poly_t_count",
    "ATAC_mean", "ATAC_max", "ATAC_p90", "ATAC_sum",
    "H3K27ac_mean", "H3K27ac_max", "H3K27ac_p90", "H3K27ac_sum",
    "H3K4me3_mean", "H3K4me3_max", "H3K4me3_p90", "H3K4me3_sum",
    "H3K27me3_mean", "H3K27me3_max", "H3K27me3_p90", "H3K27me3_sum",
    "nearest_tss_distance", "within_promoter_2kb",
]
TARGET = "efficacy_score"
GROUP_COL = "chromosome"
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
    df = pd.read_csv(INPUT)
    df = df.dropna(subset=[TARGET] + FEATURES).copy()
    print(f"N tras dropna: {len(df):,}")

    X = df[FEATURES]
    y = df[TARGET]
    groups = df[GROUP_COL]

    print("\n=== GroupKFold 5-fold ===")
    gkf = GroupKFold(n_splits=5)
    oof = np.zeros(len(X))
    fold_metrics = []

    for fold, (tr, va) in enumerate(gkf.split(X, y, groups), 1):
        model = XGBRegressor(**XGB_PARAMS)
        model.fit(X.iloc[tr], y.iloc[tr], eval_set=[(X.iloc[va], y.iloc[va])], verbose=False)
        pred = model.predict(X.iloc[va])
        oof[va] = pred
        rho = spearmanr(y.iloc[va], pred).correlation
        r2 = r2_score(y.iloc[va], pred)
        rmse = np.sqrt(mean_squared_error(y.iloc[va], pred))
        fold_metrics.append({"fold": fold, "spearman": float(rho), "r2": float(r2), "rmse": float(rmse)})
        print(f"  Fold {fold}: Spearman={rho:.4f}  R²={r2:.4f}  RMSE={rmse:.4f}")

    rho_oof = spearmanr(y, oof).correlation
    r2_oof = r2_score(y, oof)
    rmse_oof = np.sqrt(mean_squared_error(y, oof))
    print(f"\nOOF global → Spearman={rho_oof:.4f}  R²={r2_oof:.4f}  RMSE={rmse_oof:.4f}")
    print("(Antes del fix OOF era ~0.30)")

    print("\n=== Modelo final (todo el set) ===")
    final = XGBRegressor(**XGB_PARAMS)
    final.fit(X, y, verbose=False)

    joblib.dump(final, MODEL_PATH)
    with open(FEATURES_PATH, "w") as f:
        json.dump(FEATURES, f, indent=2)
    with open(METRICS_PATH, "w") as f:
        json.dump({
            "n_samples": int(len(X)),
            "oof_spearman": float(rho_oof),
            "oof_r2": float(r2_oof),
            "oof_rmse": float(rmse_oof),
            "fold_metrics": fold_metrics,
            "xgb_params": XGB_PARAMS,
        }, f, indent=2)

    print(f"Modelo: {MODEL_PATH}")
    print("=== LISTO Fix-4 ===")

if __name__ == "__main__":
    main()