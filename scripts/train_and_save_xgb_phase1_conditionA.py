#!/usr/bin/env python
# scripts/train_and_save_xgb_phase1_conditionA.py
"""
Re-entrena XGBoost Phase 1 sobre la matriz Horlbeck consolidada
y lo guarda para usarlo en Condition A (Sanson HT29/A375 + features K562).
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

# --------------------------------------------------
# Configuración
# --------------------------------------------------
INPUT = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv")
OUT_DIR = Path("models")
OUT_DIR.mkdir(exist_ok=True)

MODEL_PATH = OUT_DIR / "xgb_phase1_horlbeck_conditionA.joblib"
FEATURES_PATH = OUT_DIR / "xgb_phase1_feature_list.json"
METRICS_PATH = OUT_DIR / "xgb_phase1_train_metrics.json"

FEATURES = [
    "gc_content",
    "mfe_rnafold",
    "guide_length",
    "g_run_max",
    "poly_t_count",
    "ATAC_mean",
    "ATAC_max",
    "ATAC_p90",
    "ATAC_sum",
    "H3K27ac_mean",
    "H3K27ac_max",
    "H3K27ac_p90",
    "H3K27ac_sum",
    "H3K4me3_mean",
    "H3K4me3_max",
    "H3K4me3_p90",
    "H3K4me3_sum",
    "H3K27me3_mean",
    "H3K27me3_max",
    "H3K27me3_p90",
    "H3K27me3_sum",
    "nearest_tss_distance",
    "within_promoter_2kb",
]

TARGET = "efficacy_score"
GROUP_COL = "chromosome"          # para GroupKFold / LOCO-style
N_SPLITS = 5
RANDOM_STATE = 42

# Hiperparámetros (ajustados a lo que se usó en week8/9; si tienes los de Optuna, los cambiamos después)
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
    print("=== Cargando matriz Horlbeck consolidada ===")
    df = pd.read_csv(INPUT, low_memory=False)
    print(f"Shape original: {df.shape}")

    # Filtrar filas con target y features válidos
    df = df.dropna(subset=[TARGET] + FEATURES).copy()
    print(f"Shape tras dropna (target + features): {df.shape}")

    X = df[FEATURES]
    y = df[TARGET]
    groups = df[GROUP_COL] if GROUP_COL in df.columns else None

    print(f"\nFeatures ({len(FEATURES)}): {FEATURES}")
    print(f"Target: {TARGET}")
    print(f"N muestras: {len(X):,}")

    # --------------------------------------------------
    # Validación rápida con GroupKFold (opcional pero recomendable)
    # --------------------------------------------------
    print("\n=== Validación GroupKFold (5 folds) ===")
    gkf = GroupKFold(n_splits=N_SPLITS)
    oof_preds = np.zeros(len(X))
    fold_metrics = []

    for fold, (train_idx, val_idx) in enumerate(gkf.split(X, y, groups), 1):
        X_tr, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_tr, y_val = y.iloc[train_idx], y.iloc[val_idx]

        model = XGBRegressor(**XGB_PARAMS)
        model.fit(X_tr, y_tr, eval_set=[(X_val, y_val)], verbose=False)

        pred = model.predict(X_val)
        oof_preds[val_idx] = pred

        rho = spearmanr(y_val, pred).correlation
        r2 = r2_score(y_val, pred)
        rmse = np.sqrt(mean_squared_error(y_val, pred))
        fold_metrics.append({"fold": fold, "spearman": rho, "r2": r2, "rmse": rmse})
        print(f"  Fold {fold}: Spearman={rho:.4f}  R²={r2:.4f}  RMSE={rmse:.4f}")

    # Métricas OOF globales
    rho_oof = spearmanr(y, oof_preds).correlation
    r2_oof = r2_score(y, oof_preds)
    rmse_oof = np.sqrt(mean_squared_error(y, oof_preds))
    print(f"\nOOF global → Spearman={rho_oof:.4f}  R²={r2_oof:.4f}  RMSE={rmse_oof:.4f}")

    # --------------------------------------------------
    # Entrenar modelo final sobre TODO el dataset
    # --------------------------------------------------
    print("\n=== Entrenando modelo final sobre todo el dataset ===")
    final_model = XGBRegressor(**XGB_PARAMS)
    final_model.fit(X, y, verbose=False)

    # Guardar
    joblib.dump(final_model, MODEL_PATH)
    with open(FEATURES_PATH, "w") as f:
        json.dump(FEATURES, f, indent=2)

    metrics = {
        "n_samples": int(len(X)),
        "n_features": len(FEATURES),
        "features": FEATURES,
        "oof_spearman": float(rho_oof),
        "oof_r2": float(r2_oof),
        "oof_rmse": float(rmse_oof),
        "fold_metrics": fold_metrics,
        "xgb_params": XGB_PARAMS,
    }
    with open(METRICS_PATH, "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"\nModelo guardado en: {MODEL_PATH}")
    print(f"Lista de features:   {FEATURES_PATH}")
    print(f"Métricas:            {METRICS_PATH}")
    print("=== LISTO ===")

if __name__ == "__main__":
    main()