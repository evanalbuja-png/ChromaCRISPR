# scripts/position_only_baseline.py

import pandas as pd
import numpy as np

from sklearn.model_selection import GroupKFold
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import r2_score, mean_squared_error


# ============================================================
# DATOS
# ============================================================

df = pd.read_csv(
    "data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv"
)

df = df[df["efficacy_score"].notna()].copy()

X = df[["nearest_tss_distance"]]
y = df["efficacy_score"]
groups = df["nearest_tss_gene"]


# ============================================================
# SPLIT
# ============================================================

gkf = GroupKFold(n_splits=5)
folds = list(gkf.split(X, y, groups))


# ============================================================
# MODELOS
# ============================================================

models = {
    "Ridge": Ridge(alpha=1.0),
    "Random Forest": RandomForestRegressor(
        n_estimators=500,
        random_state=42,
        n_jobs=-1
    )
}


# ============================================================
# EVALUACIÓN
# ============================================================

results = []

for model_name, model in models.items():

    fold_results = []

    for fold, (train_idx, test_idx) in enumerate(folds, start=1):

        X_train = X.iloc[train_idx].copy()
        X_test = X.iloc[test_idx].copy()

        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]

        # Imputación usando únicamente training
        median = X_train["nearest_tss_distance"].median()

        X_train["nearest_tss_distance"] = (
            X_train["nearest_tss_distance"].fillna(median)
        )

        X_test["nearest_tss_distance"] = (
            X_test["nearest_tss_distance"].fillna(median)
        )

        model.fit(X_train, y_train)

        y_pred = model.predict(X_test)

        spearman = spearmanr(
            y_test,
            y_pred
        ).statistic

        pearson_r = pearsonr(
            y_test,
            y_pred
        ).statistic

        pearson_r2 = pearson_r ** 2

        r2 = r2_score(
            y_test,
            y_pred
        )

        rmse = np.sqrt(
            mean_squared_error(
                y_test,
                y_pred
            )
        )

        fold_results.append({
            "model": model_name,
            "fold": fold,
            "spearman": spearman,
            "pearson_r": pearson_r,
            "pearson_r2": pearson_r2,
            "r2": r2,
            "rmse": rmse
        })

        print(
            f"{model_name} | Fold {fold} | "
            f"Spearman: {spearman:.6f} | "
            f"Pearson R²: {pearson_r2:.6f} | "
            f"R²: {r2:.6f} | "
            f"RMSE: {rmse:.6f}"
        )

    results.extend(fold_results)


# ============================================================
# RESUMEN
# ============================================================

results_df = pd.DataFrame(results)

summary = (
    results_df
    .groupby("model")
    .agg({
        "spearman": "mean",
        "pearson_r": "mean",
        "pearson_r2": "mean",
        "r2": "mean",
        "rmse": "mean"
    })
    .reset_index()
)

print("\n=== POSITION-ONLY BASELINE ===")
print(summary.to_string(index=False))


# ============================================================
# GUARDAR RESULTADOS
# ============================================================

results_df.to_csv(
    "results/week9_position_baseline_folds.csv",
    index=False
)

summary.to_csv(
    "results/week9_position_baseline_summary.csv",
    index=False
)


# ============================================================
# LOG
# ============================================================

with open(
    "logs/week9_position_baseline.md",
    "w"
) as f:

    f.write("# Week 9 — Position-only Baseline\n\n")

    f.write("## Dataset\n\n")
    f.write(
        "- Dataset: Horlbeck2016\n"
        "- Target: `efficacy_score`\n"
        "- Predictor: `nearest_tss_distance`\n"
        "- Samples: efficacy_score non-null subset\n"
        "- Groups: `nearest_tss_gene`\n"
        "- Cross-validation: 5-fold GroupKFold\n\n"
    )

    f.write("## Models\n\n")
    f.write(
        "- Ridge Regression (`alpha=1.0`)\n"
        "- Random Forest (`n_estimators=500`, `random_state=42`, `n_jobs=-1`)\n\n"
    )

    f.write("## Results\n\n")
    f.write(
        "| Model | Spearman ρ | Pearson R | Pearson R² | R² | RMSE |\n"
        "|---|---:|---:|---:|---:|---:|\n"
    )

    for _, row in summary.iterrows():
        f.write(
            f"| {row['model']} | "
            f"{row['spearman']:.6f} | "
            f"{row['pearson_r']:.6f} | "
            f"{row['pearson_r2']:.6f} | "
            f"{row['r2']:.6f} | "
            f"{row['rmse']:.6f} |\n"
        )

    f.write("\n## Comparison with best complete model\n\n")
    f.write(
        "Best complete model: XGBoost F1+F2+F3+F5, "
        "Spearman ρ ≈ 0.322.\n"
    )