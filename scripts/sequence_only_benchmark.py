import pandas as pd
import numpy as np

from sklearn.model_selection import GroupKFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from scipy.stats import spearmanr
from sklearn.metrics import r2_score, mean_squared_error


# ============================================================
# MATRIZ
# ============================================================

df = pd.read_csv(
    "data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv"
)

F1 = [
    "gc_content",
    "mfe_rnafold",
    "guide_length",
    "g_run_max",
    "poly_t_count"
]

y = df["efficacy_score"]
groups = df["nearest_tss_gene"]


# ============================================================
# FOLDS — MISMA CONFIGURACIÓN QUE ABLATION
# ============================================================

gkf = GroupKFold(n_splits=5)
folds = list(gkf.split(df, y, groups))


# ============================================================
# MODELOS
# ============================================================

models = {
    "Ridge": Pipeline([
        ("scaler", StandardScaler()),
        ("model", Ridge(alpha=1.0))
    ]),

    "Random Forest": RandomForestRegressor(
        n_estimators=500,
        max_depth=None,
        min_samples_leaf=5,
        max_features="sqrt",
        bootstrap=True,
        oob_score=True,
        criterion="squared_error",
        random_state=42,
        n_jobs=-1
    )
}


# ============================================================
# BENCHMARK
# ============================================================

results = []

for model_name, model in models.items():

    print("\n========================================")
    print(model_name)
    print("========================================")

    for fold, (train_idx, test_idx) in enumerate(folds, start=1):

        X_train = df.loc[train_idx, F1]
        X_test = df.loc[test_idx, F1]

        y_train = y.loc[train_idx]
        y_test = y.loc[test_idx]

        # Misma imputación que ablation.py
        medians = X_train.median()

        X_train = X_train.fillna(medians)
        X_test = X_test.fillna(medians)

        model.fit(X_train, y_train)

        y_pred = model.predict(X_test)

        spearman = spearmanr(y_test, y_pred).statistic
        r2 = r2_score(y_test, y_pred)
        rmse = np.sqrt(mean_squared_error(y_test, y_pred))

        results.append({
            "model": model_name,
            "fold": fold,
            "spearman": spearman,
            "r2": r2,
            "rmse": rmse
        })

        print(
            f"Fold {fold} | "
            f"Spearman={spearman:.6f} | "
            f"R2={r2:.6f} | "
            f"RMSE={rmse:.6f}"
        )


# ============================================================
# RESULTADOS
# ============================================================

results_df = pd.DataFrame(results)

results_df.to_csv(
    "results/week9_sequence_only_folds.csv",
    index=False
)

summary = (
    results_df
    .groupby("model")[["spearman", "r2", "rmse"]]
    .mean()
    .reset_index()
)

summary.to_csv(
    "results/week9_sequence_only_summary.csv",
    index=False
)

print("\n========================================")
print("SEQUENCE-ONLY SUMMARY")
print("========================================")
print(summary.to_string(index=False))
