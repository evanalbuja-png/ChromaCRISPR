#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge, Lasso
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import spearmanr, pearsonr

INPUT = "data/processed/sgRNA_feature_matrix_phase1.csv"
OUTPUT = "logs/week8_baseline_models.md"

TARGET_CANDIDATES = [
    "efficacy",
    "efficacy_score",
    "efficacy_score_norm",
    "score_norm",
    "efficacy_norm",
    "activity",
    "activity_score",
    "on_target_score",
]

print("=" * 70)
print("ChromaCRISPR Phase 1 — Week 8")
print("Paso 4 — Baseline Models")
print("=" * 70)

print("\nLoading matrix...")
df = pd.read_csv(INPUT, low_memory=False)

print(f"Rows: {len(df):,}")
print(f"Columns: {len(df.columns)}")

# ------------------------------------------------------------
# 1. Detect target
# ------------------------------------------------------------

print("\n=== TARGET SEARCH ===")

target = None

for candidate in TARGET_CANDIDATES:
    if candidate in df.columns:
        target = candidate
        break

if target is None:
    print("\nNO TARGET VARIABLE FOUND.")
    print("\nAvailable columns:")
    for i, col in enumerate(df.columns):
        print(f"{i}: {col}")

    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)

    with open(OUTPUT, "w") as f:
        f.write("# ChromaCRISPR Phase 1 — Week 8\n\n")
        f.write("## Paso 4 — Baseline Models\n\n")
        f.write("### Estado: DETENIDO\n\n")
        f.write(
            "No se encontró una variable objetivo de eficacia/score "
            "en la matriz final.\n\n"
        )
        f.write("Columnas disponibles:\n\n")
        for col in df.columns:
            f.write(f"- `{col}`\n")

    print(f"\nReport saved: {OUTPUT}")
    print("\nPaso 4 detenido: no existe variable objetivo.")
    sys.exit(0)

print(f"Target detected: {target}")

# ------------------------------------------------------------
# 2. Select numeric F1/F2/F3/F5 features
# ------------------------------------------------------------

EXCLUDE = {
    target,
    "sgrna_id",
    "region_id",
    "dataset",
    "experiment",
    "guide_sequence",
    "gene_symbol",
    "gene_id",
    "chromosome",
    "coordinate",
    "strand",
}

numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

features = [
    col for col in numeric_cols
    if col not in EXCLUDE
]

print("\n=== FEATURES ===")
print(f"Numeric features selected: {len(features)}")

for col in features:
    print(f"- {col}")

if len(features) == 0:
    raise RuntimeError("No numeric F1/F2/F3/F5 features available.")

# ------------------------------------------------------------
# 3. Prepare data
# ------------------------------------------------------------

data = df[features + [target]].copy()

before = len(data)
data = data.dropna(subset=[target])
after = len(data)

print("\n=== TARGET QUALITY ===")
print(f"Rows before target filtering: {before:,}")
print(f"Rows after target filtering:  {after:,}")
print(f"Target missing:               {before - after:,}")

if after < 100:
    raise RuntimeError("Too few observations after removing missing target values.")

X = data[features]
y = data[target]

print(f"\nTarget mean: {y.mean():.6f}")
print(f"Target std:  {y.std():.6f}")
print(f"Target min:  {y.min():.6f}")
print(f"Target max:  {y.max():.6f}")

# ------------------------------------------------------------
# 4. Train/test split
# ------------------------------------------------------------

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.20,
    random_state=42,
)

print("\n=== TRAIN / TEST ===")
print(f"Train: {len(X_train):,}")
print(f"Test:  {len(X_test):,}")

# ------------------------------------------------------------
# 5. Models
# ------------------------------------------------------------

models = {
    "Ridge": Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("model", Ridge(alpha=1.0)),
    ]),

    "LASSO": Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("model", Lasso(alpha=0.001, max_iter=10000)),
    ]),

    "Random Forest": Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        (
            "model",
            RandomForestRegressor(
                n_estimators=200,
                max_depth=None,
                min_samples_leaf=2,
                random_state=42,
                n_jobs=-1,
            ),
        ),
    ]),
}

results = []

print("\n=== TRAINING ===")

for name, model in models.items():

    print(f"\nTraining {name}...")

    model.fit(X_train, y_train)
    pred = model.predict(X_test)

    rmse = np.sqrt(mean_squared_error(y_test, pred))
    r2 = r2_score(y_test, pred)
    pearson_r, _ = pearsonr(y_test, pred)
    spearman_rho, _ = spearmanr(y_test, pred)

    results.append({
        "Model": name,
        "Spearman rho": spearman_rho,
        "Pearson R2": r2,
        "Pearson R": pearson_r,
        "RMSE": rmse,
    })

    print(f"Spearman rho: {spearman_rho:.4f}")
    print(f"Pearson R:    {pearson_r:.4f}")
    print(f"Pearson R²:   {r2:.4f}")
    print(f"RMSE:         {rmse:.4f}")

results_df = pd.DataFrame(results)

# ------------------------------------------------------------
# 6. Save report
# ------------------------------------------------------------

os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)

with open(OUTPUT, "w") as f:

    f.write("# ChromaCRISPR Phase 1 — Week 8\n\n")
    f.write("## Paso 4 — Baseline Models\n\n")

    f.write("### Dataset\n\n")
    f.write(f"- Input: `{INPUT}`\n")
    f.write(f"- Total rows: {len(df):,}\n")
    f.write(f"- Target: `{target}`\n")
    f.write(f"- Rows with target: {len(data):,}\n")
    f.write(f"- Train: {len(X_train):,}\n")
    f.write(f"- Test: {len(X_test):,}\n")
    f.write("- Split: 80/20\n")
    f.write("- Random seed: 42\n\n")

    f.write("### Numeric features used\n\n")
    for col in features:
        f.write(f"- `{col}`\n")

    f.write("\n### Baseline models\n\n")
    f.write("- Ridge Regression (`alpha=1.0`)\n")
    f.write("- LASSO (`alpha=0.001`)\n")
    f.write("- Random Forest (`n_estimators=200`, `min_samples_leaf=2`)\n\n")

    f.write("### Test-set metrics\n\n")
    f.write(
        "| Model | Spearman rho | Pearson R | Pearson R² | RMSE |\n"
    )
    f.write(
        "|---|---:|---:|---:|---:|\n"
    )

    for _, row in results_df.iterrows():
        f.write(
            f"| {row['Model']} | "
            f"{row['Spearman rho']:.4f} | "
            f"{row['Pearson R']:.4f} | "
            f"{row['Pearson R2']:.4f} | "
            f"{row['RMSE']:.4f} |\n"
        )

    f.write("\n### Notes\n\n")
    f.write(
        "- Missing numeric feature values were median-imputed within "
        "the training pipeline.\n"
    )
    f.write(
        "- Ridge and LASSO used standardized features.\n"
    )
    f.write(
        "- No feature selection or hyperparameter optimization was performed.\n"
    )
    f.write(
        "- This is an initial baseline evaluation, not the final model.\n"
    )

print("\n=== REPORT ===")
print(OUTPUT)
print("\nBaseline analysis completed.")

