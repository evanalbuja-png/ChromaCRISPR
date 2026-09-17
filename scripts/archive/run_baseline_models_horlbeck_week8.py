#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge, Lasso
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import spearmanr

INPUT = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy.csv")
LOG = Path("logs/week8_baseline_models_horlbeck.md")

RANDOM_STATE = 42
TEST_SIZE = 0.20

# -------------------------------------------------------------------------
# 1. CARGAR Y FILTRAR HORLBECK CON EFFICACY
# -------------------------------------------------------------------------

df = pd.read_csv(INPUT, low_memory=False)

df = df[df["efficacy_score"].notna()].copy()

if len(df) == 0:
    raise ValueError("ERROR: no hay filas con efficacy_score.")

# -------------------------------------------------------------------------
# 2. SELECCIONAR FEATURES NUMÉRICAS F1 + F2 + F3 + F5
# -------------------------------------------------------------------------

# F1: sequence features
F1 = [
    "gc_content",
    "mfe_rnafold",
    "guide_length",
    "g_run_max",
    "poly_t_count",
]

# F2: ATAC
F2 = [
    "ATAC_mean",
    "ATAC_max",
    "ATAC_p90",
    "ATAC_sum",
]

# F3: histone marks
F3 = [
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
]

# F5: genomic context
F5 = [
    "nearest_tss_distance",
    "within_promoter_2kb",
]

FEATURES = F1 + F2 + F3 + F5

missing_features = [f for f in FEATURES if f not in df.columns]

if missing_features:
    raise ValueError(
        f"ERROR: faltan features requeridas: {missing_features}"
    )

# Convertir exclusivamente las features a numérico
X = df[FEATURES].apply(pd.to_numeric, errors="coerce")
y = pd.to_numeric(df["efficacy_score"], errors="coerce")

valid_y = y.notna()

X = X.loc[valid_y].copy()
y = y.loc[valid_y].copy()

# -------------------------------------------------------------------------
# 3. TRAIN / TEST
# -------------------------------------------------------------------------

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=TEST_SIZE,
    random_state=RANDOM_STATE,
)

# -------------------------------------------------------------------------
# 4. MODELOS
# -------------------------------------------------------------------------

models = {
    "Ridge Regression": Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("model", Ridge(alpha=1.0)),
    ]),

    "LASSO": Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("model", Lasso(alpha=0.01, max_iter=10000)),
    ]),

    "Random Forest": Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("model", RandomForestRegressor(
            n_estimators=200,
            random_state=RANDOM_STATE,
            n_jobs=-1,
        )),
    ]),
}

# -------------------------------------------------------------------------
# 5. ENTRENAMIENTO + EVALUACIÓN
# -------------------------------------------------------------------------

results = []

for name, model in models.items():

    model.fit(X_train, y_train)

    pred = model.predict(X_test)

    spearman_rho, spearman_p = spearmanr(y_test, pred)
    pearson_r2 = r2_score(y_test, pred)
    rmse = np.sqrt(mean_squared_error(y_test, pred))

    results.append({
        "model": name,
        "spearman_rho": spearman_rho,
        "pearson_R2": pearson_r2,
        "RMSE": rmse,
    })

# -------------------------------------------------------------------------
# 6. REPORTE
# -------------------------------------------------------------------------

results_df = pd.DataFrame(results)

print("=" * 80)
print("CHROMACRISPR PHASE 1 — BASELINE MODELS — HORLBECK2016")
print("=" * 80)
print(f"Dataset: {INPUT}")
print(f"Rows with efficacy_score : {len(df):,}")
print(f"Features                 : {len(FEATURES)}")
print(f"Train rows               : {len(X_train):,}")
print(f"Test rows                : {len(X_test):,}")
print(f"Random state             : {RANDOM_STATE}")
print()

print(results_df.to_string(index=False))
print()

# -------------------------------------------------------------------------
# 7. GUARDAR LOG
# -------------------------------------------------------------------------

LOG.parent.mkdir(parents=True, exist_ok=True)

with LOG.open("w", encoding="utf-8") as f:

    f.write("# Week 8 — Baseline Models sobre Horlbeck2016\n\n")

    f.write("## Dataset\n\n")
    f.write(f"- Input: `{INPUT}`\n")
    f.write(f"- Filas con `efficacy_score`: {len(df):,}\n")
    f.write(f"- Features F1 + F2 + F3 + F5: {len(FEATURES)}\n")
    f.write(f"- Train/test: 80/20\n")
    f.write(f"- Random state: {RANDOM_STATE}\n")
    f.write(f"- Train rows: {len(X_train):,}\n")
    f.write(f"- Test rows: {len(X_test):,}\n\n")

    f.write("## Features utilizadas\n\n")

    for feature in FEATURES:
        f.write(f"- `{feature}`\n")

    f.write("\n## Modelos baseline\n\n")

    f.write(
        "| Model | Spearman rho | Pearson R² | RMSE |\n"
        "|---|---:|---:|---:|\n"
    )

    for _, row in results_df.iterrows():
        f.write(
            f"| {row['model']} | "
            f"{row['spearman_rho']:.6f} | "
            f"{row['pearson_R2']:.6f} | "
            f"{row['RMSE']:.6f} |\n"
        )

    f.write("\n## Configuración\n\n")
    f.write("- Ridge: `alpha=1.0`, StandardScaler + median imputation.\n")
    f.write("- LASSO: `alpha=0.01`, StandardScaler + median imputation.\n")
    f.write("- Random Forest: `n_estimators=200`, `random_state=42`, `n_jobs=-1`.\n")
    f.write("- Métrica de correlación: Spearman rho.\n")
    f.write("- Métrica de ajuste: R² de Pearson mediante `r2_score`.\n")
    f.write("- Error: RMSE.\n")

print(f"Reporte guardado: {LOG}")
