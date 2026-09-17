#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer
from sklearn.ensemble import RandomForestRegressor

INPUT = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy.csv")
LOG = Path("logs/week8_feature_importance.md")

RANDOM_STATE = 42

F1 = [
    "gc_content",
    "mfe_rnafold",
    "guide_length",
    "g_run_max",
    "poly_t_count",
]

F2 = [
    "ATAC_mean",
    "ATAC_max",
    "ATAC_p90",
    "ATAC_sum",
]

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

F5 = [
    "nearest_tss_distance",
    "within_promoter_2kb",
]

FEATURES = F1 + F2 + F3 + F5

# 1. Cargar y filtrar subset Horlbeck con efficacy
df = pd.read_csv(INPUT, low_memory=False)

df = df[
    (df["dataset"] == "Horlbeck2016") &
    (df["efficacy_score"].notna())
].copy()

if df.empty:
    raise ValueError("ERROR: no hay filas Horlbeck2016 con efficacy_score.")

# 2. Preparar X/y
missing = [f for f in FEATURES if f not in df.columns]
if missing:
    raise ValueError(f"ERROR: faltan features: {missing}")

X = df[FEATURES].apply(pd.to_numeric, errors="coerce")
y = pd.to_numeric(df["efficacy_score"], errors="coerce")

valid = y.notna()
X = X.loc[valid]
y = y.loc[valid]

# 3. Imputación
imputer = SimpleImputer(strategy="median")
X_imp = pd.DataFrame(
    imputer.fit_transform(X),
    columns=FEATURES,
    index=X.index
)

# 4. Split idéntico al baseline
X_train, X_test, y_train, y_test = train_test_split(
    X_imp,
    y,
    test_size=0.20,
    random_state=RANDOM_STATE
)

# 5. Random Forest
rf = RandomForestRegressor(
    n_estimators=200,
    random_state=RANDOM_STATE,
    n_jobs=-1
)

rf.fit(X_train, y_train)

# 6. Feature importance
importance = pd.DataFrame({
    "feature": FEATURES,
    "importance": rf.feature_importances_
}).sort_values(
    "importance",
    ascending=False
).reset_index(drop=True)

importance.index = importance.index + 1

# 7. Mostrar top 10
print("=" * 80)
print("CHROMACRISPR PHASE 1 — RANDOM FOREST FEATURE IMPORTANCE")
print("=" * 80)
print(f"Subset: Horlbeck2016 con efficacy_score")
print(f"Filas: {len(df):,}")
print(f"Features: {len(FEATURES)}")
print(f"Random state: {RANDOM_STATE}")
print()
print("TOP 10 FEATURES")
print(importance.head(10).to_string(index=True))
print()

# 8. Guardar ranking completo en log
LOG.parent.mkdir(parents=True, exist_ok=True)

with LOG.open("w", encoding="utf-8") as f:
    f.write("# Week 8 — Random Forest Feature Importance\n\n")
    f.write("## Dataset\n\n")
    f.write(f"- Input: `{INPUT}`\n")
    f.write("- Dataset: `Horlbeck2016`\n")
    f.write("- Target: `efficacy_score`\n")
    f.write(f"- Filas utilizadas: {len(df):,}\n")
    f.write(f"- Features: {len(FEATURES)}\n")
    f.write("- Train/test: 80/20\n")
    f.write(f"- Random state: {RANDOM_STATE}\n")
    f.write("- Random Forest: `n_estimators=200`, `n_jobs=-1`\n\n")

    f.write("## Top 10 features\n\n")
    f.write("| Rank | Feature | Importance |\n")
    f.write("|---:|---|---:|\n")

    for rank, row in importance.head(10).iterrows():
        f.write(
            f"| {rank} | `{row['feature']}` | "
            f"{row['importance']:.8f} |\n"
        )

    f.write("\n## Ranking completo\n\n")
    f.write("| Rank | Feature | Importance |\n")
    f.write("|---:|---|---:|\n")

    for rank, row in importance.iterrows():
        f.write(
            f"| {rank} | `{row['feature']}` | "
            f"{row['importance']:.8f} |\n"
        )

print(f"Ranking guardado: {LOG}")
