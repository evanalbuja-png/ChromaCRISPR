import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import RidgeCV
import xgboost as xgb
import joblib
import json
from pathlib import Path

m = pd.read_csv("data/processed/K562_training_matrix_v1.csv", low_memory=False)
man = pd.read_csv("data/processed/K562_feature_manifest_v1.tsv", sep="\t")
feats = [c for c in man.loc[man["accepted_for_modeling"].astype(str).str.lower().isin(["true","1"]),"column"]
         if c != "EEP_percentile" and c in m.columns]

chrom_raw = m["chromosome"].astype(str)
chrom = chrom_raw.str.replace("^chr", "", regex=True)
# 23 folds: 1-22 + X (exclude Y if tiny / optional — Sec often includes autosomes+X)
chroms = [str(i) for i in range(1, 23)] + ["X"]
print("Available chrom counts:")
print(chrom.value_counts().reindex(chroms + ["Y"]).to_string())

y = m["EEP_percentile"].values.astype(float)
X = m[feats].copy()
# global median impute from full data for LOCO simplicity (or per-fold train median)
for c in X.columns:
    if X[c].isna().any():
        X[c] = X[c].fillna(X[c].median())

# Hyperparams from Semana L (fixed, no retune)
# Ridge: RidgeCV alphas same grid as L
# XGB: best from 100 trials
xgb_params = {
    "n_estimators": 289, "max_depth": 4, "learning_rate": 0.02972886488008395,
    "subsample": 0.7663021691191627, "colsample_bytree": 0.5668763022299779,
    "min_child_weight": 4, "reg_lambda": 0.31785318816701275,
    "objective": "reg:squarederror", "n_jobs": 4, "random_state": 42,
}

def spearman(a, b):
    r, p = spearmanr(a, b)
    return (float(r) if r == r else np.nan), (float(p) if p == p else np.nan)

rows = []
for h in chroms:
    te = (chrom == h).values
    tr = ~te
    n_te = int(te.sum())
    if n_te < 20:
        rows.append({"chrom": h, "n_test": n_te, "ridge_rho": np.nan, "xgb_rho": np.nan, "note": "too few"})
        print(f"chr{h}: skip n={n_te}")
        continue
    Xtr, ytr = X.loc[tr], y[tr]
    Xte, yte = X.loc[te], y[te]
    # per-fold median already done globally; OK

    ridge = Pipeline([("sc", StandardScaler()),
                      ("m", RidgeCV(alphas=np.logspace(-2, 4, 30)))])
    ridge.fit(Xtr, ytr)
    pr = ridge.predict(Xte)
    rr, rp = spearman(yte, pr)

    xgbm = xgb.XGBRegressor(**xgb_params)
    xgbm.fit(Xtr, ytr)
    px = xgbm.predict(Xte)
    xr, xp = spearman(yte, px)

    rows.append({"chrom": h, "n_test": n_te, "ridge_rho": rr, "ridge_p": rp,
                 "xgb_rho": xr, "xgb_p": xp, "note": ""})
    print(f"chr{h:>2}: n={n_te:5d}  ridge={rr:.4f}  xgb={xr:.4f}")

df = pd.DataFrame(rows)
print("\n=== LOCO SUMMARY ===")
for model in ["ridge", "xgb"]:
    col = f"{model}_rho"
    valid = df[col].dropna()
    print(f"{model}: mean={valid.mean():.4f}  sd={valid.std(ddof=1):.4f}  "
          f"median={valid.median():.4f}  min={valid.min():.4f}  max={valid.max():.4f}  n_folds={len(valid)}")

# Outliers: >2 SD from mean
for model in ["ridge", "xgb"]:
    col = f"{model}_rho"
    mu, sd = df[col].mean(), df[col].std(ddof=1)
    df[f"{model}_outlier"] = (df[col] - mu).abs() > 2 * sd
    print(f"\n{model} outliers (|z|>2):")
    print(df.loc[df[f'{model}_outlier'], ["chrom","n_test",col]].to_string(index=False))

# Correlation of rho with N
for model in ["ridge", "xgb"]:
    r, p = spearmanr(df["n_test"], df[f"{model}_rho"], nan_policy="omit")
    print(f"Spearman(n_test, {model}_rho) = {r:.3f} p={p:.3f}")

df.to_csv("logs/weekM_loco_folds.tsv", sep="\t", index=False)
summary = {
    "ridge": {"mean": float(df["ridge_rho"].mean()), "sd": float(df["ridge_rho"].std(ddof=1)),
              "median": float(df["ridge_rho"].median())},
    "xgb": {"mean": float(df["xgb_rho"].mean()), "sd": float(df["xgb_rho"].std(ddof=1)),
            "median": float(df["xgb_rho"].median())},
    "ref_val_chr20_ridge": 0.335,
    "ref_val_chr20_xgb": 0.403,
}
Path("logs/weekM_loco_summary.json").write_text(json.dumps(summary, indent=2))
print("\nSaved logs/weekM_loco_folds.tsv")
print(summary)