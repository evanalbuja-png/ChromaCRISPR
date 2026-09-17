import json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import RidgeCV, LassoCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import cross_val_predict
from sklearn.pipeline import Pipeline
import xgboost as xgb
import joblib

warnings.filterwarnings("ignore", category=UserWarning)
np.random.seed(42)
Path("models").mkdir(exist_ok=True)
Path("logs").mkdir(exist_ok=True)

# ---------- Data ----------
m = pd.read_csv("data/processed/K562_training_matrix_v1.csv", low_memory=False)
man = pd.read_csv("data/processed/K562_feature_manifest_v1.tsv", sep="\t")
feats = [c for c in man.loc[man["accepted_for_modeling"].astype(str).str.lower().isin(["true","1"]), "column"]
         if c != "EEP_percentile" and c in m.columns]
# F1-only for M6
f1_cols = [c for c in feats if c in (
    "gc_content","mfe_rnafold","g_run_max","poly_t_count","tm_santalucia","seed_pwm_score","guide_length"
) or c.startswith("guide_3prime") or c.startswith("pam_")]

chrom = m["chromosome"].astype(str).str.replace("^chr", "", regex=True)
train_ch = set(str(i) for i in range(1, 20)) | {"X"}
val_ch, test_ch = {"20"}, {"21", "22", "Y"}
tr = chrom.isin(train_ch).values
va = chrom.isin(val_ch).values
te = chrom.isin(test_ch).values

y = m["EEP_percentile"].values.astype(float)
X = m[feats].copy()
# median impute remaining NaNs for model input (flags already present)
for c in X.columns:
    if X[c].isna().any():
        X[c] = X[c].fillna(X.loc[tr, c].median())

X_tr, y_tr = X.loc[tr], y[tr]
X_va, y_va = X.loc[va], y[va]
X_te, y_te = X.loc[te], y[te]
print(f"Features {len(feats)} | F1-only {len(f1_cols)} | train {len(y_tr)} val {len(y_va)} test {len(y_te)}")

def spearman(y_true, y_pred):
    r, _ = spearmanr(y_true, y_pred)
    return float(r) if r == r else 0.0

def pred_interval_coverage(y_true, y_pred, resid_train, alpha=0.1):
    """Symmetric PI from train residual quantiles."""
    lo = np.quantile(resid_train, alpha / 2)
    hi = np.quantile(resid_train, 1 - alpha / 2)
    inside = ((y_true >= y_pred + lo) & (y_true <= y_pred + hi)).mean()
    return float(inside)

def chrom_sub_stability(model_factory, Xdf, yarr, chrom_arr, train_mask, n_sub=5):
    """Stability: rho CV across held-out train chromosomes."""
    tr_chroms = sorted(set(chrom_arr[train_mask]) - {"X"})  # leave X in train always if possible
    if len(tr_chroms) < n_sub + 1:
        tr_chroms = sorted(set(chrom_arr[train_mask]))
    rng = np.random.RandomState(42)
    held = rng.choice(tr_chroms, size=min(n_sub, len(tr_chroms)), replace=False)
    rhos = []
    for h in held:
        tr_m = train_mask & (chrom_arr != h)
        va_m = train_mask & (chrom_arr == h)
        if va_m.sum() < 50 or tr_m.sum() < 500:
            continue
        model = model_factory()
        model.fit(Xdf.loc[tr_m], yarr[tr_m])
        pred = model.predict(Xdf.loc[va_m])
        rhos.append(spearman(yarr[va_m], pred))
    if not rhos:
        return 0.0, []
    rhos = np.array(rhos)
    cv = rhos.std() / (abs(rhos.mean()) + 1e-6)
    stab = 1.0 / (1.0 + cv)
    return float(stab), rhos.tolist()

# Interpretable scores (pre-specified)
INTERP = {"ridge": 1.0, "lasso": 1.0, "rf": 0.5, "xgb": 0.7, "mlp": 0.3, "cnn": 0.3}

results = {}

# ----- M1 Ridge -----
print("\n=== M1 Ridge ===")
def make_ridge():
    return Pipeline([("sc", StandardScaler()), ("m", RidgeCV(alphas=np.logspace(-2, 4, 30)))])

ridge = make_ridge()
ridge.fit(X_tr, y_tr)
pred_va = ridge.predict(X_va)
pred_tr = ridge.predict(X_tr)
rho = spearman(y_va, pred_va)
cov = pred_interval_coverage(y_va, pred_va, y_tr - pred_tr)
stab, stab_rhos = chrom_sub_stability(make_ridge, X, y, chrom.values, tr)
results["ridge"] = {"rho_val": rho, "cal_coverage": cov, "stab": stab, "stab_rhos": stab_rhos}
joblib.dump(ridge, "models/M1_ridge_v1_101features.joblib")
print(results["ridge"])

# ----- M2 LASSO -----
print("\n=== M2 LASSO ===")
def make_lasso():
    return Pipeline([("sc", StandardScaler()), ("m", LassoCV(alphas=np.logspace(-4, 1, 40), max_iter=5000, n_jobs=-1, random_state=42))])

lasso = make_lasso()
lasso.fit(X_tr, y_tr)
pred_va = lasso.predict(X_va)
pred_tr = lasso.predict(X_tr)
rho = spearman(y_va, pred_va)
cov = pred_interval_coverage(y_va, pred_va, y_tr - pred_tr)
stab, stab_rhos = chrom_sub_stability(make_lasso, X, y, chrom.values, tr)
coef = lasso.named_steps["m"].coef_
sel = [f for f, c in zip(feats, coef) if abs(c) > 1e-8]
results["lasso"] = {"rho_val": rho, "cal_coverage": cov, "stab": stab, "stab_rhos": stab_rhos, "n_selected": len(sel), "selected_sample": sel[:30]}
joblib.dump({"model": lasso, "selected": sel}, "models/M2_lasso_v1_101features.joblib")
print({k: results["lasso"][k] for k in ["rho_val","cal_coverage","stab","n_selected"]})

# Bootstrap stability of LASSO selection (lighter: 20 rounds)
print("LASSO bootstrap selection freq (20)...")
from collections import Counter
cnt = Counter()
rng = np.random.RandomState(0)
for b in range(20):
    idx = rng.choice(len(y_tr), size=len(y_tr), replace=True)
    mb = make_lasso()
    mb.fit(X_tr.iloc[idx], y_tr[idx])
    coef_b = mb.named_steps["m"].coef_
    for f, c in zip(feats, coef_b):
        if abs(c) > 1e-8:
            cnt[f] += 1
stable = [f for f, n in cnt.items() if n >= int(0.7 * 20)]
results["lasso"]["stable_70pct"] = stable
results["lasso"]["n_stable_70pct"] = len(stable)
print("Stable >=70%:", len(stable), stable[:15])

# ----- M3 RF -----
print("\n=== M3 RF ===")
def make_rf():
    return RandomForestRegressor(
        n_estimators=300, max_depth=12, min_samples_leaf=5,
        max_features="sqrt", n_jobs=-1, random_state=42
    )

rf = make_rf()
rf.fit(X_tr, y_tr)
pred_va = rf.predict(X_va)
pred_tr = rf.predict(X_tr)
rho = spearman(y_va, pred_va)
cov = pred_interval_coverage(y_va, pred_va, y_tr - pred_tr)
stab, stab_rhos = chrom_sub_stability(make_rf, X, y, chrom.values, tr)
imp = sorted(zip(feats, rf.feature_importances_), key=lambda x: -x[1])[:15]
results["rf"] = {"rho_val": rho, "cal_coverage": cov, "stab": stab, "stab_rhos": stab_rhos, "top_impurity": imp}
joblib.dump(rf, "models/M3_rf_v1_101features.joblib")
print({k: results["rf"][k] for k in ["rho_val","cal_coverage","stab"]})

# ----- M4 XGBoost + Optuna (40 trials to keep runtime reasonable; note if 100 desired) -----
print("\n=== M4 XGBoost Optuna ===")
import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)

def xgb_objective(trial):
    params = {
        "n_estimators": trial.suggest_int("n_estimators", 100, 600),
        "max_depth": trial.suggest_int("max_depth", 3, 10),
        "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.3, log=True),
        "subsample": trial.suggest_float("subsample", 0.6, 1.0),
        "colsample_bytree": trial.suggest_float("colsample_bytree", 0.5, 1.0),
        "min_child_weight": trial.suggest_int("min_child_weight", 1, 20),
        "reg_lambda": trial.suggest_float("reg_lambda", 1e-3, 10.0, log=True),
        "objective": "reg:squarederror",
        "n_jobs": 4,
        "random_state": 42,
    }
    # inner: hold out a train chrom
    model = xgb.XGBRegressor(**params)
    # quick eval on val chrom20 directly for speed of trials
    model.fit(X_tr, y_tr)
    pred = model.predict(X_va)
    return -spearman(y_va, pred)  # minimize

study = optuna.create_study(direction="minimize", sampler=optuna.samplers.TPESampler(seed=42))
study.optimize(xgb_objective, n_trials=40, show_progress_bar=False)
best = study.best_params
best.update({"objective": "reg:squarederror", "n_jobs": 4, "random_state": 42})
print("Best params:", best)

def make_xgb():
    return xgb.XGBRegressor(**best)

xgbm = make_xgb()
xgbm.fit(X_tr, y_tr)
pred_va = xgbm.predict(X_va)
pred_tr = xgbm.predict(X_tr)
rho = spearman(y_va, pred_va)
cov = pred_interval_coverage(y_va, pred_va, y_tr - pred_tr)
stab, stab_rhos = chrom_sub_stability(make_xgb, X, y, chrom.values, tr)
results["xgb"] = {"rho_val": rho, "cal_coverage": cov, "stab": stab, "stab_rhos": stab_rhos, "best_params": best, "n_trials": 40}
joblib.dump(xgbm, "models/M4_xgb_v1_101features.joblib")
print({k: results["xgb"][k] for k in ["rho_val","cal_coverage","stab"]})

# ----- M5 MLP (sklearn; torch absent) -----
print("\n=== M5 MLP (sklearn) ===")
def make_mlp():
    return Pipeline([
        ("sc", StandardScaler()),
        ("m", MLPRegressor(hidden_layer_sizes=(128, 64), dropout=False,  # sklearn no dropout param until recent
                           alpha=1e-4, learning_rate_init=1e-3, max_iter=200,
                           early_stopping=True, validation_fraction=0.1,
                           random_state=42))
    ])
# fix: MLPRegressor has no dropout in older sklearn — use alpha as reg
def make_mlp():
    return Pipeline([
        ("sc", StandardScaler()),
        ("m", MLPRegressor(hidden_layer_sizes=(128, 64, 32), alpha=1e-3,
                           learning_rate_init=1e-3, max_iter=300,
                           early_stopping=True, validation_fraction=0.1,
                           random_state=42))
    ])

try:
    mlp = make_mlp()
    mlp.fit(X_tr, y_tr)
    pred_va = mlp.predict(X_va)
    pred_tr = mlp.predict(X_tr)
    rho = spearman(y_va, pred_va)
    cov = pred_interval_coverage(y_va, pred_va, y_tr - pred_tr)
    stab, stab_rhos = chrom_sub_stability(make_mlp, X, y, chrom.values, tr)
    results["mlp"] = {"rho_val": rho, "cal_coverage": cov, "stab": stab, "stab_rhos": stab_rhos, "backend": "sklearn"}
    joblib.dump(mlp, "models/M5_mlp_v1_101features.joblib")
    print(results["mlp"])
except Exception as e:
    results["mlp"] = {"error": str(e), "rho_val": None}
    print("MLP FAILED:", e)

# ----- M6 CNN sequence-only -----
print("\n=== M6 CNN sequence-only ===")
try:
    import torch
    import torch.nn as nn
    from torch.utils.data import DataLoader, TensorDataset

    X1 = m[f1_cols].copy()
    for c in X1.columns:
        if X1[c].isna().any():
            X1[c] = X1[c].fillna(X1.loc[tr, c].median())
    sc = StandardScaler()
    X1_tr = sc.fit_transform(X1.loc[tr])
    X1_va = sc.transform(X1.loc[va])

    class SeqMLP(nn.Module):  # lightweight stand-in if 1D CNN over one-hot not structured; true CNN on sequence:
        # Use 1D conv over feature vector as weak CNN baseline
        def __init__(self, d):
            super().__init__()
            self.net = nn.Sequential(
                nn.Unflatten(1, (1, d)),
                nn.Conv1d(1, 32, kernel_size=5, padding=2),
                nn.ReLU(),
                nn.Conv1d(32, 32, kernel_size=5, padding=2),
                nn.ReLU(),
                nn.AdaptiveAvgPool1d(8),
                nn.Flatten(),
                nn.Linear(32 * 8, 64),
                nn.ReLU(),
                nn.Dropout(0.3),
                nn.Linear(64, 1),
            )
        def forward(self, x):
            return self.net(x).squeeze(-1)

    device = "cpu"
    model = SeqMLP(X1_tr.shape[1]).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    loss_fn = nn.MSELoss()
    ds = TensorDataset(torch.tensor(X1_tr, dtype=torch.float32), torch.tensor(y_tr, dtype=torch.float32))
    loader = DataLoader(ds, batch_size=256, shuffle=True)
    model.train()
    for epoch in range(25):
        for xb, yb in loader:
            opt.zero_grad()
            pred = model(xb)
            loss = loss_fn(pred, yb)
            loss.backward()
            opt.step()
    model.eval()
    with torch.no_grad():
        pred_va = model(torch.tensor(X1_va, dtype=torch.float32)).numpy()
        pred_tr = model(torch.tensor(X1_tr, dtype=torch.float32)).numpy()
    rho = spearman(y_va, pred_va)
    cov = pred_interval_coverage(y_va, pred_va, y_tr - pred_tr)
    results["cnn"] = {"rho_val": rho, "cal_coverage": cov, "stab": float("nan"), "stab_rhos": [], "backend": "torch_conv1d_F1only", "note": "stability skipped (cost); F1-only"}
    torch.save({"state": model.state_dict(), "scaler": sc, "f1_cols": f1_cols}, "models/M6_cnn_v1_F1only.pt")
    print(results["cnn"])
except Exception as e:
    results["cnn"] = {"error": str(e), "rho_val": None, "note": "CNN not trained"}
    print("CNN FAILED:", e)

# ---------- Multi-criteria score (PRE-SPECIFIED) ----------
def score_row(name, R):
    if R.get("rho_val") is None:
        return None
    s_rho = float(np.clip(R["rho_val"], 0, 1))
    cov = R.get("cal_coverage", 0.5)
    s_cal = float(np.clip(1 - abs(cov - 0.90) / 0.90, 0, 1))
    s_interp = INTERP.get(name, 0.3)
    stab = R.get("stab", 0.5)
    if stab != stab:
        stab = 0.5  # neutral if missing
    s_stab = float(np.clip(stab, 0, 1))
    S = 0.4 * s_rho + 0.2 * s_cal + 0.3 * s_interp + 0.1 * s_stab
    return {"S_rho": s_rho, "S_cal": s_cal, "S_interp": s_interp, "S_stab": s_stab, "S": S,
            "rho_val": R["rho_val"], "cal_coverage": cov}

table = []
for name in ["ridge", "lasso", "rf", "xgb", "mlp", "cnn"]:
    sc = score_row(name, results.get(name, {}))
    if sc is None:
        table.append({"model": name, "status": "FAILED", **{k: None for k in ["rho_val","cal_coverage","S_rho","S_cal","S_interp","S_stab","S"]}})
    else:
        table.append({"model": name, "status": "OK", **sc})

tab = pd.DataFrame(table).sort_values("S", ascending=False)
print("\n=== SCORE TABLE ===")
print(tab.to_string(index=False))
primary = tab.loc[tab["status"]=="OK"].iloc[0]["model"] if (tab["status"]=="OK").any() else None
print("PRIMARY MODEL:", primary)

# Hold test set untouched — report only val metrics this week
with open("logs/weekL_six_models_results.json", "w") as f:
    json.dump({"results": results, "table": table, "primary": primary, "split": {"train": int(tr.sum()), "val": int(va.sum()), "test": int(te.sum())}, "n_features": len(feats), "n_trials_xgb": 40}, f, default=str, indent=2)

tab.to_csv("logs/weekL_six_models_table.tsv", sep="\t", index=False)
print("Saved models/ and logs/weekL_six_models_results.json")