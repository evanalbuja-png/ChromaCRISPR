import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import xgboost as xgb
from scipy.stats import spearmanr
import joblib

Path("logs/weekP_plots").mkdir(parents=True, exist_ok=True)

man = pd.read_csv("data/processed/K562_feature_manifest_v1.tsv", sep="\t")
feats = [c for c in man.loc[man["accepted_for_modeling"].astype(str).str.lower().isin(["true","1"]),"column"]
         if c != "EEP_percentile"]

m = pd.read_csv("data/processed/K562_training_matrix_v1.csv", low_memory=False)
chrom = m["chromosome"].astype(str).str.replace("^chr","", regex=True)
tr = chrom.isin(set(str(i) for i in range(1,20))|{"X"}).values
va = chrom.isin({"20"}).values
te = chrom.isin({"21","22","Y"}).values

X = m[feats].copy()
for c in X.columns:
    if X[c].isna().any():
        X[c] = X[c].fillna(X.loc[tr, c].median())
y = m["EEP_percentile"].values

params = dict(
    n_estimators=289, max_depth=4, learning_rate=0.02972886488008395,
    subsample=0.7663021691191627, colsample_bytree=0.5668763022299779,
    min_child_weight=4, reg_lambda=0.31785318816701275,
    objective="reg:squarederror", n_jobs=4, random_state=42,
    base_score=0.5,  # avoid [5E1] serialization issue if retraining
)
model = xgb.XGBRegressor(**params)
model.fit(X.loc[tr], y[tr])
print("val rho", spearmanr(y[va], model.predict(X.loc[va])).correlation)
print("test rho", spearmanr(y[te], model.predict(X.loc[te])).correlation)
joblib.dump(model, "models/M4_xgb_EEP_train_for_SHAP.joblib")

X_te = X.loc[te]
y_te = y[te]
# TreeSHAP via XGBoost pred_contribs (last column = bias)
contrib = model.get_booster().predict(xgb.DMatrix(X_te), pred_contribs=True)
shap_vals = contrib[:, :-1]
bias = contrib[:, -1]
print("shap shape", shap_vals.shape, "bias mean", bias.mean())

mean_abs = np.abs(shap_vals).mean(axis=0)
imp = pd.DataFrame({"feature": feats, "mean_abs_shap": mean_abs})
# bootstrap CI
rng = np.random.RandomState(42)
n = shap_vals.shape[0]
boot = np.zeros((200, len(feats)))
for b in range(200):
    idx = rng.choice(n, size=n, replace=True)
    boot[b] = np.abs(shap_vals[idx]).mean(axis=0)
imp["boot_lo"] = np.percentile(boot, 2.5, axis=0)
imp["boot_hi"] = np.percentile(boot, 97.5, axis=0)
imp = imp.sort_values("mean_abs_shap", ascending=False)
print("\nTOP 15:")
print(imp.head(15).to_string(index=False))
imp.to_csv("logs/weekP_shap_global_importance.tsv", sep="\t", index=False)

def modality(f):
    if f.startswith(("ATAC","dnase","atac_")): return "F2"
    if f.startswith(("H3K","CTCF","chromhmm")): return "F3"
    if f.startswith("hic_"): return "F4"
    if f.startswith(("log_dist","ccre","phylop","dist_")): return "F5"
    if "missing" in f: return "flag"
    return "F1"
imp["mod"] = imp["feature"].map(modality)
print("\nSum mean|SHAP| by modality:")
print(imp.groupby("mod")["mean_abs_shap"].sum().sort_values(ascending=False))

# plots top 20
top = imp.head(20)
fig, ax = plt.subplots(figsize=(8,6))
ax.barh(list(top["feature"])[::-1], list(top["mean_abs_shap"])[::-1], color="steelblue")
ax.set_xlabel("mean |SHAP|"); ax.set_title("Global SHAP (test chr21/22)")
plt.tight_layout(); fig.savefig("logs/weekP_plots/global_importance.png", dpi=120); plt.close()

for feat in imp["feature"].head(5):
    j = feats.index(feat)
    fig, ax = plt.subplots(figsize=(5,4))
    ax.scatter(X_te[feat].values, shap_vals[:,j], s=8, alpha=0.5, c="teal")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel(feat); ax.set_ylabel("SHAP")
    plt.tight_layout()
    fig.savefig(f"logs/weekP_plots/dependence_{feat.replace('/','_')}.png", dpi=120)
    plt.close()

# Waterfall high/low
pred_te = model.predict(X_te)
order = np.argsort(pred_te)
for label, i in [("low", order[0]), ("high", order[-1])]:
    sv = shap_vals[i]
    topj = np.argsort(np.abs(sv))[::-1][:8]
    print(f"\n{label} pred={pred_te[i]:.2f} true={y_te[i]:.2f}")
    for j in topj:
        print(f"  {feats[j]}: x={X_te.iloc[i][feats[j]]:.4g} shap={sv[j]:+.3f}")

# I1 top EEP in test
te_df = m.loc[te].copy().reset_index(drop=True)
te_df["pred"] = pred_te
print("\n=== I1 top 5 EEP in test ===")
top5 = te_df.nlargest(5, "EEP_percentile")
for i, row in top5.iterrows():
    # i is position after reset_index if we use nlargest on te_df with default index 0..n-1
    pass
# nlargest keeps original index from te_df after reset
Xte_pos = {idx: i for i, idx in enumerate(te_df.index)}
for idx, row in top5.iterrows():
    pos = list(te_df.index).index(idx)
    bits=[]
    for c in ["ATAC_pm100_mean","ATAC_pm500_mean","H3K27ac_pm500_mean","log_dist_target_tss","gc_content","poly_t_count"]:
        if c in feats:
            bits.append(f"{c}={row.get(c,np.nan):.3g}/sh={shap_vals[pos,feats.index(c)]:+.2f}")
    print(f"{row['gene_target']} EEP={row['EEP_percentile']:.1f} pred={row['pred']:.1f} | " + "; ".join(bits))

# I2 high H3K9me3
if "H3K9me3_pm500_mean" in feats:
    j = feats.index("H3K9me3_pm500_mean")
    print("\n=== I2 top 5 H3K9me3 in test ===")
    hi = te_df.nlargest(5, "H3K9me3_pm500_mean")
    for idx, row in hi.iterrows():
        pos = list(te_df.index).index(idx)
        print(f"{row['gene_target']} EEP={row['EEP_percentile']:.1f} H3K9me3={row['H3K9me3_pm500_mean']:.3f} SHAP={shap_vals[pos,j]:+.3f}")

print("OK")