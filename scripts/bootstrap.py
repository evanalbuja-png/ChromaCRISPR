import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
from xgboost import XGBRegressor

ROOT = Path(".")
RNG = np.random.default_rng(42)
N_BOOT = 1000

def boot_spearman(y, p, n_boot=N_BOOT, seed=42):
    rng = np.random.default_rng(seed)
    y, p = np.asarray(y, float), np.asarray(p, float)
    m = np.isfinite(y) & np.isfinite(p)
    y, p = y[m], p[m]
    n = len(y)
    rhos = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        rhos[i] = spearmanr(y[idx], p[idx]).correlation
    point = float(spearmanr(y, p).correlation)
    lo, hi = np.percentile(rhos, [2.5, 97.5])
    return point, float(lo), float(hi), n, rhos

def boot_mean_of_values(values, n_boot=N_BOOT, seed=42):
    """Bootstrap CI for mean of fold-level rhos (n folds small)."""
    rng = np.random.default_rng(seed)
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    means = np.empty(n_boot)
    for i in range(n_boot):
        means[i] = rng.choice(v, size=len(v), replace=True).mean()
    return float(v.mean()), float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5)), len(v)

def fmt(point, lo, hi):
    return f"{point:.3f} [{lo:.3f}, {hi:.3f}]"

params = dict(
    max_depth=4, learning_rate=0.02972886488008395, n_estimators=289,
    min_child_weight=4, subsample=0.7663021691191627,
    colsample_bytree=0.5668763022299779, reg_lambda=0.31785318816701275,
    objective="reg:squarederror", n_jobs=-1, random_state=42, tree_method="hist",
)
feats = json.loads((ROOT/"models/M4_feature_list_101.json").read_text())
mat = pd.read_csv(ROOT/"data/processed/K562_training_matrix_v1.csv", low_memory=False)
chrom = mat["chromosome"].astype(str).str.replace("^chr", "", regex=True).values
X = mat[feats].apply(pd.to_numeric, errors="coerce").fillna(0).values
folds = [str(i) for i in range(1, 23)] + ["X"]

# ========== 1) LOCO EEP OOF (regenerate if needed) ==========
eep_oof_path = ROOT/"results/xgb_oof_predictions_eep_m4.csv"
if eep_oof_path.exists():
    eep_oof = pd.read_csv(eep_oof_path)
    print("loaded", eep_oof_path)
else:
    print("Regenerating EEP LOCO OOF (one-time)...")
    y = mat["EEP_percentile"].astype(float).values
    oof = np.full(len(mat), np.nan)
    fold_rhos = []
    for h in folds:
        te = chrom == h
        if te.sum() < 20:
            continue
        model = XGBRegressor(**params)
        model.fit(X[~te], y[~te])
        pred = model.predict(X[te])
        oof[te] = pred
        fold_rhos.append(spearmanr(y[te], pred).correlation)
        print(f"  EEP chr{h}: {fold_rhos[-1]:.4f}")
    eep_oof = mat[["guide_sequence","chromosome"]].copy()
    eep_oof["EEP_percentile"] = y
    eep_oof["xgb_oof_pred_eep"] = oof
    eep_oof.to_csv(eep_oof_path, index=False)
    pd.DataFrame({"fold": [f for f in folds if (chrom==f).sum()>=20], "rho": fold_rhos}).to_csv(
        ROOT/"results/weekU_eep_loco_folds.csv", index=False)
    print("WROTE", eep_oof_path, "mean fold", np.mean(fold_rhos))

# z OOF already exists
z_oof = pd.read_csv(ROOT/"results/xgb_oof_predictions_zscore_continuous.csv")

# Bootstrap on OOF pairs
eep_pt, eep_lo, eep_hi, eep_n, _ = boot_spearman(
    eep_oof["EEP_percentile"], eep_oof["xgb_oof_pred_eep"], seed=1)
z_pt, z_lo, z_hi, z_n, _ = boot_spearman(
    z_oof["z_score_within_study"], z_oof["xgb_oof_pred_z"], seed=2)

# Also CI on mean of fold rhos if available
def fold_cis_from_file(path, col="rho"):
    if not Path(path).exists():
        return None
    df = pd.read_csv(path)
    c = col if col in df.columns else df.columns[-1]
    return boot_mean_of_values(df[c].values, seed=3)

eep_folds = ROOT/"results/weekU_eep_loco_folds.csv"
z_folds = ROOT/"results/weekS_loco_zscore_folds.csv"
eep_fold_ci = fold_cis_from_file(eep_folds) or fold_cis_from_file("logs/weekM_loco_folds.tsv")
z_fold_ci = fold_cis_from_file(z_folds)

print("\n=== 1 LOCO OOF bootstrap Spearman ===")
print("EEP OOF:", fmt(eep_pt, eep_lo, eep_hi), "n=", eep_n)
print("z  OOF:", fmt(z_pt, z_lo, z_hi), "n=", z_n)
if eep_fold_ci:
    print("EEP mean-fold:", fmt(*eep_fold_ci[:3]), "nfolds=", eep_fold_ci[3])
if z_fold_ci:
    print("z  mean-fold:", fmt(*z_fold_ci[:3]), "nfolds=", z_fold_ci[3])

# ========== 2) Condition A ==========
ca = pd.read_csv(ROOT/"results/weekS_conditionA_zscore_vs_eep_preds.csv")
print("\n=== 2 Condition A ===")
ca_rows = {}
for name, ycol, pcol in [
    ("EEP_HT29", "lfc_HT29", "pred_eep"),
    ("EEP_A375", "lfc_A375", "pred_eep"),
    ("Z_HT29", "lfc_HT29", "pred_z"),
    ("Z_A375", "lfc_A375", "pred_z"),
]:
    # ρ(pred, -LFC)
    pt, lo, hi, n, _ = boot_spearman(-ca[ycol].values, ca[pcol].values, seed=10+hash(name)%100)
    ca_rows[name] = (pt, lo, hi, n)
    print(name, fmt(pt, lo, hi), "n=", n)

# ========== 3) Ablation Semana N (fold-level if available) ==========
print("\n=== 3 Ablation sequential ===")
abl_files = [
    ROOT/"logs/weekN_ablation_sequential.tsv",
    ROOT/"logs/weekN_ablation_seq_folds.tsv",
    ROOT/"results/ablation_results_RF.csv",
]
abl_summary = {}
for p in abl_files:
    if p.exists():
        df = pd.read_csv(p, sep="\t" if p.suffix==".tsv" else ",")
        print(p.name, "cols", list(df.columns), "shape", df.shape)
        print(df.head(8).to_string())

# Try sequential summary means
seq = ROOT/"logs/weekN_ablation_sequential.tsv"
if seq.exists():
    df = pd.read_csv(seq, sep="\t")
    # expect columns like step, mean_rho or similar
    print("sequential full:\n", df.to_string())

folds_abl = ROOT/"logs/weekN_ablation_seq_folds.tsv"
if folds_abl.exists():
    df = pd.read_csv(folds_abl, sep="\t")
    print("fold-level cols", list(df.columns))
    # group by step/ablation
    step_col = None
    for c in df.columns:
        if c.lower() in ("step","ablation","stage","features","config"):
            step_col = c
            break
    rho_col = None
    for c in df.columns:
        if "rho" in c.lower() or "spearman" in c.lower():
            rho_col = c
            break
    if step_col and rho_col:
        for step, g in df.groupby(step_col):
            pt, lo, hi, nf = boot_mean_of_values(g[rho_col].values, seed=20)
            abl_summary[str(step)] = (pt, lo, hi, nf)
            print("ABL", step, fmt(pt, lo, hi))

# ========== 4) Benchmarking RS3 + position ==========
print("\n=== 4 Benchmarking ===")
# RS3: need per-guide scores + EEP
rs3 = pd.read_csv(ROOT/"results/ruleset3_horlbeck_scored.csv")
eep_h = pd.read_csv(ROOT/"results/horlbeck_guides_eep.csv")
rs3["guide_sequence"] = rs3["guide_sequence"].astype(str).str.upper()
eep_h["guide_sequence"] = eep_h["guide_sequence"].astype(str).str.upper()
# find score col
sc = [c for c in rs3.columns if "score" in c.lower() or "rs3" in c.lower() or "ruleset" in c.lower()]
print("rs3 score cols", sc, "rs3", rs3.shape)
score_col = sc[0] if sc else rs3.columns[-1]
j = eep_h.merge(rs3[["guide_sequence", score_col]], on="guide_sequence", how="inner")
print("rs3 merge", j.shape)
rs3_pt, rs3_lo, rs3_hi, rs3_n, _ = boot_spearman(j["EEP_percentile"], j[score_col], seed=30)
print("RS3 global", fmt(rs3_pt, rs3_lo, rs3_hi), "n=", rs3_n)

# LOCO mean-fold CI from weekQ
for label, path in [("RS3_loco_folds", "results/weekQ_rs3_loco.csv"),
                    ("pos_loco_folds", "results/weekQ_position_baseline_loco.csv")]:
    p = ROOT/path
    if p.exists():
        df = pd.read_csv(p)
        print(label, df.head(), list(df.columns))
        rcol = "rho" if "rho" in df.columns else [c for c in df.columns if "rho" in c.lower()][0]
        pt, lo, hi, nf = boot_mean_of_values(df[rcol].values, seed=40)
        print(label, fmt(pt, lo, hi), "nfolds", nf)

# Position baseline nonlinear: recompute if we have distance + EEP
# use matrix nearest_tss_distance if present
if "nearest_tss_distance" in mat.columns or any("tss" in c.lower() for c in mat.columns):
    dist_col = "nearest_tss_distance" if "nearest_tss_distance" in mat.columns else \
        [c for c in mat.columns if "tss" in c.lower() and "dist" in c.lower()][0]
    # GBR LOCO would need retrain; use fold file if only means
    pass

# XGB full LOCO for benchmarking table = EEP OOF / fold CIs already computed

# ---- Write markdown ----
lines = []
lines.append("# Week U — Bootstrap 95% CIs (final, preprint-ready)\n")
lines.append(f"**N_boot:** {N_BOOT}  |  **Method:** resample pairs (y, pred) for OOF/Condition A; resample fold-ρ for LOCO-mean CIs when only fold table exists.\n")
lines.append("## 1. LOCO K562 (XGB 101 features)\n")
lines.append("| Target | Point (OOF Spearman) | 95% CI (bootstrap pairs) | Ready-to-cite |")
lines.append("|--------|----------------------|---------------------------|---------------|")
lines.append(f"| EEP | {eep_pt:.4f} | [{eep_lo:.4f}, {eep_hi:.4f}] | {fmt(eep_pt,eep_lo,eep_hi)} |")
lines.append(f"| z-score continuo | {z_pt:.4f} | [{z_lo:.4f}, {z_hi:.4f}] | {fmt(z_pt,z_lo,z_hi)} |")
if eep_fold_ci:
    lines.append(f"\nLOCO **mean of 23 fold ρ** EEP: {fmt(*eep_fold_ci[:3])} (n_folds={eep_fold_ci[3]})")
if z_fold_ci:
    lines.append(f"\nLOCO **mean of 23 fold ρ** z: {fmt(*z_fold_ci[:3])} (n_folds={z_fold_ci[3]})")
lines.append("\nNote: OOF pair bootstrap ≈ overall Spearman; fold-mean bootstrap reflects LOCO protocol variance across chromosomes.\n")
lines.append("## 2. Condition A (ρ pred vs −LFC, n≈100942)\n")
lines.append("| Model | Cohort | Ready-to-cite |")
lines.append("|-------|--------|---------------|")
for k, (pt, lo, hi, n) in ca_rows.items():
    lines.append(f"| {k} | — | {fmt(pt,lo,hi)} (n={n}) |")
lines.append("\n## 3. Ablation sequential (EEP, Semana N)\n")
if abl_summary:
    lines.append("| Step | Ready-to-cite |")
    lines.append("|------|---------------|")
    for k, (pt, lo, hi, nf) in abl_summary.items():
        lines.append(f"| {k} | {fmt(pt,lo,hi)} (n_folds={nf}) |")
else:
    lines.append("_Fold-level ablation table not auto-parsed; see printed head above — fill manually if needed._\n")
lines.append("\n## 4. Benchmarking\n")
lines.append("| Method | Ready-to-cite |")
lines.append("|--------|---------------|")
lines.append(f"| XGB EEP (OOF) | {fmt(eep_pt,eep_lo,eep_hi)} |")
lines.append(f"| Rule Set 3 (global Spearman vs EEP) | {fmt(rs3_pt,rs3_lo,rs3_hi)} (n={rs3_n}) |")
lines.append("\nOverlap check (XGB EEP CI vs RS3 CI): "
             f"{'NO overlap' if eep_lo > rs3_hi else 'OVERLAP'} "
             f"(XGB lo={eep_lo:.3f} vs RS3 hi={rs3_hi:.3f})\n")

out = ROOT/"logs/weekU_bootstrap_final.md"
out.write_text("\n".join(lines))
print("\nWROTE", out)
print("\n".join(lines))