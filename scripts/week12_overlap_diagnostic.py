#!/usr/bin/env python
# scripts/week12_overlap_diagnostic.py
"""
Por qué en las 526 guías del overlap el modelo full (ρ≈-0.13)
rinde peor que efficacy_score crudo (ρ≈-0.35) vs LFC Sanson.
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

HORLBECK = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_FIXED.csv")
SANSON = Path("data/interim/conditionA/sanson_feature_matrix_conditionA.csv")
PRED = Path("results/conditionA_sanson_predictions_FIXED.csv")
OUT_MD = Path("logs/week12_overlap_diagnostic.md")
OUT_CSV = Path("results/week12_overlap_diagnostic_table.csv")

FEATURES = [
    "gc_content", "mfe_rnafold", "guide_length", "g_run_max", "poly_t_count",
    "ATAC_mean", "ATAC_max", "ATAC_p90", "ATAC_sum",
    "H3K27ac_mean", "H3K27ac_max", "H3K27ac_p90", "H3K27ac_sum",
    "H3K4me3_mean", "H3K4me3_max", "H3K4me3_p90", "H3K4me3_sum",
    "H3K27me3_mean", "H3K27me3_max", "H3K27me3_p90", "H3K27me3_sum",
    "nearest_tss_distance", "within_promoter_2kb",
]
CHROM_FEATURES = [f for f in FEATURES if f.startswith(("ATAC", "H3K"))]

def main():
    h = pd.read_csv(HORLBECK)
    s = pd.read_csv(SANSON)
    p = pd.read_csv(PRED)

    # Unir overlap
    ov = (
        h[["guide_sequence", "efficacy_score", "chromosome"] + FEATURES]
        .merge(s[["guide_sequence", "lfc_HT29", "lfc_A375"] + FEATURES],
               on="guide_sequence", suffixes=("_h", "_s"))
        .merge(p[["guide_sequence", "pred_efficacy"]], on="guide_sequence", how="left")
    )
    print(f"Overlap n={len(ov)}")

    # Missingness en features Sanson (las usadas para predecir)
    feat_s = [f + "_s" if f + "_s" in ov.columns else f for f in FEATURES]
    # Tras merge con suffixes, columnas de Sanson features
    sanson_feat_cols = []
    for f in FEATURES:
        if f + "_s" in ov.columns:
            sanson_feat_cols.append(f + "_s")
        elif f in ov.columns and f + "_h" not in ov.columns:
            sanson_feat_cols.append(f)

    # Reconstruir de forma más clara desde SANSON + PRED + HORLBECK
    ov2 = (
        p.merge(h[["guide_sequence", "efficacy_score"]], on="guide_sequence", how="inner")
         .merge(s[["guide_sequence"] + FEATURES], on="guide_sequence", how="left",
                suffixes=("", "_sanson_raw"))
    )
    # p ya tiene features + pred; h trae efficacy_score
    # Verificar missing en p (matriz usada para predecir)
    missing_flags = p[FEATURES].isna()
    p["_n_missing"] = missing_flags.sum(axis=1)
    p["_has_any_missing"] = p["_n_missing"] > 0
    p["_n_missing_chrom"] = missing_flags[CHROM_FEATURES].sum(axis=1)

    ov3 = p.merge(h[["guide_sequence", "efficacy_score"]], on="guide_sequence", how="inner")
    print(f"Overlap vía predicciones: {len(ov3)}")

    lines = []
    lines.append("# Week 12 — Diagnóstico overlap Horlbeck↔Sanson (526 guías)\n")
    lines.append("## Objetivo\n")
    lines.append(
        "Explicar por qué el modelo full (ρ≈−0.13) rinde peor que "
        "efficacy_score crudo (ρ≈−0.35) vs LFC en el overlap de secuencia idéntica.\n"
    )

    # --- 1. Correlaciones centrales ---
    lines.append("## 1. Correlaciones en el overlap\n")
    rows = []
    for cell in ["lfc_HT29", "lfc_A375"]:
        rho_pred, _ = spearmanr(ov3["pred_efficacy"], ov3[cell], nan_policy="omit")
        rho_raw, _ = spearmanr(ov3["efficacy_score"], ov3[cell], nan_policy="omit")
        rho_model_vs_target, _ = spearmanr(ov3["pred_efficacy"], ov3["efficacy_score"], nan_policy="omit")
        rows.append((cell, rho_pred, rho_raw))
        lines.append(f"- **{cell}**: modelo vs LFC ρ={rho_pred:.4f} | "
                     f"efficacy_score vs LFC ρ={rho_raw:.4f}")
    rho_mt, _ = spearmanr(ov3["pred_efficacy"], ov3["efficacy_score"], nan_policy="omit")
    lines.append(f"- **Modelo vs efficacy_score (target Horlbeck)**: ρ={rho_mt:.4f}")
    lines.append("")
    lines.append(
        "Si ρ(modelo, efficacy_score) es alto pero ρ(modelo, LFC) es bajo → "
        "el modelo predice bien el target de training pero ese target no se "
        "transfiere igual de bien que el score crudo (poco probable si son el mismo número).\n"
    )
    lines.append(
        "Si ρ(modelo, efficacy_score) ya es bajo en el overlap → degradación "
        "del modelo *dentro* del subset (generalización / features).\n"
    )

    # --- 2. Missingness ---
    lines.append("## 2. Missingness / imputación\n")
    # Full Sanson
    full_miss = p["_n_missing"]
    full_miss_chrom = p["_n_missing_chrom"]
    ov_miss = ov3["_n_missing"]
    ov_miss_chrom = ov3["_n_missing_chrom"]

    lines.append("| Set | N | % con ≥1 missing | mean n_missing | % con ≥1 missing cromatina |")
    lines.append("|-----|---|------------------|----------------|----------------------------|")
    lines.append(
        f"| Sanson full | {len(p):,} | {100*p['_has_any_missing'].mean():.2f}% | "
        f"{full_miss.mean():.3f} | {100*(full_miss_chrom>0).mean():.2f}% |"
    )
    lines.append(
        f"| Overlap 526 | {len(ov3):,} | {100*(ov_miss>0).mean():.2f}% | "
        f"{ov_miss.mean():.3f} | {100*(ov_miss_chrom>0).mean():.2f}% |"
    )
    lines.append("")

    # Por feature
    lines.append("### Missing rate por feature (overlap vs full)\n")
    lines.append("| Feature | % missing full | % missing overlap |")
    lines.append("|---------|----------------|-------------------|")
    for f in FEATURES:
        mf = 100 * p[f].isna().mean()
        mo = 100 * ov3[f].isna().mean()
        if mf > 0.01 or mo > 0.01:
            lines.append(f"| {f} | {mf:.2f}% | {mo:.2f}% |")
    lines.append("")

    # --- 3. ¿El modelo predice bien el target en el overlap? ---
    lines.append("## 3. Modelo vs efficacy_score en el overlap\n")
    lines.append(f"- Spearman(pred, efficacy_score) = **{rho_mt:.4f}**")
    # En training OOF era ~0.30; si aquí es mucho menor → subset difícil
    resid = ov3["pred_efficacy"] - ov3["efficacy_score"]
    lines.append(f"- Residual (pred − score): mean={resid.mean():.4f}, "
                 f"std={resid.std():.4f}, MAE={resid.abs().mean():.4f}")
    lines.append("")

    # --- 4. Top 10 discrepancias ---
    lines.append("## 4. Top 10 |predicción − efficacy_score|\n")
    ov3 = ov3.copy()
    ov3["abs_err"] = (ov3["pred_efficacy"] - ov3["efficacy_score"]).abs()
    top = ov3.nlargest(10, "abs_err")[
        ["guide_sequence", "efficacy_score", "pred_efficacy", "abs_err",
         "lfc_HT29", "lfc_A375", "ATAC_mean", "H3K27ac_mean", "H3K4me3_mean",
         "nearest_tss_distance", "_n_missing"]
    ]
    lines.append("```")
    lines.append(top.to_string(index=False))
    lines.append("```\n")

    # --- 5. Clasificación ---
    lines.append("## 5. Clasificación de causa (con evidencia)\n")

    # Evidencia automática
    miss_diff = abs((ov_miss > 0).mean() - p["_has_any_missing"].mean())
    model_vs_target_weak = abs(rho_mt) < 0.20
    model_vs_lfc_much_worse = True  # already known ~0.13 vs 0.35

    if miss_diff > 0.05:
        cause = "imputación diferencial"
        evidence = (
            f"Tasa de missing en overlap ({100*(ov_miss>0).mean():.1f}%) "
            f"difiere del full ({100*p['_has_any_missing'].mean():.1f}%)."
        )
    elif model_vs_target_weak:
        cause = "generalización del modelo / features en el subset"
        evidence = (
            f"Spearman(pred, efficacy_score) en overlap = {rho_mt:.3f}, "
            f"bajo respecto al OOF global (~0.30). El modelo ya falla en "
            f"reproducir el target Horlbeck dentro de este subset."
        )
    else:
        cause = "generalización del ranking / optimización global"
        evidence = (
            f"Spearman(pred, efficacy_score) en overlap = {rho_mt:.3f} "
            f"(razonable vs OOF ~0.30), pero pred vs LFC ({rows[0][1]:.3f}) "
            f"<< score crudo vs LFC ({rows[0][2]:.3f}). El score crudo "
            f"conserva señal fenotípica punto-a-punto que el modelo, "
            f"optimizado para el ranking global de ~18k guías, no preserva "
            f"igual de bien en este subset de transferencia."
        )

    lines.append(f"**Causa más probable: {cause}**\n")
    lines.append(f"Evidencia: {evidence}\n")
    lines.append(
        "Nota: el score crudo *es* el target de esas 526 guías; cualquier "
        "compresión/regresión a la media del modelo reduce la correlación "
        "con un segundo fenotipo (LFC Sanson) si la relación efficacy→LFC "
        "es aproximadamente lineal y monotónica.\n"
    )

    # Guardar
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    ov3.to_csv(OUT_CSV, index=False)

    print("\n".join(lines))
    print(f"\n=== Guardado: {OUT_MD} ===")
    print(f"=== Tabla: {OUT_CSV} ===")

if __name__ == "__main__":
    main()