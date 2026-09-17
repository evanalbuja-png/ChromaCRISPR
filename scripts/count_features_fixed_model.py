#!/usr/bin/env python
# scripts/count_features_fixed_model.py
import json
from pathlib import Path
import pandas as pd

feat_path = Path("models/xgb_phase1_feature_list_FIXED.json")
mat_path = Path("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_FIXED.csv")

with open(feat_path) as f:
    feats = json.load(f)

print(f"N features en lista del modelo: {len(feats)}")
print("Lista completa:")
for i, c in enumerate(feats, 1):
    group = (
        "F1" if c in ["gc_content", "mfe_rnafold", "guide_length", "g_run_max", "poly_t_count"] else
        "F2" if c.startswith("ATAC") else
        "F3" if c.startswith("H3K") else
        "F5" if c in ["nearest_tss_distance", "within_promoter_2kb"] else
        "?"
    )
    print(f"  {i:2d}. [{group}] {c}")

from collections import Counter
groups = []
for c in feats:
    if c in ["gc_content", "mfe_rnafold", "guide_length", "g_run_max", "poly_t_count"]:
        groups.append("F1")
    elif c.startswith("ATAC"):
        groups.append("F2")
    elif c.startswith("H3K"):
        groups.append("F3")
    elif c in ["nearest_tss_distance", "within_promoter_2kb"]:
        groups.append("F5")
    else:
        groups.append("?")
print("\nPor grupo:", dict(Counter(groups)))
print("F4 (Hi-C): no incluido")
print("Proposal total previsto: 104 (F1=24, F2=18, F3=42, F4=8, F5=12)")
print(f"Usado en Phase 1 FIXED: {len(feats)}")