#!/usr/bin/env python
# scripts/check_oof_split_type.py
"""Inspeccionar scripts de entrenamiento usados para el OOF 0.28/0.30."""

from pathlib import Path
import re

scripts = [
    "scripts/fix4_retrain_xgb_on_fixed.py",
    "scripts/train_and_save_xgb_phase1_conditionA.py",
    "scripts/baseline_sequence_only_conditionA.py",
    "scripts/shap_xgb_best_phase1.py",
    "scripts/run_baseline_models_horlbeck_week8.py",
    "scripts/baseline_models_phase1_week8.py",
]

print("=== Tipo de split en scripts de entrenamiento ===\n")
for sp in scripts:
    p = Path(sp)
    if not p.exists():
        print(f"[missing] {sp}")
        continue
    text = p.read_text(encoding="utf-8", errors="ignore")
    has_groupkfold = "GroupKFold" in text
    has_chromosome = "chromosome" in text.lower() and ("GroupKFold" in text or "groups" in text)
    has_random_split = "train_test_split" in text or "KFold" in text and "GroupKFold" not in text
    has_loco = "leave" in text.lower() and "chrom" in text.lower()

    print(f"FILE: {sp}")
    print(f"  GroupKFold:     {has_groupkfold}")
    print(f"  groups=chromosome likely: {has_chromosome}")
    print(f"  train_test_split / plain KFold: {has_random_split}")
    print(f"  LOCO explicit:  {has_loco}")

    # Extract relevant lines
    for i, line in enumerate(text.splitlines(), 1):
        if any(k in line for k in ["GroupKFold", "train_test_split", "KFold", "groups", "chromosome", "LOCO", "split"]):
            if line.strip() and not line.strip().startswith("#"):
                print(f"  L{i}: {line.strip()[:120]}")
    print()