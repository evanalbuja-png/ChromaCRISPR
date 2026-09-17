# ChromaCRISPR

Multi-modal prediction of CRISPRi guide efficacy in **K562**, combining sequence and chromatin features, with a formal comparison of two target designs:

1. **EEP** — within-gene rank (ensemble efficacy percentile)  
2. **Continuous harmonized z-score** — within-study standardized efficacy (no within-gene rank)

Cross-cell-type transfer is evaluated on **Sanson et al. 2018** CRISPRi screens (**HT29**, **A375**).

> Preprint: *bioRxiv link to be added*  
> Archival: Zenodo (tag `v0.1.0`) · Intended track: PCI Mathematical & Computational Biology

---

## Main results

| Setting | Target | Spearman ρ (95% CI) |
|--------|--------|---------------------|
| K562 LOCO (23 folds) | EEP | 0.375 [0.360, 0.389] |
| K562 LOCO (23 folds) | Continuous z-score | 0.429 [0.412, 0.447] |
| Transfer HT29 | EEP | 0.008 [0.002, 0.014] |
| Transfer HT29 | Continuous z-score | 0.136 [0.129, 0.142] |
| Transfer A375 | EEP | 0.012 [0.006, 0.018] |
| Transfer A375 | Continuous z-score | 0.158 [0.152, 0.164] |

Full bootstrap tables: `logs/development_history/weekU_bootstrap_final.md` and CSVs under `results/`.

**Takeaway:** continuous z-score matches or beats EEP on K562 LOCO and **substantially** improves external transfer; within-gene ranking is a poor label for cross-cell generalization in this setting.

---

## Repository contents

| Path | What it is |
|------|------------|
| `data/processed/K562_training_matrix_v1.csv` | 19,446 guides × 101 features + labels (author-derived) |
| `data/processed/K562_training_labels_harmonized_v1.csv` | Harmonized EEP + z-score labels |
| `data/README.md` | ENCODE / screen accessions (**raw third-party data not shipped**) |
| `models/` | Trained models + `M4_feature_list_101.json` |
| `results/` | LOCO folds, ablation, Condition A, SHAP aggregates, OOF predictions |
| `figures/` | Figures 2–3 (ablation, SHAP by modality, target-design comparison) |
| `scripts/` | Pipeline scripts used for the analyses above |
| `logs/development_history/` | Week-by-week methodological log (Semanas A–U) |
| `environment.yml` | Conda environment |

**Not included (by design):** ENCODE bigWigs, Hi-C `.hic`/`.mcool`, full GENCODE GTF, third-party screen Excel/raw tables. See `data/README.md` to re-download if you need to rebuild features from scratch.

---

## Environment

```bash
conda env create -f environment.yml
conda activate chromacrispr-phase1
