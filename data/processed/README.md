# Processed matrices (author-derived)

## K562_training_matrix_v1.csv
- N = 19,446 guides (Horlbeck + Gasperini high-confidence)
- 101 predictive features (F1–F5) + identifiers
- Targets: `EEP_percentile`, `z_score_within_study` (continuous harmonized z-score)
- Feature list: `models/M4_feature_list_101.json`

## K562_training_labels_harmonized_v1.csv
- Harmonization: within-study z-score → optional within-gene percentile (EEP)
- Columns include: guide_sequence, chromosome, coordinate, gene_target, dataset_origin, raw_score, z_score_within_study, EEP_percentile

Rebuild from raw sources is documented in `logs/development_history/` and `scripts/`.
