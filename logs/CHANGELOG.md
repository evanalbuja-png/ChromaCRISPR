# ChromaCRISPR Phase 1 — CHANGELOG

Formato: YYYY-MM-DD | Acción | Archivo/versión | Motivo

2026-09-09 | Inventario inicial | data/processed/, models/, interim/, raw/ | Punto de partida real documentado. Matriz activa: sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_FIXED.csv. Solo XGBoost entrenado. DS1/DS3/DS4 referencias existen pero integración rank-based y QC formal Sec. 7.1 pendiente de auditoría.
2026-09-09 | Semana A cerrada | horlbeck2016_reference_FINAL.csv + *_FIXED.csv | QC Sec. 7.1 parcial: hg38 y qc_pass presentes, pero re-mapeo Bowtie2 zero-mismatch y estadísticas de unique mapping no documentados. efficacy_score continuo (no EEP rank-based). Limitación heredada aceptada.
2026-09-09 | Semana E cerrada | sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_FIXED.csv | 23 features vs 104 especificadas (F1=5/24, F2=4/18, F3=12/42, F4=0/8, F5=2/12). Tabla de auditoría completa en este mensaje.
2026-09-09 | Semana B cerrada (definitiva) | data/interim/DS3_gasperini/DS3_gasperini_guides_efficacy_highconf_polarity_fixed.csv | 1.328 guías high-confidence. Score = -log2(proportion_remaining + 0.001). Polaridad alineada con Horlbeck. Desviación ≥50 células y limitación de mapeo documentadas.
2026-09-10 | Semana C cerrada (Opción 2) | DS4 Replogle | No hay tabla de efficacy scores oficial y manejable en el proyecto ni identificada en búsqueda acotada. DS4 excluido de Semana D. Coordenadas reference conservadas.
2026-09-10 | Semana D verificación nomenclatura | K562_training_labels_harmonized_v1.csv | Confirmado: 260 genes ENSG en Horlbeck no fragmentan genes existentes (intersección = 0). Sin v2. Bloque 1 cerrado definitivamente.
2026-09-10 | Semana F cerrada | data/interim/F1_sequence_features_v1.csv | F1: 6/24 oficiales (gc, mfe, g_run, poly_t, tm_santalucia, seed_pwm). PAM context real pendiente. Proxy 3' documentado. Spearman PWM-EEP=0.07 (sin leakage).
2026-09-10 | Semana G cerrada | data/interim/F2_accessibility_features_v1.csv | F2: 18/18 oficiales (6 ventanas × mean/max/p90). Spearman ±500 vs FIXED = 1.0. DNase + peak features bloqueados (archivos ausentes).
2026-09-11 | Semana G cerrada (completa) | data/interim/F2_accessibility_features_v1.csv | F2: 18 oficiales + 3 extras (DNase ENCFF452XDU, peaks ENCFF993BAP). Spearman ±500=1.0. 0 NaNs.
2026-09-11 | Semana H cerrada | data/interim/F3_histone_chromhmm_features_v1.csv | F3: 31 oficiales (15 mean marks + CTCF + ChromHMM15). n=42 del proposal documentado como discrepancia. Accessions H3K4me1/CTCF corregidos vs Sec. 7.2. Spearman vs FIXED ≈1.0.
2026-09-11 | Semana H imputación corregida | F3 | NaN→mediana + H3K4me3_missing/H3K27me3_missing (Sec. 7.3). Accessions proposal incorrectos documentados vía metadata ENCODE.
2026-09-11 | Semana H CERRADA | data/interim/F3_histone_chromhmm_features_v1.csv | 19446 filas, 31 oficiales, mediana+flags Sec.7.3, accessions proposal incorrectos documentados vía ENCODE metadata.
2026-09-11 | Semana I cerrada (parcial) | data/interim/F4_hic_features_v1.csv | F4: hic_contact_log OK (ρ_dist=-0.90); compartment E1 calculado pero sanity ATAC falló (no forzar). TAD/loops pendiente.
2026-09-11 | Semana I corrección final | data/interim/F4_hic_features_v1.csv | Fix: compartment_A propaga NaN de E1 correctamente (786). Semana I cerrada definitivamente.
2026-09-12 | Semana J cerrada (parcial) | data/interim/F5_genomic_context_features_v1.csv | F5: 7 oficiales (dists, cCRE×4, PhyloP). GERP/JASPAR pendientes. PAM extraído pero no validado (F1 sigue pendiente).
2026-09-12 | PAM F1 validado | F5/F1 | Match guía/RC ±80bp → NGG 99.99% (n=18433). 3/3 controles OK. NaN F5 solo strand Gasperini.
2026-09-12 | Fix NaN condicional | F5 PAM + nota general | pam_is_NGG/one-hot: 1013 NaN propagados (no 0). Patrón recurrente: nunca codificar "no calculado" como False/0.
2026-09-12 | Semana K cerrada | data/processed/K562_training_matrix_v1.csv + K562_feature_manifest_v1.tsv | 19446×120, sin fan-out. Manifest con accepted_for_modeling.
2026-09-12 | Manifest v1 ajustado | distance_to_nearest_atac_peak → F2_extra accepted; hic_contact_raw nota corregida.
2026-09-12 | Bloque 2 CERRADO (verificado) | Matriz 19446×120; PAM NaN OK; multi-gen F4/F5 distintos; anti-leakage OK.
2026-09-12 | Week K.5 correlation audit | logs/weekK5_correlation_audit.md | 101 features vs EEP; max|rho|~0.11; no blockers; Semana L OK.
2026-09-12 | Semana L cerrada | models/M*_v1_101features.* | Primario=Ridge (S=0.727). XGB mejor rho=0.40. Split cromosoma Sec.7.3. Test intacto.
2026-09-12 | Semana L re-cerrada | MLP torch BN+dropout; CNN seq+flanks; XGB 100 trials. Primario=Ridge (S=0.727) por protocolo; S_interp categórico documentado.
2026-09-12 | Semana M cerrada | LOCO 23-fold Ridge 0.303±0.037; XGB 0.377±0.039. Métrica oficial = LOCO mean±SD. Sin outliers Ridge.
2026-09-13 | Semana N cerrada | Ablation XGB LOCO: F2 +0.032 p≈0; F5 +0.018 p≈0; F3/F4 null marginal. RQ1 sí cromatina; RQ2 F1>F2≈F5>>F3/F4.
2026-09-13 | Semana O cerrada | Condition A 101 feat: all-set |ρ|~0.01 (full < seq); overlap 526 |ρ|~0.15 full>seq. Brecha = transferencia cell-type, no dimensionalidad.
2026-09-13 | Semana O control target | EEP→Sanson ALL |ρ|~0.01 vs efficacy_score continuo |ρ|~0.15–0.17. Causa dominante: target within-gene, no domain shift. 
2026-09-13 | Semana P cerrada | TreeSHAP via pred_contribs; F1>F2>F5>>F3>F4; coherente ablation N; I3 N/A.
2026-09-13 | Semana P extensión | LOCO-OOF SHAP n=19446 (ρ vs test490=0.98); I1 top100: 57% chromatin+ / 61% seq-dominated. Cierre robusto.
2026-09-13 | Semana Q cerrada | Position LOCO ρ=0.021; RS3 LOCO ρ=0.259 (proxy RS2); DeepCRISPR bloqueante TF1; XGB 0.377 lidera.
2026-09-13 | Semana Q position-baseline diagnostic | target vs nearest_expressed equivalentes; linear ρ≈0.02 real; GBR univariado ρ≈0.12; no contradice F5 ablation/SHAP.
2026-09-13 | Semana R cerrada (definitivo) | RS2/DeepCRISPR bloqueantes de infraestructura confirmados tras segundo intento en entorno aislado. Benchmarking Sec.12 cierra con RS3 + position-baseline como comparadores disponibles.

## 2026-09-14 — Week S target redesign
- LOCO z-score continuo: ρ=0.429±0.044 (vs EEP 0.377±0.039)
- Condition A z vs EEP: HT29 0.136 vs 0.008; A375 0.158 vs 0.012; perm p≈0.0001
- Ablation F1 0.371 → FULL 0.429 (Δ+0.058) under continuous target
- Log: logs/weekS_target_redesign.md

- Sanity: LOCO z Horlbeck-only 0.448 > H+G 0.429 (no inflation by study ID); Cond.A z 0.136/0.158 slightly below Week-O Horlbeck-only control ~0.15/0.17

## 2026-09-14 — Week U bootstrap CIs
- LOCO EEP 0.374 [0.362, 0.386]; z 0.430 [0.418, 0.442] (no overlap)
- Cond.A z HT29/A375 0.136/0.158 CIs vs EEP 0.008/0.012 (no overlap)
- Ablation fold-mean CIs; XGB vs RS3 no CI overlap
- Log: logs/weekU_bootstrap_final.md
