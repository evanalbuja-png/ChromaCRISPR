# Week K.5 — Correlation audit (accepted features vs EEP)

**N accepted features scored:** 101
**Matrix:** K562_training_matrix_v1.csv (19446 rows)

## Top 10 |rho|
            feature       rho         pval     n
       poly_t_count -0.113654 6.386119e-57 19446
       guide_length -0.107249 7.564330e-51 19446
guide_3prime_pos4_C -0.102965 5.510782e-47 19446
     ATAC_pm100_p90  0.083726 1.346720e-31 19446
    ATAC_pm100_mean  0.082727 6.939330e-31 19446
     seed_pwm_score  0.075427 6.114615e-26 19446
     ATAC_pm100_max  0.073858 6.168779e-25 19446
guide_3prime_pos4_G  0.069874 1.753765e-22 19446
           ccre_PLS  0.064672 1.760250e-19 19446
    ATAC_pm200_mean  0.060853 2.010870e-17 19446

## Top 10 closest to zero
                    feature       rho     pval     n
           H3K27me3_missing -0.000009 0.998964 19446
            H3K4me3_missing -0.000009 0.998964 19446
          phylop100_missing -0.000018 0.997987 19446
log_dist_target_tss_missing -0.000037 0.995833 19446
       H3K27me3_pm1000_mean -0.000531 0.940952 19446
              chromhmm_Enh1 -0.000818 0.909170 19446
        H3K4me3_pm2000_mean  0.000886 0.901631 19446
       H3K27me3_pm2000_mean -0.001433 0.841608 19446
               chromhmm_Biv -0.001508 0.833480 19446
        H3K4me3_pm1000_mean -0.002668 0.709891 19446

## Expected-signal features (ATAC / H3K27ac / H3K4me3 / DNase / PLS / chromhmm_Tss)
            feature       rho         pval     n note
     ATAC_pm100_p90  0.083726 1.346720e-31 19446     
    ATAC_pm100_mean  0.082727 6.939330e-31 19446     
     ATAC_pm100_max  0.073858 6.168779e-25 19446     
           ccre_PLS  0.064672 1.760250e-19 19446     
    ATAC_pm200_mean  0.060853 2.010870e-17 19446     
     ATAC_pm200_p90  0.048487 1.332753e-11 19446     
     ATAC_pm200_max  0.044885 3.803944e-10 19446     
    ATAC_pm500_mean  0.040718 1.346589e-08 19446     
     ATAC_pm500_p90  0.034245 1.783449e-06 19446     
 H3K27ac_pm500_mean  0.033974 2.151482e-06 19446     
   ATAC_pm1000_mean  0.033415 3.151125e-06 19446     
     ATAC_pm500_max  0.030929 1.604689e-05 19446     
H3K27ac_pm1000_mean  0.028805 5.885800e-05 19446     
          ccre_pELS -0.027084 1.585435e-04 19446     
    ATAC_pm1000_p90  0.026564 2.116522e-04 19446     
    ATAC_pm1000_max  0.025007 4.875291e-04 19446     
   ATAC_pm2000_mean  0.021609 2.582488e-03 19446     
    ATAC_pm2000_p90  0.020954 3.476226e-03 19446     
       chromhmm_Tss  0.020526 4.204392e-03 19446     
H3K27ac_pm2000_mean  0.020463 4.321274e-03 19446     
    ATAC_pm2000_max  0.017149 1.678495e-02 19446     
   dnase_200bp_mean  0.013973 5.135071e-02 19446     
  chromhmm_TssFlnkU  0.011993 9.446050e-02 19446     
    ATAC_pm5000_p90  0.011306 1.149012e-01 19446     
  chromhmm_TssFlnkD -0.010568 1.405957e-01 19446     
    ATAC_pm5000_max  0.009138 2.025741e-01 19446     
 H3K4me3_pm500_mean -0.005721 4.250204e-01 19446     
   chromhmm_TssFlnk -0.004420 5.376416e-01 19446     
   ATAC_pm5000_mean  0.003312 6.442310e-01 19446     
H3K4me3_pm1000_mean -0.002668 7.098912e-01 19446     
H3K4me3_pm2000_mean  0.000886 9.016315e-01 19446     
    H3K4me3_missing -0.000009 9.989641e-01 19446     

## Multi-gen stability (F4/F5)
(see console output: ALL vs NO_MULTI vs ONLY_MULTI)

## Missingness flags
(see console)

## Features meriting review
- No features with |rho| > 0.3
### Expected-signal with |rho| < 0.01
            feature       rho     pval     n note
    ATAC_pm5000_max  0.009138 0.202574 19446     
 H3K4me3_pm500_mean -0.005721 0.425020 19446     
   chromhmm_TssFlnk -0.004420 0.537642 19446     
   ATAC_pm5000_mean  0.003312 0.644231 19446     
H3K4me3_pm1000_mean -0.002668 0.709891 19446     
H3K4me3_pm2000_mean  0.000886 0.901631 19446     
    H3K4me3_missing -0.000009 0.998964 19446     

## Conclusion
(fill after reviewing console: anomalies or clean audit)
## Multi-gen stability (números)
| Feature | ALL | NO_MULTI | ONLY_MULTI |
|---|---|---|---|
| hic_contact_log | -0.0131 | -0.0146 | 0.0615 (n=130, NS) |
| log_dist_target_tss | -0.0196 | -0.0205 | 0.1015 (n=144, NS) |

## Missingness flags
Todas |rho| < 1e-4 vs EEP.

## Features meriting review
- **Ninguna con |rho|>0.3.**
- **H3K4me3_* ~0:** no indica pipeline roto (QC Semana H OK); plausible por diseño EEP within-gene + sesgo promotor. No bloquea entrenamiento; el modelo multivariado puede usar interacciones.
- **pam_is_NGG / columnas casi constantes:** ConstantInputWarning esperado; poca utilidad univariada, retenidas solo si el manifest las marca accepted.

## Conclusión
**Auditoría completa, sin anomalías bloqueantes.** Señal univariada débil pero direccionalmente coherente (ATAC, H3K27ac, PLS, poly_t). Sin leakage moderado-alto. F4/F5 estables fuera de multi-gen. **Semana L habilitada.**
