# Semana L — Seis arquitecturas (re-cierre Sec. 6.4)

**Split:** train chr1–19+X, val chr20, test chr21/22/Y (test intacto).

## Cumplimiento de especificación
| Modelo | Spec Sec. 6.4 | Implementación final |
|--------|---------------|----------------------|
| Ridge/LASSO/RF | OK | sin cambio material |
| XGBoost | Optuna 100 trials | **100 trials** (antes 40; ρ 0.401→0.403) |
| MLP | dropout 0.3 + batch norm | **torch**, BN + Dropout(0.3) |
| CNN | 20-nt + flancos | **one-hot seq_ctx ±10**, pad 42; no tabular |

## S_interp (importante)
Asignación **categórica fija** por arquitectura (Ridge/LASSO=1.0, XGB=0.7, RF=0.5, MLP/CNN=0.3), **no medida empíricamente**. Con peso 30%, el protocolo **favorece estructuralmente** a modelos lineales si su ρ no es mucho peor. Ridge como primario es “correcto según el protocolo pre-especificado”, no “mejor que XGBoost en ranking univariado”.

## Tabla final (val)
| model | rho_val | cal_90 | S |
|-------|---------|--------|---|
| ridge | 0.335 | 0.914 | **0.727** PRIMARY |
| lasso | 0.327 | 0.907 | 0.724 |
| xgb   | 0.403 | 0.881 | 0.662 |
| rf    | 0.371 | 0.814 | 0.573 |
| mlp   | 0.365 | 0.853 | 0.520 |
| cnn   | 0.259 | 0.891 | 0.441 |

## Limitaciones
- LASSO bootstrap 20× (no 100): set “estable ≥70%” puede ser ruidoso / más amplio que el fit full.
- CNN stab no calculada → S_stab neutral 0.5.
- CNN ρ bajo refuerza que cromatina (F2–F5) aporta frente a secuencia sola.

## Artefactos
- models/M1_ridge_v1_101features.joblib (**primario**)
- models/M4_xgb_v1_101features.joblib
- models/M5_mlp_v1_101features.pt
- models/M6_cnn_v1_seqflanks.pt
