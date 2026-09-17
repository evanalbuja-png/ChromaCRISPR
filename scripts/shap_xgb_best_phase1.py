# Librerías
import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt
from sklearn.model_selection import GroupKFold
from scipy.stats import spearmanr
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from xgboost import XGBRegressor
from xgboost import DMatrix

# Importar matriz
df = pd.read_csv(
    'data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv'
)

# Separar por grupos
y = df["efficacy_score"]
groups = df["nearest_tss_gene"]

# Crear objeto GroupKFold
gkf = GroupKFold(n_splits=5)
folds = gkf.split(df, y, groups)
folds = list(folds)


# Crear listas de las columnas

F1 = [
    "gc_content",
    "mfe_rnafold",
    "guide_length",
    "g_run_max",
    "poly_t_count"
]

F2 = [
    "ATAC_mean",
    "ATAC_max",
    "ATAC_p90",
    "ATAC_sum"
]

F3 = [
    "H3K27ac_mean",
    "H3K27ac_max",
    "H3K27ac_p90",
    "H3K27ac_sum",
    "H3K4me3_mean",
    "H3K4me3_max",
    "H3K4me3_p90",
    "H3K4me3_sum",
    "H3K27me3_mean",
    "H3K27me3_max",
    "H3K27me3_p90",
    "H3K27me3_sum"
]

F5 = [
    "nearest_tss_distance",
    "within_promoter_2kb"
]


# Creación de combinaciones
F1_F2_F3_F5 = F1 + F2 + F3 + F5

# ============================================================
# XGBOOST
# ============================================================

resultados_xgb = []

X_global = []

shap_global = []

for fold, (train_idx, test_idx) in enumerate(folds, start=1):

    print("\nFold:", fold)

    X_train = df.loc[train_idx, F1_F2_F3_F5]
    X_test = df.loc[test_idx, F1_F2_F3_F5]

    y_train = y.loc[train_idx]
    y_test = y.loc[test_idx]

    # Rellenar los datos
    medianas = X_train.median()
    X_train = X_train.fillna(medianas)
    X_test = X_test.fillna(medianas)

    # Modelo
    xgb = XGBRegressor(
        n_estimators=500,
        max_depth=3,
        learning_rate=0.05,
        min_child_weight=5,
        random_state=42,
        n_jobs=-1
    )

    # Entrenamiento
    xgb.fit(X_train, y_train)

    y_test_pred_xgb = xgb.predict(X_test)

    booster = xgb.get_booster()
    X_test_dmatrix = DMatrix(X_test)

    # TreeSHAP 
    
    shap_values = booster.predict(
        X_test_dmatrix,
        pred_contribs=True
    )

    X_global.append(X_test.copy())
    shap_global.append(shap_values)

    # Predicción

    pred_from_shap = shap_values.sum(axis=1)

    print(
        "SHAP reconstruction max error:",
        np.max(np.abs(pred_from_shap - y_test_pred_xgb))
    )

    # Métricas
    test_spearman_xgb = spearmanr(
        y_test,
        y_test_pred_xgb
    ).statistic

    test_r2_xgb = r2_score(
        y_test,
        y_test_pred_xgb
    )

    test_rmse_xgb = np.sqrt(
        mean_squared_error(
            y_test,
            y_test_pred_xgb
        )
    )

    print("XGB Test Spearman:", test_spearman_xgb)
    print("XGB Test R²:", test_r2_xgb)
    print("XGB Test RMSE:", test_rmse_xgb)

    # Guardar resultados
    resultados_xgb.append({
        "fold": fold,
        "spearman": test_spearman_xgb,
        "r2": test_r2_xgb,
        "rmse": test_rmse_xgb
    })

    pd.DataFrame(resultados_xgb).to_csv(
    "results/xgb_shap_results.csv",
    index=False
    )

# ============================================================
# RANKING GLOBAL SHAP
# ============================================================

shap_global = np.vstack(shap_global)

# Última columna = bias
shap_features = shap_global[:, :-1]

feature_names = F1_F2_F3_F5

shap_importance = pd.DataFrame({
    "feature": feature_names,
    "mean_abs_shap": np.mean(
        np.abs(shap_features),
        axis=0
    )
})

shap_importance = shap_importance.sort_values(
    "mean_abs_shap",
    ascending=False
)

print("\nTop 15 SHAP features:")
print(shap_importance.head(15))

shap_importance.to_csv(
    "results/xgb_shap_feature_importance.csv",
    index=False
)

# ============================================================
# DATOS GLOBALES PARA DEPENDENCE PLOTS
# ============================================================

X_global = pd.concat(
    X_global,
    axis=0,
    ignore_index=True
)

shap_global = np.vstack(shap_global)

shap_features = shap_global[:, :-1]

print("X_global shape:", X_global.shape)
print("SHAP features shape:", shap_features.shape)

# ============================================================
# SHAP DEPENDENCE PLOTS
# ============================================================

features_dependence = [
    "nearest_tss_distance",
    "ATAC_mean",
    "ATAC_max",
    "H3K27ac_mean",
    "mfe_rnafold"
]

for feature in features_dependence:

    feature_idx = feature_names.index(feature)

    shap_feature = shap_features[:, feature_idx]
    x_feature = X_global[feature].values

    plt.figure(figsize=(7, 5))

    print(
        feature,
        "X:",
        len(x_feature),
        "SHAP:",
        len(shap_feature)
    )

    plt.scatter(
        x_feature,
        shap_feature,
        alpha=0.25,
        s=10
    )

    plt.axhline(
        0,
        linewidth=1
    )

    plt.xlabel(feature)
    plt.ylabel("SHAP value")
    plt.title(f"SHAP dependence: {feature}")

    plt.tight_layout()

    plt.savefig(
        f"results/shap_dependence_{feature}.png",
        dpi=300
    )

    plt.close()