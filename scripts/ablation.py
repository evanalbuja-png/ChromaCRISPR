# Librerías
import pandas as pd
import numpy as np
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import spearmanr
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from xgboost import XGBRegressor


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
F1_F2 = F1 + F2
F1_F2_F3 = F1 + F2 + F3
F1_F2_F3_F5 = F1 + F2 + F3 + F5

ablations = {
    "F1": F1,
    "F1+F2": F1_F2,
    "F1+F2+F3": F1_F2_F3,
    "F1+F2+F3+F5": F1_F2_F3_F5
}


# ============================================================
# RANDOM FOREST — ABLATION
# ============================================================

rf = RandomForestRegressor(
    n_estimators=500,
    max_depth=None,
    min_samples_leaf=5,
    max_features="sqrt",
    bootstrap=True,
    oob_score=True,
    criterion="squared_error",
    random_state=42,
    n_jobs=-1
)

resultados = []

for nombre, features in ablations.items():

    print("\nAblation:", nombre)

    resultados_spearman = []
    resultados_r2 = []
    resultados_rmse = []

    for fold, (train_idx, test_idx) in enumerate(folds, start=1):

        print("\nFold:", fold)

        X_train = df.loc[train_idx, features]
        X_test = df.loc[test_idx, features]

        y_train = y.loc[train_idx]
        y_test = y.loc[test_idx]

        # Rellenar los datos
        medianas = X_train.median()
        X_train = X_train.fillna(medianas)
        X_test = X_test.fillna(medianas)

        # Entrenamiento
        rf.fit(X_train, y_train)

        y_train_pred = rf.predict(X_train)
        y_test_pred = rf.predict(X_test)

        # Índice de correlación de Spearman
        train_spearman = spearmanr(
            y_train,
            y_train_pred
        ).statistic

        test_spearman = spearmanr(
            y_test,
            y_test_pred
        ).statistic

        resultados_spearman.append(test_spearman)

        # Cálculo del R2
        train_r2 = r2_score(
            y_train,
            y_train_pred
        )

        test_r2 = r2_score(
            y_test,
            y_test_pred
        )

        resultados_r2.append(test_r2)

        # Cálculo de RMSE
        train_rmse = np.sqrt(
            mean_squared_error(
                y_train,
                y_train_pred
            )
        )

        test_rmse = np.sqrt(
            mean_squared_error(
                y_test,
                y_test_pred
            )
        )

        resultados_rmse.append(test_rmse)

        print("Train Spearman:", train_spearman)
        print("Test Spearman:", test_spearman)
        print("Train R²:", train_r2)
        print("Test R²:", test_r2)
        print("Train RMSE:", train_rmse)
        print("Test RMSE:", test_rmse)
        print("Train size:", len(y_train))
        print("Test size:", len(y_test))

    # Imprimir promedios
    resultados.append({
        "ablation": nombre,
        "spearman": np.mean(resultados_spearman),
        "r2": np.mean(resultados_r2),
        "rmse": np.mean(resultados_rmse)
    })

    resultados_df = pd.DataFrame(resultados)
    print(resultados_df)


resultados_df.to_csv(
    "results/ablation_results_RF.csv",
    index=False
)


# ============================================================
# XGBOOST — TUNING
# ============================================================

param_grid = [
    {"max_depth": 3, "learning_rate": 0.05, "min_child_weight": 5},
    {"max_depth": 3, "learning_rate": 0.1, "min_child_weight": 5},
    {"max_depth": 5, "learning_rate": 0.05, "min_child_weight": 5},
    {"max_depth": 5, "learning_rate": 0.1, "min_child_weight": 5},
]

resultados_xgb = []

for nombre, features in ablations.items():
    
    print("\nAblation:", nombre)
    for params in param_grid:

        print("\nParámetros:", params)

        for fold, (train_idx, test_idx) in enumerate(folds, start=1):

            print("\nFold:", fold)

            X_train = df.loc[train_idx, features]
            X_test = df.loc[test_idx, features]

            y_train = y.loc[train_idx]
            y_test = y.loc[test_idx]

            # Rellenar los datos
            medianas = X_train.median()
            X_train = X_train.fillna(medianas)
            X_test = X_test.fillna(medianas)

            # Modelo
            xgb = XGBRegressor(
                n_estimators=500,
                max_depth=params["max_depth"],
                learning_rate=params["learning_rate"],
                min_child_weight=params["min_child_weight"],
                random_state=42,
                n_jobs=-1
            )

            # Entrenamiento
            xgb.fit(X_train, y_train)

            # Predicción
            y_test_pred_xgb = xgb.predict(X_test)

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
                "ablation": nombre,
                "max_depth": params["max_depth"],
                "learning_rate": params["learning_rate"],
                "min_child_weight": params["min_child_weight"],
                "fold": fold,
                "spearman": test_spearman_xgb,
                "r2": test_r2_xgb,
                "rmse": test_rmse_xgb
            })

            pd.DataFrame(resultados_xgb).to_csv(
            "results/xgb_tuning_results.csv",
            index=False
            )