import pandas as pd

# Importar dataframe
df = pd.read_csv("data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck.csv")

# Comprobar dimensiones y nombres de las columnas del dataframe
print(df.shape)
print(df.columns)

# Conteo de NaN por columna

nan_counts = df.isna().sum()
print(nan_counts)

# Cálculo de porcentaje de NaN por columna
nan_percentage = (nan_counts / len(df)) * 100
print(nan_percentage)
