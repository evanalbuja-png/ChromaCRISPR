import pandas as pd

# Importar dataframe
df = pd.read_csv("data/processed/temporal/sgRNA_feature_matrix_phase1_with_efficacy.csv")

# Extraer matriz
horlbeck_df = df[df['dataset'] == 'Horlbeck2016']
horlbeck_df_final = horlbeck_df[horlbeck_df["efficacy_score"].notna()]
# Contar cuantos NaN hay en el efficacy_score

nan_counts = horlbeck_df_final['efficacy_score'].isna().sum()
print(f"Cantidad de valores NaN en efficacy_score: {nan_counts}")
# Porcentaje
nan_percentage = (nan_counts / len(horlbeck_df_final)) * 100
print(f"Porcentaje de valores NaN en efficacy_score: {nan_percentage:.2f}")

# Verificar si gene_symbol y strand son lo mismo

gene_nan = horlbeck_df_final["gene_symbol"].isna()
strand_nan = horlbeck_df_final["strand"].isna()
gene_nan.equals(strand_nan)
print(f"¿Los valores NaN en gene_symbol y strand son iguales? {gene_nan.equals(strand_nan)}")

# Eliminar columnas innecesarias
horlbeck_df_final = horlbeck_df_final.drop(
    columns=["gene_id", "experiment", "gene_symbol", "strand"]
)

# Guardar el dataframe final
horlbeck_df_final.to_csv(
    "data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck.csv",
    index=False
)