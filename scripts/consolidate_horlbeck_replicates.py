import pandas as pd

input_file = "data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck.csv"
output_file = "data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck_consolidated.csv"

df = pd.read_csv(input_file)

# Columnas que identifican el mismo locus
group_cols = [
    "guide_sequence",
    "coordinate",
    "region_id",
    "nearest_tss_gene"
]

# Promediar efficacy_score únicamente cuando existe
# más de un registro para la misma guía/locus.
df_consolidated = (
    df.groupby(group_cols, as_index=False)
      .agg({
          **{
              col: "first"
              for col in df.columns
              if col not in group_cols + ["efficacy_score"]
          },
          "efficacy_score": "mean"
      })
)

print("Filas originales:", len(df))
print("Filas consolidadas:", len(df_consolidated))
print("Filas eliminadas:", len(df) - len(df_consolidated))

df_consolidated.to_csv(
    output_file,
    index=False
)

print("\nGuardado en:")
print(output_file)
