import pandas as pd

# ============================================================
# Archivos
# ============================================================

working_file = "data/processed/sgRNA_feature_matrix_phase1_with_efficacy_horlbeck.csv"
original_file = "data/processed/temporal/sgRNA_feature_matrix_phase1_with_efficacy.csv"

columns = [
    "sgrna_id",
    "guide_sequence",
    "coordinate",
    "region_id",
    "nearest_tss_gene",
    "experiment",
    "efficacy_score"
]

# ============================================================
# 1. Leer la matriz de trabajo
#    Aquí ya tenemos solamente Horlbeck con efficacy_score
# ============================================================

df = pd.read_csv(
    working_file,
    usecols=[
        "sgrna_id",
        "guide_sequence",
        "coordinate",
        "region_id",
        "nearest_tss_gene",
        "efficacy_score"
    ]
)

# ============================================================
# 2. Identificar secuencias repetidas
# ============================================================

counts = df["guide_sequence"].value_counts()

duplicated_sequences = counts[counts > 1].index

df_duplicadas = df[
    df["guide_sequence"].isin(duplicated_sequences)
].copy()

# ============================================================
# 3. Quedarnos con secuencias que tienen exactamente
#    dos registros y scores diferentes
# ============================================================

score_counts = (
    df_duplicadas
    .groupby("guide_sequence")["efficacy_score"]
    .nunique()
)

duplicados_scores_distintos = score_counts[
    score_counts > 1
].index

df_duplicados_distintos = df_duplicadas[
    df_duplicadas["guide_sequence"].isin(
        duplicados_scores_distintos
    )
].copy()

# ============================================================
# 4. Identificar los 48 candidatos:
#    mismo coordinate en los dos registros
# ============================================================

coordinate_counts = (
    df_duplicados_distintos
    .groupby("guide_sequence")["coordinate"]
    .nunique()
)

replicas_candidatas = coordinate_counts[
    coordinate_counts == 1
].index

df_48 = df_duplicados_distintos[
    df_duplicados_distintos["guide_sequence"].isin(
        replicas_candidatas
    )
].copy()

# ============================================================
# 5. Leer el archivo original por bloques para recuperar
#    la columna "experiment"
# ============================================================

results = []

for chunk in pd.read_csv(
    original_file,
    usecols=columns,
    chunksize=10000
):
    subset = chunk[
        chunk["guide_sequence"].isin(replicas_candidatas)
    ]

    if not subset.empty:
        results.append(subset)

df_48_original = pd.concat(
    results,
    ignore_index=True
)

# ============================================================
# 6. Mostrar resultado
# ============================================================

print("============================================================")
print("48 CANDIDATOS A REPLICAS")
print("============================================================")

print(f"Filas: {len(df_48_original)}")
print(
    f"Secuencias: "
    f"{df_48_original['guide_sequence'].nunique()}"
)

print("\nNúmero de experimentos por secuencia:")

experiment_counts = (
    df_48_original
    .groupby("guide_sequence")["experiment"]
    .nunique()
)

print(experiment_counts.value_counts())

# ============================================================
# 7. Guardar para inspección
# ============================================================

output_file = (
    "data/processed/temporal/"
    "horlbeck_48_replicates_with_experiment.csv"
)

df_48_original.to_csv(
    output_file,
    index=False
)

print("\nResultado guardado en:")
print(output_file)
