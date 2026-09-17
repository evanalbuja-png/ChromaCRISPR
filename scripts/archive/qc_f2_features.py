import pandas as pd
import numpy as np


INPUT = "data/interim/features/sgRNA_with_f2.csv"
OUTPUT = "logs/f2_features_qc_report.md"


print("Loading dataset...")

df = pd.read_csv(
    INPUT,
    low_memory=False
)


features = [
    "ATAC_mean",
    "ATAC_max",
    "ATAC_p90",
    "ATAC_sum"
]


with open(OUTPUT, "w") as report:

    report.write("# ChromaCRISPR Phase 1 - F2 ATAC Feature QC\n\n")
    report.write("Dataset: sgRNA_with_f2.csv\n\n")

    report.write(f"Total guides: {len(df)}\n\n")


    report.write("## Feature distributions\n\n")

    report.write(
        "| Feature | Mean | Median | Min | Max | Std |\n"
    )
    report.write(
        "|---|---:|---:|---:|---:|---:|\n"
    )

    for feature in features:

        values = df[feature].dropna()

        report.write(
            f"| {feature} | "
            f"{values.mean():.4f} | "
            f"{values.median():.4f} | "
            f"{values.min():.4f} | "
            f"{values.max():.4f} | "
            f"{values.std():.4f} |\n"
        )


    report.write("\n## Zero / low signal percentage\n\n")

    report.write(
        "| Feature | Zero (%) | <1 (%) | <5 (%) |\n"
    )
    report.write(
        "|---|---:|---:|---:|\n"
    )

    for feature in features:

        values = df[feature].dropna()

        zero = (values == 0).mean() * 100
        low1 = (values < 1).mean() * 100
        low5 = (values < 5).mean() * 100

        report.write(
            f"| {feature} | "
            f"{zero:.3f} | "
            f"{low1:.3f} | "
            f"{low5:.3f} |\n"
        )


    if "dataset" in df.columns:

        report.write(
            "\n## Feature distribution by dataset\n\n"
        )

        for dataset, group in df.groupby("dataset"):

            report.write(
                f"### {dataset}\n\n"
            )

            report.write(
                "| Feature | Mean | Median |\n"
            )
            report.write(
                "|---|---:|---:|\n"
            )

            for feature in features:

                values = group[feature].dropna()

                report.write(
                    f"| {feature} | "
                    f"{values.mean():.4f} | "
                    f"{values.median():.4f} |\n"
                )

            report.write("\n")


print("Saved:")
print(OUTPUT)

print("\nSummary:")

for feature in features:

    values = df[feature].dropna()

    print(
        feature,
        "mean=",
        round(values.mean(),4),
        "median=",
        round(values.median(),4),
        "zero%=",
        round((values==0).mean()*100,3)
    )
