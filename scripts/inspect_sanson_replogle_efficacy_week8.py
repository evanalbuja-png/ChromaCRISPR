#!/usr/bin/env python3

import pandas as pd

FILES = {
    "Sanson2018": "data/raw/crispr_datasets/Sanson2018/41467_2018_7901_MOESM6_ESM.xlsx",
    "Replogle2022": "data/raw/crispr_datasets/Replogle2022/NIHMS1812939-supplement-11.xlsx",
}

for dataset, path in FILES.items():
    print("\n" + "=" * 80)
    print(dataset)
    print("=" * 80)

    xl = pd.ExcelFile(path)

    for sheet in xl.sheet_names:
        print("\n--- SHEET:", sheet, "---")

        df = pd.read_excel(path, sheet_name=sheet, header=None)

        print("Shape:", df.shape)

        # Mostrar primeras 15 filas completas
        print("\nPrimeras 15 filas:")
        print(df.head(15).to_string(index=False, header=False))

        print("\nFilas que contienen palabras relacionadas con score/eficacia:")
        mask = df.astype(str).apply(
            lambda col: col.str.contains(
                r"score|efficacy|activity|effect|log2|fold|fitness|depletion|enrichment|read|count|phenotype",
                case=False,
                regex=True,
                na=False
            )
        ).any(axis=1)

        hits = df[mask]

        if len(hits):
            print(hits.head(30).to_string(index=False, header=False))
        else:
            print("NONE")

