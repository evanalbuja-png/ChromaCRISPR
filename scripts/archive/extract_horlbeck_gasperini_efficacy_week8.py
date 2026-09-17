#!/usr/bin/env python3

import os
import pandas as pd

RAW = "data/raw/crispr_datasets"

OUT_HORLBECK = "data/interim/horlbeck2016_efficacy.csv"
OUT_GASPERINI = "data/interim/gasperini2019_efficacy.csv"
REPORT = "logs/week8_horlbeck_gasperini_efficacy.md"

os.makedirs("data/interim", exist_ok=True)
os.makedirs("logs", exist_ok=True)

# ============================================================
# HORLBECK 2016
# ============================================================

horlbeck_file = os.path.join(
    RAW, "Horlbeck2016", "elife-19760-supp1-v2.xlsx"
)

print("=" * 70)
print("HORLBECK2016")
print("=" * 70)

h_crispri = pd.read_excel(
    horlbeck_file,
    sheet_name="CRISPRi"
)

h_crispra = pd.read_excel(
    horlbeck_file,
    sheet_name="CRISPRa"
)

h_crispri = h_crispri.rename(columns={
    "sgRNA sequence": "guide_sequence",
    "CRISPRi activity score [Horlbeck et al., eLife 2016]":
        "efficacy_score"
})

h_crispri["dataset"] = "Horlbeck2016"
h_crispri["experiment"] = "CRISPRi"

h_crispri = h_crispri[
    [
        "guide_sequence",
        "gene symbol",
        "chromosome",
        "PAM genomic coordinate [hg19]",
        "strand targeted",
        "efficacy_score",
        "dataset",
        "experiment"
    ]
]

h_crispra = h_crispra.rename(columns={
    "sgRNA sequence": "guide_sequence",
    "CRISPRa activity score":
        "efficacy_score"
})

h_crispra["dataset"] = "Horlbeck2016"
h_crispra["experiment"] = "CRISPRa"

h_crispra = h_crispra[
    [
        "guide_sequence",
        "gene symbol",
        "chromosome",
        "PAM genomic coordinate [hg19]",
        "strand targeted",
        "efficacy_score",
        "dataset",
        "experiment"
    ]
]

horlbeck = pd.concat(
    [h_crispri, h_crispra],
    ignore_index=True
)

horlbeck["guide_sequence"] = (
    horlbeck["guide_sequence"]
    .astype(str)
    .str.upper()
    .str.strip()
)

horlbeck["efficacy_score"] = pd.to_numeric(
    horlbeck["efficacy_score"],
    errors="coerce"
)

horlbeck = horlbeck.dropna(
    subset=["guide_sequence", "efficacy_score"]
)

horlbeck = horlbeck.drop_duplicates(
    subset=["guide_sequence", "experiment"],
    keep="first"
)

horlbeck.to_csv(
    OUT_HORLBECK,
    index=False
)

print("CRISPRi rows:", len(h_crispri))
print("CRISPRa rows:", len(h_crispra))
print("Combined rows:", len(horlbeck))
print("Unique guides:", horlbeck["guide_sequence"].nunique())

print("\nScore distribution:")
print(horlbeck["efficacy_score"].describe())


# ============================================================
# GASPERINI 2019
# ============================================================

print("\n" + "=" * 70)
print("GASPERINI2019")
print("=" * 70)

s1_file = os.path.join(
    RAW, "Gasperini2019",
    "NIHMS1038673-supplement-TableS1.xlsx"
)

s2_file = os.path.join(
    RAW, "Gasperini2019",
    "NIHMS1038673-supplement-TableS2.xlsx"
)

# -----------------------------
# S1
# -----------------------------

s1_guides = pd.read_excel(
    s1_file,
    sheet_name="A_Pilot_gRNA_library"
)

s1_scores = pd.read_excel(
    s1_file,
    sheet_name="B_Pilot_145_enhancergene_pairs"
)

s1_guides = s1_guides.rename(columns={
    "Spacer": "guide_sequence"
})

s1_scores = s1_scores.rename(columns={
    "Diff_expression_test_fold_change":
        "efficacy_score"
})

s1 = s1_guides.merge(
    s1_scores[
        [
            "Target_Site",
            "efficacy_score"
        ]
    ],
    on="Target_Site",
    how="left"
)

s1["dataset"] = "Gasperini2019"
s1["experiment"] = "S1"

# -----------------------------
# S2
# -----------------------------

s2_guides = pd.read_excel(
    s2_file,
    sheet_name="S2A_AtScale_library_gRNA.cs"
)

s2_scores = pd.read_excel(
    s2_file,
    sheet_name="B_AtScale_664_enhancergenepairs"
)

s2_guides = s2_guides.rename(columns={
    "Spacer": "guide_sequence"
})

s2_scores = s2_scores.rename(columns={
    "Diff_expression_test_fold_change":
        "efficacy_score"
})

s2 = s2_guides.merge(
    s2_scores[
        [
            "Target_Site",
            "efficacy_score"
        ]
    ],
    on="Target_Site",
    how="left"
)

s2["dataset"] = "Gasperini2019"
s2["experiment"] = "S2"

# -----------------------------
# Combine
# -----------------------------

gasperini = pd.concat(
    [s1, s2],
    ignore_index=True
)

gasperini["guide_sequence"] = (
    gasperini["guide_sequence"]
    .astype(str)
    .str.upper()
    .str.strip()
)

gasperini["efficacy_score"] = pd.to_numeric(
    gasperini["efficacy_score"],
    errors="coerce"
)

gasperini_valid = gasperini.dropna(
    subset=["efficacy_score"]
).copy()

gasperini_valid = gasperini_valid.drop_duplicates(
    subset=["guide_sequence", "experiment"],
    keep="first"
)

gasperini_valid.to_csv(
    OUT_GASPERINI,
    index=False
)

print("S1 guides:", len(s1_guides))
print("S1 scored guide rows:", s1["efficacy_score"].notna().sum())

print("S2 guides:", len(s2_guides))
print("S2 scored guide rows:", s2["efficacy_score"].notna().sum())

print("Combined scored rows:", len(gasperini_valid))
print("Unique guides:", gasperini_valid["guide_sequence"].nunique())

print("\nScore distribution:")
print(gasperini_valid["efficacy_score"].describe())


# ============================================================
# REPORT
# ============================================================

with open(REPORT, "w", encoding="utf-8") as f:

    f.write("# ChromaCRISPR Phase 1 — Week 8\n\n")
    f.write("## Extracción de eficacia — Horlbeck2016 y Gasperini2019\n\n")

    f.write("## Horlbeck2016\n\n")
    f.write(f"- Filas CRISPRi: {len(h_crispri)}\n")
    f.write(f"- Filas CRISPRa: {len(h_crispra)}\n")
    f.write(f"- Filas combinadas: {len(horlbeck)}\n")
    f.write(
        f"- Guías únicas: "
        f"{horlbeck['guide_sequence'].nunique()}\n"
    )
    f.write(
        "- Variable: `CRISPRi activity score` / "
        "`CRISPRa activity score`\n"
    )
    f.write(
        "- Normalización: no aplicada; se conservan "
        "los scores originales.\n\n"
    )

    f.write("## Gasperini2019\n\n")
    f.write(
        f"- S1 guías: {len(s1_guides)}\n"
    )
    f.write(
        f"- S1 guías con score: "
        f"{s1['efficacy_score'].notna().sum()}\n"
    )
    f.write(
        f"- S2 guías: {len(s2_guides)}\n"
    )
    f.write(
        f"- S2 guías con score: "
        f"{s2['efficacy_score'].notna().sum()}\n"
    )
    f.write(
        f"- Filas combinadas con score: "
        f"{len(gasperini_valid)}\n"
    )
    f.write(
        f"- Guías únicas con score: "
        f"{gasperini_valid['guide_sequence'].nunique()}\n"
    )
    f.write(
        "- Variable: `Diff_expression_test_fold_change`\n"
    )
    f.write(
        "- La asociación se realizó mediante "
        "`Spacer → Target_Site → efficacy_score`.\n"
    )

print("\n" + "=" * 70)
print("REPORT:", REPORT)
print("HORLBECK:", OUT_HORLBECK)
print("GASPERINI:", OUT_GASPERINI)
print("=" * 70)

