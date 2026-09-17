import pandas as pd
from pathlib import Path


input_file = Path(
    "data/raw/sanson2018_addgene_reference/"
    "broadgpp-brunello-library-contents.txt"
)

output_file = Path(
    "data/interim/sanson2018_brunello_reference.csv"
)


# Cargar Brunello
df = pd.read_csv(input_file, sep="\t")


# Conversión RefSeq GRCh38 -> chromosome notation
refseq_to_chr = {
    "NC_000001.11": "chr1",
    "NC_000002.12": "chr2",
    "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",
    "NC_000005.10": "chr5",
    "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",
    "NC_000008.11": "chr8",
    "NC_000009.12": "chr9",
    "NC_000010.11": "chr10",
    "NC_000011.10": "chr11",
    "NC_000012.12": "chr12",
    "NC_000013.11": "chr13",
    "NC_000014.9": "chr14",
    "NC_000015.10": "chr15",
    "NC_000016.10": "chr16",
    "NC_000017.11": "chr17",
    "NC_000018.10": "chr18",
    "NC_000019.10": "chr19",
    "NC_000020.11": "chr20",
    "NC_000021.9": "chr21",
    "NC_000022.11": "chr22",
    "NC_000023.11": "chrX",
    "NC_000024.10": "chrY"
}


out = pd.DataFrame()

out["guide_sequence"] = (
    df["sgRNA Target Sequence"]
    .str.upper()
    .str.strip()
)

out["gene_symbol"] = df["Target Gene Symbol"]

out["gene_id"] = (
    df["Target Gene ID"]
    .astype("Int64")
)

out["chromosome"] = (
    df["Genomic Sequence"]
    .map(refseq_to_chr)
)

out["coordinate"] = (
    df["Position of Base After Cut (1-based)"]
    .astype("Int64")
)

out["strand"] = df["Strand"]

out["genome_build"] = "hg38"

out["source"] = "Addgene"

out["dataset"] = "Sanson2018"

out["library"] = "Brunello"

out["is_non_targeting_control"] = (
    df["Target Gene Symbol"]
    == "Non-Targeting Control"
)

print("Output shape:")
print(out.shape)

print("\nMissing chromosomes:")
print(out["chromosome"].isna().sum())

print("\nChromosomes:")
print(out["chromosome"].value_counts())


output_file.parent.mkdir(
    parents=True,
    exist_ok=True
)

out.to_csv(
    output_file,
    index=False
)

print("\nSaved:")
print(output_file)
