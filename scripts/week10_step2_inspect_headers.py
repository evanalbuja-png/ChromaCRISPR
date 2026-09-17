import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)

path = "data/raw/crispr_datasets/Sanson2018/41467_2018_7901_MOESM6_ESM.xlsx"

for sheet in ["SetA raw reads", "SetB raw reads"]:
    print(f"\n{'='*80}\n{sheet}\n{'='*80}")
    df = pd.read_excel(path, sheet_name=sheet, header=None, nrows=10)
    print(df.to_string())
    full = pd.read_excel(path, sheet_name=sheet, header=None)
    print(f"\nShape completo de la hoja: {full.shape}")

for sheet in ["SetA sgRNA annotations", "SetB sgRNA annotations"]:
    print(f"\n{'='*80}\n{sheet}\n{'='*80}")
    df = pd.read_excel(path, sheet_name=sheet)
    print(df.shape)
    print(df.head(3))
