import pandas as pd
import numpy as np

df = pd.read_csv("data/interim/DS3_gasperini/gasperini_highconf_efficacy.csv")  # ajustar nombre real

PSEUDOCOUNT = 0.001
df["effect_size_log2"] = -np.log2(df["proportion_remaining"] + PSEUDOCOUNT)

print(f"Rango proportion_remaining: {df['proportion_remaining'].min():.3f} - {df['proportion_remaining'].max():.3f}")
print(f"Rango effect_size_log2: {df['effect_size_log2'].min():.3f} - {df['effect_size_log2'].max():.3f}")
print(f"Correlación (debe ser -1.0 exacta, es transformación monótona): "
      f"{df['proportion_remaining'].corr(df['effect_size_log2'], method='spearman'):.4f}")

df.to_csv("data/interim/DS3_gasperini/gasperini_highconf_efficacy_polarity_fixed.csv", index=False)
print("Guardado con columna effect_size_log2 (mayor = más efectivo, consistente con Horlbeck).")
