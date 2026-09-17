# Semana B — Integración DS3 Gasperini 2019 (cerrada)

**Fecha:** 2026-09-09  
**Objetivo:** Derivar efficacy scores a nivel de guía individual desde Gasperini et al. 2019 (Cell) cumpliendo el espíritu de Sec. 7.1, con polaridad corregida para armonización futura.

## Archivos usados
- `data/raw/crispr_datasets/Gasperini2019/NIHMS1038673-supplement-TableS2.xlsx`
  - Hoja `Key` → definición oficial de columnas
  - Hoja `B_AtScale_664_enhancergenepairs` → 664 pares + fold_change + high_confidence
  - Hoja `S2A_AtScale_library_gRNA.cs` → librería de guías
- `data/processed/gasperini2019_reference_FINAL.csv` → coordenadas hg38 (limitación heredada)

## Qué representa el score original
Según la hoja Key del suplementario:
> `Diff_expression_test_fold_change` = **Proportion remaining transcript of target gene**

No es log2FC ni p-value. Es la fracción de transcrito que permanece tras CRISPRi (0–1).  
Valores más bajos = mayor represión = guía más efectiva.

## Transformación de polaridad (Opción A)
Se aplicó:
effect_size_log2 = -log2(proportion_remaining + 0.001)
- Dirección corregida: valor más alto = mayor eficacia (consistente con Horlbeck).
- Familia logarítmica alineada conceptualmente con “log2 fold-change” de Sec. 7.1.
- Spearman(proportion_remaining, effect_size_log2) = −1.000000 (verificado).

**Nota histórica:** En Phase 1 (week8/week9) ya se había decidido no usar esta columna como target de entrenamiento precisamente porque no era equivalente al score de Horlbeck. La transformación actual resuelve la polaridad para permitir su uso controlado en la armonización rank-based de Semana D.

## Desviación documentada – criterio ≥50 células
Sec. 7.1 exige: “guides with ≥ 50 cells perturbed and a clear expression readout for the proximal gene”.  
Las tablas suplementarias públicas **no contienen** columna de n_cells por guía.  
Se usó `high_confidence_subset = TRUE` del paper (470 pares → 1.328 guías) como proxy. Este filtro combina significancia empírica ajustada + criterios internos de calidad del estudio (incluyendo cobertura).  
**Desviación explícita**, mismo estándar de documentación que Jurkat/MCF7 → HT29/A375.

## Limitación heredada – coordenadas hg38
`gasperini2019_reference_FINAL.csv` (creado 18 Jul 2026) aporta las coordenadas hg38.  
No existe log del pipeline de mapeo Bowtie2 zero-mismatch ni estadísticas de unique/multimap.  
Mismo estatus de “limitación heredada” que Horlbeck en Semana A.  
Overlap con S2A: ~12.3k / 13.2k (≈900 guías sin match exacto; no bloquea).

## Cobertura final (set recomendado)
- Archivo: `data/interim/DS3_gasperini/DS3_gasperini_guides_efficacy_highconf_polarity_fixed.csv`
- Guías high-confidence con effect_size_log2 + coordenada hg38: **1.328**
- Rango effect_size_log2: 0.117 – 5.272 (media 0.620)
- Unique Target_Sites: ~470 (high-confidence)

## Decisión de uso
- Usar **high-confidence + effect_size_log2** como fuente de entrenamiento adicional en Semana D.
- No se integra todavía a la matriz consolidada.
- No se modificó ningún archivo de Horlbeck ni las matrices FIXED.

## Criterio de éxito
CUMPLIDO — score derivado, polaridad corregida, desviaciones y limitaciones documentadas explícitamente.
