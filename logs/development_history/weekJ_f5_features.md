# Semana J — F5 Genomic Context (cerrada, parcial)

**Fecha:** 2026-09-12

## Conteo n=12 vs texto
- Lectura expandida justificada: dist target TSS + dist expressed TSS + 4 clases cCRE + PhyloP (+ GERP + JASPAR) ≈ 9–11.
- n=12 declarado: discrepancia residual documentada.
- Oficiales implementadas ahora: **7** (2 dist + 4 cCRE + PhyloP).
- Pendientes: GERP++, JASPAR motif count.

## Fuentes
| Track | Accession / URL |
|---|---|
| cCRE | ENCSR800VNX → ENCFF420VPZ |
| K562 gene quant | ENCSR615EEK → ENCFF421TJX |
| PhyloP 100way | UCSC hg38.phyloP100way.bw |
| TSS | GENCODE v49 |

## PAM context (F1 pendiente)
- Extraído 6mer + one-hot 4 posiciones para 18.318 guías con strand (Horlbeck).
- Gasperini: sin strand en reference → sin PAM.
- Secuencias ejemplo no muestran NGG evidente → **convención de coordenada no validada**.
- **Decisión: PAM real de F1 sigue PENDIENTE** hasta verificación manual de 2–3 loci.

## Imputación Sec. 7.3
- `log_dist_target_tss`: mediana + `log_dist_target_tss_missing` (60)
- `phylop100_mean_pm100`: mediana + `phylop100_missing` (2)

## Archivo
`data/interim/F5_genomic_context_features_v1.csv` (19.446 filas)

## Corrección post-diagnóstico (2026-09-12)

### NaN
- 2256 NaN = solo `strand` (1128) + antiguo `pam_seq_6mer` (1128) de Gasperini.
- No había NaN ocultos en cCRE/PhyloP/dists. Continuas ya imputadas.

### PAM real (F1) — VALIDADO
- Bug: offset desde `coordinate` sin anclar la guía al genoma.
- Fix: match exacto de guía o RC en ±80 bp; PAM = 3 nt 3′ en hebra del match.
- 3/3 controles NGG; 99.99% NGG en 18.433 recuperadas; 1.013 sin match en ventana.
- Columnas: `pam_3nt`, `pam_is_NGG`, `pam_pos{1-3}_{ACGT}` (12 one-hot).
- **F1 PAM context real: CERRADO** (con cobertura parcial 18.4k/19.4k).

## Fix NaN propagation PAM (mismo patrón que compartment_A en Semana I)
- 1.013 guías sin match: `pam_is_NGG` y 12 one-hot pasaron de 0/False implícito → **NaN**.
- Rate NGG solo sobre matched (skipna) debe ser ≥99.9%.
- Regla general: features derivadas condicionalmente nunca deben codificar "ausencia de cálculo" como 0/False.
