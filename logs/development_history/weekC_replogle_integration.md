# Semana C — Integración DS4 Replogle 2022 (cerrada — Opción 2)

**Fecha:** 2026-09-10  
**Objetivo:** Derivar efficacy scores a nivel de guía individual desde Replogle et al. 2022 (Cell) cumpliendo Sec. 7.1.

## Archivos presentes en el proyecto
- `data/raw/crispr_datasets/Replogle2022/NIHMS1812939-supplement-11.xlsx`  
  → Solo tablas de diseño de librería (pares de guías). Sin scores, sin n_cells, sin log2FC.
- `data/processed/replogle2022_reference_FINAL.csv`  
  → 22.576 filas de coordenadas hg38 (anotación). Limitación heredada de mapeo (mismo estatus que Horlbeck/Gasperini).

## Búsqueda acotada de tabla de scores (límite estricto)
Se buscó:
1. Data Availability del paper → SRA PRJNA831566 + portal https://gwps.wi.mit.edu + Figshare.
2. Figshare oficial:
   - “commonly requested supplemental files” (503 MB): Anderson-Darling p-values de DE, embeddings, datos Z-norm para figuras. No contiene un score de knockdown del gen target por guía listo para usar.
   - Processed AnnData (≈160 GB) y MTX (≈97 GB): demasiado grandes; requerirían re-derivar el effect size.
3. Supplement-11 ya descargado: solo librería.
4. No hay menciones en logs del proyecto de haber descargado scores de Replogle.

**Resultado:** No se identificó una tabla oficial, pequeña (<500 MB) y lista de “guía → effect size del gen target + n_cells + p-adj”.  
Se aplica el límite acordado → **Opción 2**.

## Decisión
- DS4 (Replogle 2022) **queda fuera** de la harmonización de Semana D.
- Se conserva el reference de coordenadas por si se retoma más adelante.
- No se descarga ningún objeto grande sin confirmación explícita posterior.
- No se inventa ni se usa fuente no oficial (TruthSeq u otras).

## Criterio de éxito
CUMPLIDO en la modalidad “documentación inequívoca de por qué no es viable en este entorno”.
