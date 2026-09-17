## Tarea 1 - Recuperación oficial Sanson2018 Brunello (Addgene)

Fecha: 2026-07-18

### Fuente oficial
- Addgene plasmid:
  - #73179 Human sgRNA library Brunello in lentiCRISPRv2
- Archivo:
  - broadgpp-brunello-library-contents.txt
- Ubicación:
  - data/raw/sanson2018_addgene_reference/

### Validación del archivo
- Filas:
  - 77,441
- Columnas:
  - 11
- Campo de guía:
  - sgRNA Target Sequence
- Longitud:
  - 20 nt
- Genome build detectado:
  - GRCh38 / hg38

### Join preliminar por secuencia
Método:
- Match exacto entre:
  - Brunello: sgRNA Target Sequence
  - sgRNA_unified: guide_sequence

Resultados:
- Unified guides:
  - 183,723
- Brunello guides:
  - 77,441
- Matches:
  - 2,142

Métricas:
- Match rate sobre unified:
  - 1.17%
- Cobertura Brunello:
  - 2.68%

### Interpretación
El join directo por secuencia recupera únicamente una fracción pequeña de Brunello.
Se confirma la necesidad de utilizar la referencia oficial Addgene como fuente primaria para recuperación de Sanson2018.

## Tarea 2 - Normalización referencia Sanson2018 Brunello

Fecha: 2026-07-18

### Objetivo
Generar una referencia estandarizada de Brunello para integración con sgRNA_unified.

### Archivo generado
- data/interim/sanson2018_brunello_reference.csv

### Transformaciones realizadas
- Renombrado de campos a esquema unificado.
- Conversión de accesiones RefSeq GRCh38:
  - NC_000001.11 → chr1
  - ...
  - NC_000024.10 → chrY
- Normalización de secuencias sgRNA:
  - uppercase
  - longitud 20 nt

### Resultado
- Guías totales:
  - 77,441

- Guías con coordenadas hg38:
  - 76,441

- Controles non-targeting:
  - 1,000

### Genome build
- GRCh38 / hg38

### Decisiones
Los controles Non-Targeting Control se mantienen sin coordenadas genómicas y se marcan mediante:
- is_non_targeting_control = True

## Tarea 3 - Diagnóstico recuperación Sanson2018 en sgRNA_unified

Fecha: 2026-07-18

### Método
Join exacto por secuencia entre:
- Brunello referencia Addgene normalizada
- data/interim/sgRNA_unified.csv

Se excluyeron controles Non-Targeting Control.

### Resultados
- Brunello genómico:
  - 76,441 guías

- Matches:
  - 1,150 registros
  - 1,085 guías únicas

### Distribución de datasets encontrados
- Sanson2018:
  - 1,033

- Replogle2022:
  - 83

- Horlbeck2016:
  - 34

### Observaciones
- El join por secuencia recupera una fracción pequeña del catálogo Brunello.
- La secuencia no es suficiente para resolver procedencia porque existen guías compartidas entre datasets.
- Brunello Addgene está en GRCh38/hg38.
- Las coordenadas de sgRNA_unified no son directamente comparables debido a diferencias de genome build (principalmente hg19).

### Decisión
La referencia Addgene Brunello será considerada fuente de verdad para recuperación Sanson2018.

## Tarea 4 - Evaluación de conflictos Brunello vs sgRNA_unified

Fecha: 2026-07-18

### Comparación
Join por secuencia entre:
- sanson2018_brunello_reference.csv
- sgRNA_unified.csv

### Resultados

Matches:
- 1150 registros
- 1085 guías únicas

Concordancia:
- gene_symbol:
  - 91.6%

- chromosome:
  - 3.0%

Genome build:
- Brunello:
  - 100% hg38

- Unified:
  - 1067 hg19
  - 83 hg38

Duplicados:
- Brunello:
  - 0 secuencias duplicadas

- Unified matches:
  - 65 secuencias duplicadas por presencia en múltiples datasets

### Decisión
La referencia Addgene Brunello será utilizada como fuente de verdad para Sanson2018.
La recuperación no se realizará mediante secuencia aislada.
La armonización hg19/hg38 será necesaria antes de integrar coordenadas.

## Tarea 5 - Generación tabla oficial Sanson2018 recuperada

Fecha: 2026-07-18

### Archivo generado

data/interim/sanson2018_recovered_official.csv

### Resultado

- Registros:
  - 77,441

- Guías genómicas:
  - 76,441

- Non-Targeting Controls:
  - 1,000

### Validación

- Genome build:
  - hg38 (77,441 registros)

- Longitud guías:
  - 20 nt

- Secuencias duplicadas:
  - 0

### Campos faltantes

Los 1,000 controles Non-Targeting no poseen:
- gene_id
- chromosome
- coordinate
- strand

Los scores originales no fueron recuperados en esta tabla porque la referencia intermedia normalizada no conserva Rule Set 2 score.

### Decisión

La tabla oficial Brunello queda preparada como fuente primaria para Sanson2018.
No se modifica todavía sgRNA_unified.csv.

## Tarea 6 - Estimación cobertura recuperación Sanson2018

Fecha: 2026-07-18

### Comparación

Fuente actual:
- sgRNA_unified.csv

Referencia:
- sanson2018_recovered_official.csv

### Resultados

Unified actual:
- filas:
  - 183,723

- guías únicas:
  - 167,470


Brunello oficial:
- total:
  - 77,441

- guías genómicas:
  - 76,441


Solapamiento por secuencia:
- presentes actualmente:
  - 1,085

- nuevas guías esperadas:
  - 75,356


### Conclusión

La recuperación oficial incrementará significativamente la cobertura de Sanson2018.

No se realizará reemplazo directo de registros existentes.
La estrategia será integración controlada preservando trazabilidad de fuente.

## Tarea 7-8: Recuperación e integración oficial Sanson2018 Brunello

Fecha: 2026-07-18

Fuente:
- Addgene plasmid #73179
- broadgpp-brunello-library-contents.txt

Resultados:
- Total Brunello:
  - 77,441 sgRNAs
- Guías genómicas:
  - 76,441
- Non-targeting controls:
  - 1,000

Integración:
- Archivo generado:
  - data/interim/sgRNA_unified_with_sanson2018_v2.csv

Antes:
- 183,723 registros

Después:
- 259,079 registros

Nuevos registros:
- 75,356

Trazabilidad:
- source=Addgene_Brunello

Validación:
- Los nuevos registros no generan duplicados dentro de Sanson2018.
- Los duplicados guide_sequence+dataset observados provienen exclusivamente de datasets previamente existentes:
  - Replogle2022
  - Gasperini2019
  - Horlbeck2016

Conclusión:
La recuperación oficial de Brunello aumenta significativamente la cobertura de Sanson2018 y queda lista para continuar con la estandarización de datasets.

## Tarea 9: Evaluación genome build tras integración Sanson2018

Archivo:
- data/interim/sgRNA_unified_with_sanson2018_v2.csv

Distribución:
- total:
  - 259079

Genome build:
- hg19:
  - 151177
- hg38:
  - 107902

Por dataset:

Gasperini2019:
- hg19: 16307

Horlbeck2016:
- hg19: 20809

Replogle2022:
- hg38: 32546

Sanson2018:
- hg19 original: 114061
- hg38 Addgene Brunello: 75356

Coordenadas faltantes:
- chromosome:
  - 148361
- coordinate:
  - 148361

Conclusión:
La recuperación Brunello fue integrada correctamente.
La siguiente etapa será normalización de genome build y recuperación de coordenadas faltantes.

## Tarea 10: Diagnóstico Sanson2018 genome build

Resultado:

Sanson2018 total:
- 189417 registros

Fuentes:
- existing_unified:
  - 114061
- Addgene_Brunello:
  - 75356

Genome build:
- hg19:
  - 114061
- hg38:
  - 75356

Coordenadas:
- Addgene_Brunello:
  - coordenadas completas
- existing_unified:
  - sin chromosome/coordinate

Conclusión:
Los registros antiguos de Sanson2018 no pueden ser sometidos a liftover directo.
La recuperación de anotaciones debe realizarse mediante matching contra la referencia oficial Brunello.

## Paso 1 - Diagnóstico Brunello bug

Fecha: 2026-07-18

### Objetivo
Investigar el match rate anómalo entre Sanson2018 en sgRNA_unified.csv y Brunello Addgene oficial.

### Archivos evaluados

Input unified:
- data/interim/sgRNA_unified.csv

Sanson source:
- data/raw/crispr_datasets/Sanson2018/41467_2018_7901_MOESM6_ESM.xlsx

Brunello reference:
- data/raw/sanson2018_addgene_reference/broadgpp-brunello-library-contents.txt

### Hallazgos

sgRNA_unified.csv contiene Sanson2018 generado desde:
- SetA sgRNA annotations (57050 filas)
- SetB sgRNA annotations (57011 filas)

Fuente:
- suplemento MOESM6 de Sanson2018

Estas tablas contienen:
- sgRNA Sequence
- Annotated Gene Symbol
- Annotated Gene ID

No corresponden directamente a la librería Brunello Addgene.

Brunello Addgene contiene:
- 77441 guías oficiales

Comparación inicial:
- Match por secuencia: ~2.6%

Interpretación:
- El bajo match rate no se debe a corrupción de secuencias.
- La comparación inicial mezclaba dos fuentes distintas.
- No se evidencia bug de parsing en scripts/parsers/sanson.py.

Conclusión:
- Paso 1 completado.
- La referencia Brunello Addgene debe tratarse separadamente del dataset Sanson2018 Dolcetto/CRISPRi.

## Corrección crítica: Sanson2018 vs Brunello

Fecha: 2026-07-18

Diagnóstico:
La comparación inicial entre Sanson2018 y Addgene Brunello era incorrecta.

Sanson2018 Supplementary Data:
- Archivo: 41467_2018_7901_MOESM6_ESM.xlsx
- Librería: Dolcetto CRISPRi
- SetA: 57,050 guías únicas
- SetB: 57,011 guías únicas
- Longitud: 20 nt
- Duplicados: 0

Addgene Brunello:
- Librería distinta
- CRISPR knockout Cas9 activo

El bajo match rate (~2.6%) no corresponde a corrupción de secuencias.
La discrepancia se debe a comparación entre librerías diferentes.

Decisión:
No usar Brunello como referencia de recuperación para Sanson2018.

## Final Sanson2018 TSS resolution

Dataset:
- 419 sgRNAs

Resolution:
- GENCODE: 190
- FANTOM5: 174
- unresolved_TSS_review: 55

Total resolved TSS:
364/419 (86.9%)

CRISPRi promoter candidates:
131

Unresolved cases:
55
Main categories:
- OR gene family
- ZNF family
- repetitive/complex loci

Diagnostic steps completed:
- exact hg38 sequence search
- reverse complement search
- ±2 kb TSS window search
- FANTOM5 CAGE integration
- 1 mismatch rescue attempt

Decision:
Remaining unresolved cases retained with explicit metadata
(coordinate_source = unresolved_TSS_review)
and excluded from promoter-level feature extraction.