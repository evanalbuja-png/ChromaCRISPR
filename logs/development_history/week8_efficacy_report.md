# ChromaCRISPR Phase 1 — Week 8

## Paso 5 — Identificación de variable de eficacia

- Matriz base: `data/processed/sgRNA_feature_matrix_phase1.csv`
- Guías: 155,454

### Resultado

Se inspeccionaron los archivos tabulares disponibles en los directorios de datos originales/intermedios.

Se identificaron candidatos de columnas de eficacia y de llaves de unión, pero **no se realizó una integración automática en esta ejecución** para evitar asignar incorrectamente un score a una guía.

### Candidatos detectados

#### sanson2018
- Archivo: `data/interim/sanson2018_brunello_reference.csv`
- Targets candidatos: ``
- Llaves candidatas: `guide_sequence`

#### sanson2018
- Archivo: `data/interim/sanson2018_recovered_official.csv`
- Targets candidatos: `score_raw;score_norm;score_type`
- Llaves candidatas: `guide_sequence;guide_length`

#### sanson2018
- Archivo: `data/interim/sanson_features_reference.tsv`
- Targets candidatos: ``
- Llaves candidatas: `guide_sequence;guide_position`

#### Gasperini2019
- Archivo: `data/interim/gasperini_coordinate_validation_sample.csv`
- Targets candidatos: ``
- Llaves candidatas: `guide_sequence`

#### Gasperini2019
- Archivo: `data/interim/gasperini_validation_window10kb_ensembl_intersect.txt`
- Targets candidatos: ``
- Llaves candidatas: ``

### Estado

- Matriz `sgRNA_feature_matrix_phase1.csv`: no modificada.
- `efficacy_score`: todavía no incorporado.
- Se requiere verificar la correspondencia exacta entre cada dataset, su score y la llave `guide_sequence`/ID antes de efectuar el merge.
