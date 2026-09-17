# Semana R — Intento final CRISPOR RS2 / DeepCRISPR (cerrado, bloqueante confirmado)

**Fecha:** 2026-09-13

## Resumen
Segundo intento, en entorno conda aislado (chromacrispr-bench, Python 3.8) 
separado de chromacrispr-phase1, para no arriesgar el pipeline principal.

## Rule Set 2
No se logró ejecutar — ni vía azimuth en el entorno nuevo (Python 3.8) ni vía 
extracción standalone de coeficientes publicados.

## DeepCRISPR
No se logró ejecutar — sin reimplementación moderna mantenida disponible, y 
sin Docker viable en el entorno de trabajo para el contenedor legacy oficial.

## Decisión final
Ambos quedan cerrados como **bloqueantes de infraestructura definitivos** para 
esta fase del proyecto. Se agotaron: (1) instalación directa, (2) entorno 
aislado con versión de Python más compatible, (3) búsqueda de alternativas 
mantenidas por la comunidad, (4) contenedor oficial. No se reabre esta línea 
sin un cambio real de entorno de cómputo (acceso a Docker, o una máquina con 
Python 3.6 legacy disponible).

## Benchmarking Sec. 12 — estado final
Comparación de ChromaCRISPR queda establecida contra: sequence-only (F1), 
Rule Set 3 (sustituto documentado de RS2), y position-baseline. CRISPOR nativo 
y DeepCRISPR quedan fuera, con esta limitación explícita en cualquier reporte 
o manuscrito derivado.
