#!/usr/bin/env python3
"""
Script: run_azimuth_rs2.py
Descripción: Inferencia de scores Rule Set 2 (Azimuth V3) aplicando monkey-patching 
             de módulos legacy de scikit-learn sobre sys.modules.
Entorno: chromacrispr-bench
"""

import sys
import os
import pandas as pd
import numpy as np

# ==============================================================================
# MONKEY-PATCHING DE SCikit-LEARN LEGACY
# ==============================================================================
# Intercepta las llamadas de pickle.load() dentro de Azimuth redireccionando
# los módulos legacy de scikit-learn <= 0.20 a la API actual.
try:
    import sklearn.ensemble._gb as _gb
    import sklearn.ensemble._forest as _forest
    import sklearn.tree._classes as _tree_classes
    
    # Inyección de módulos en sys.modules
    sys.modules['sklearn.ensemble.gradient_boosting'] = _gb
    sys.modules['sklearn.ensemble.forest'] = _forest
    sys.modules['sklearn.tree.tree'] = _tree_classes
    
    # Mapeo de la función de pérdida renombrada/movida
    if not hasattr(_gb, 'LeastSquaresError'):
        _gb.LeastSquaresError = getattr(_gb, 'HalfSquaredError', getattr(_gb, 'LossFunction', None))
    
    print("[INFO] Monkey-patch de scikit-learn aplicado correctamente en sys.modules.")
except Exception as patch_err:
    print(f"[WARN] Inconveniente al aplicar monkey-patch de sklearn: {patch_err}")

# Cargar ruta de Azimuth local al PYTHONPATH
azimuth_path = os.path.expanduser("~/Azimuth")
if azimuth_path not in sys.path:
    sys.path.append(azimuth_path)


def generate_rs2_scores():
    input_csv = "results/horlbeck_for_rs2_export.csv"
    output_csv = "results/rs2_horlbeck_scored.csv"
    
    if not os.path.exists(input_csv):
        print(f"[ERROR] No existe el archivo de entrada: {input_csv}")
        sys.exit(1)
        
    print(f"[INFO] Cargando dataset exportado: {input_csv}...")
    df = pd.read_csv(input_csv)
    print(f"[INFO] Total de registros: {len(df)}")
    
    try:
        import azimuth.model_comparison as mc
        print("[INFO] Módulo Azimuth importado exitosamente.")
        
        # Invocación de la predicción Azimuth V3 usando la secuencia de 30 nucleótidos
        sequences = df['context_30nt'].values
        print("[INFO] Ejecutando mc.predict() sobre secuencias 30-mer...")
        predictions = mc.predict(sequences, None, None)
        
        df['rs2_score'] = predictions
        
        # Exportar únicamente las columnas requeridas
        out_df = df[['guide_sequence', 'rs2_score']].dropna()
        out_df.to_csv(output_csv, index=False)
        print(f"[ÉXITO] Matriz de scores guardada en: {output_csv} ({len(out_df)} guías).")
        
    except Exception as e:
        print(f"[ERROR CRÍTICO] Azimuth falló durante la inferencia: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    generate_rs2_scores()