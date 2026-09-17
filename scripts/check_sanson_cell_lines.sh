#!/bin/bash
# scripts/check_sanson_cell_lines.sh
set -euo pipefail

echo "=== Archivos en data/raw relacionados con Sanson ==="
find data/raw -iname "*sanson*" -o -iname "*dolcetto*" -o -iname "*41467*" 2>/dev/null | head -40

echo ""
echo "=== Buscar Jurkat / MCF7 / HT29 / A375 en nombres y contenido (muestra) ==="
find data -iname "*sanson*" -o -iname "*dolcetto*" 2>/dev/null | while read f; do
  echo "--- $f ---"
  case "$f" in
    *.xlsx|*.xls)
      python - << PY
import pandas as pd
from pathlib import Path
p = Path("$f")
try:
    xl = pd.ExcelFile(p)
    print("  sheets:", xl.sheet_names)
    for sh in xl.sheet_names[:8]:
        df = pd.read_excel(p, sheet_name=sh, nrows=3)
        cols = [str(c).lower() for c in df.columns]
        print(f"  sheet={sh!r} cols_sample={list(df.columns)[:12]}")
        joined = " ".join(cols)
        for term in ["jurkat", "mcf7", "mcf-7", "ht29", "ht-29", "a375", "k562"]:
            if term in joined:
                print(f"    FOUND term in columns: {term}")
except Exception as e:
    print("  excel error:", e)
PY
      ;;
    *.csv|*.tsv|*.txt|*.md)
      grep -i -n -E "jurkat|mcf7|mcf-7|ht29|a375|k562" "$f" 2>/dev/null | head -5 || true
      ;;
  esac
done

echo ""
echo "=== Columnas de sanson2018_lfc_scores.csv ==="
python - << 'EOF'
import pandas as pd
df = pd.read_csv("data/processed/sanson2018_lfc_scores.csv", nrows=2)
print(df.columns.tolist())
print(df.head(1).T)
EOF

echo ""
echo "=== Referencias en paper / logs ==="
grep -r -i -n -E "jurkat|mcf7|ht29|a375" logs/ data/raw/crispr_datasets/Sanson2018/ 2>/dev/null | head -30 || true