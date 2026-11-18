#!/usr/bin/env bash
set -euo pipefail
TASK_ROOT="$(cd "$(dirname "$0")" && pwd)"
if [[ ! -x "$TASK_ROOT/.env/bin/python" ]]; then
  export CONDA_PKGS_DIRS="$TASK_ROOT/cache/conda-pkgs"
  conda env create --prefix "$TASK_ROOT/.env" --file "$TASK_ROOT/config/environment.yml"
fi
export PATH="$TASK_ROOT/.env/bin:$PATH"
"$TASK_ROOT/.env/bin/python" "$TASK_ROOT/scripts/install_prodigal.py"
"$TASK_ROOT/.env/bin/python" "$TASK_ROOT/scripts/fetch_public.py"
"$TASK_ROOT/.env/bin/python" "$TASK_ROOT/scripts/build_reference.py" --workbook "${1:?Pass the function reference XLSX path}"
"$TASK_ROOT/.env/bin/python" -m unittest discover -s "$TASK_ROOT/tests" -v
"$TASK_ROOT/.env/bin/python" "$TASK_ROOT/scripts/validate_workflow.py"
