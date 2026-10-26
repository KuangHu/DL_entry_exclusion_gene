#!/usr/bin/env bash
#SBATCH --job-name=gtdb-background
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=gtdb-background-%j.log

# Starting resource requests only; adjust for your sampling size and cluster policy.
# Run from modules/gtdb_background after the prepare stage has completed.
# Add your site's account/partition via sbatch command-line flags if required.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit with sbatch from the module directory}"
PYTHON="${GTDB_PYTHON:-$PWD/.venv/bin/python}"
OUT="${GTDB_OUT:-$PWD/results/pilot}"
CACHE="${GTDB_CACHE:-$PWD/cache/gtdb_r232}"
export MPLBACKEND=Agg
export PYTHONUTF8=1
export MPLCONFIGDIR="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/gtdb-matplotlib-${SLURM_JOB_ID}"
mkdir -p "$MPLCONFIGDIR"
# Extra arguments after the script name are forwarded, including matching extraction
# settings and --known-positive-panel. This stage never downloads data.
"$PYTHON" run_pipeline.py --stage analyze --out "$OUT" --cache "$CACHE" "$@"
