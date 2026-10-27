#!/usr/bin/env bash
#SBATCH --job-name=esm2-centers
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=esm2-centers-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from modules/protein_embeddings}"
: "${ESM_PYTHON:?Set to the GPU environment Python}"
: "${ESM_FASTA:?Set to center_proteins.faa}"
: "${ESM_OUT:?Set to the prepared embedding output directory}"
: "${ESM_MODEL_DIR:?Set to a locally downloaded ESM2-3B checkpoint}"
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export PYTHONUTF8=1
"$ESM_PYTHON" embed.py --stage infer --fasta "$ESM_FASTA" --out "$ESM_OUT" \
  --model-dir "$ESM_MODEL_DIR" --model-format "${ESM_MODEL_FORMAT:-esm2}" "$@"
