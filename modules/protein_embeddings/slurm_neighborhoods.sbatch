#!/usr/bin/env bash
#SBATCH --job-name=esm2-neighbors
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=esm2-neighbors-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from modules/protein_embeddings}"
: "${ESM_PYTHON:?Set GPU-environment Python}"
: "${ESM_OUT:?Set extracted handoff directory}"
: "${ESM_MODEL_DIR:?Set local ESM2-3B or ESMFold-v1 checkpoint}"
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONUTF8=1
"$ESM_PYTHON" -u embed_neighborhoods.py --out "$ESM_OUT" \
  --model-dir "$ESM_MODEL_DIR" --model-format "${ESM_MODEL_FORMAT:-esmfold-stem}" \
  --batch-size "${ESM_BATCH_SIZE:-16}" --token-budget "${ESM_TOKEN_BUDGET:-4096}" "$@"
