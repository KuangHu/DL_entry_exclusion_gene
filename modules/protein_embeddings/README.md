# Protein embeddings: reuse, exact Atlas download, frozen ESM2 inference

This module produces sequence representations; it does not assign protein functions
or negative labels. The representation is ESM2 t36 3B UR50D, final layer 36,
2560 dimensions, averaged over the full residue sequence with BOS/EOS and padding
excluded. Inference uses fp16 on CUDA and fp32 pooling/output. Existing ESM C or
other ESM2 model sizes cannot be mixed into this vector space.

## Environment

Use a dedicated Python 3.11 GPU environment. Install a CUDA-enabled PyTorch >=2.1
compatible with your cluster, then `python -m pip install -r requirements.txt`.
The sampling-only environment does not need to be modified. Download/stage the
official `facebook/esm2_t36_3B_UR50D` Hugging Face checkpoint directory on a network
node before running inference. All model loading uses `local_files_only=True`.

An existing `facebook/esmfold_v1` checkpoint can supply its frozen ESM2 stem using
`--model-format esmfold-stem`. This is the format used by the earlier project
positive vectors. Only the language model is evaluated.

## 1. Reuse existing vectors

From this module directory:

```bash
python embed.py --stage prepare --fasta /path/to/center_proteins.faa --out /path/to/esm_centers \
  --reuse /path/to/panel.json /path/to/embeddings.npz /path/to/embedding_manifest.json
```

`--reuse` can be repeated. The panel must map its `id` to full-sequence
`sequence_sha256`; NPZ `ids` and `embeddings` must agree with the panel and the
manifest must specify the compatible model/full-sequence mean. Matches use the
sequence SHA256, not the FASTA identifier. The earlier v2 panel is directly
supported. Without a compatible cache, run prepare without `--reuse`.

## 2. Optional Atlas downloads

ESM Atlas provides individual precomputed embeddings under MGYP IDs. It does not
provide a direct GTDB-accession-to-embedding mapping. The project sequence-search
probe returned HTTP 503 on 2026-10-26; this does not establish sequence absence.

If a mapping is available from an independently obtained exact-sequence lookup,
provide a TSV with columns `sequence_sha256` and `atlas_id`:

```bash
python embed.py --stage prepare --fasta /path/to/center_proteins.faa --out /path/to/esm_centers \
  --atlas-map /path/to/atlas_map.tsv
```

Each candidate Atlas sequence is fetched and its full SHA256 must exactly match
before downloading its vector. Similar sequences are rejected. Failed downloads
are logged and remain pending; Atlas metadata coverage is not inferred from API
failure. This stage needs network access. Without an Atlas index/map, skip this
optional stage and compute pending sequences directly.

## 3. Infer pending proteins offline

```bash
python embed.py --stage infer --fasta /path/to/center_proteins.faa --out /path/to/esm_centers \
  --model-dir /path/to/esm2_t36_3B_UR50D --batch-size 16 --token-budget 4096
```

Or submit from this module directory:

```bash
export ESM_PYTHON=/path/to/gpu_env/bin/python
export ESM_FASTA=/path/to/center_proteins.faa
export ESM_OUT=/path/to/esm_centers
export ESM_MODEL_DIR=/path/to/esm2_t36_3B_UR50D
sbatch --account=... --partition=... slurm_embed.sbatch
```

For an ESMFold checkpoint, set `ESM_MODEL_FORMAT=esmfold-stem`.
The Slurm file requests one GPU, 8 CPUs and 32 GiB host RAM; adjust site settings
and resources. For long proteins a separate policy is needed: the default maximum
is 1022 residues and overlength inputs fail explicitly, with no silent truncation.
The current 40–300 aa centers fit this policy. A CUDA out-of-memory error reduces
batch size automatically; a single-sequence failure stops the job.

Vectors are committed to `embedding_cache.sqlite` after every batch. Rerun the
same command to resume; only missing sequence hashes are inferred. Use one writer
per output directory. `--limit 128` is suitable for a smoke test; rerun without
the limit for the full input. The signature, full input checksum, scope limit and
per-sequence provenance are recorded, so a limited run is not a full completion.

## Outputs and verification

- `embeddings.npz`: input-order `ids`, full `sequence_sha256`, float32 `embeddings`.
- `sequence_manifest.tsv`: ID, sequence hash, length and source for every input row.
- `embedding_cache.sqlite`: resumable vectors and detailed source records.
- `pending.faa`: unique sequences still missing after the last successful stage.
- `run_<stage>.json`: configuration, status, source counts and input checksum.
- `atlas_download_audit.json`: optional download outcomes including rejected matches.

Before downstream comparisons, check input order/hashes, finite vectors and model
signature, and compare a small computed set against existing cached vectors. Atlas
and GPU results can have small floating point differences. This module computes
no classifier accuracy or protein-function inference.

```bash
python verify.py --fasta /path/to/center_proteins.faa --out /path/to/esm_centers
```

`verify.py` checks the full input order, hashes, vector integrity and per-sequence
cache/provenance, and writes `embedding_validation.json`. For numerical checkpoint
calibration against a prior panel, use `validate_model.py --model-dir ...
--model-format esmfold-stem --panel .../panel.json --vectors .../embeddings.npz
--out .../model_calibration.json`. Optional `--atlas-fasta` and `--atlas-vector`
accept a previously downloaded, matching Atlas example.

Official Atlas download/API information:
https://github.com/facebookresearch/esm/blob/main/scripts/atlas/README.md

## Resume the complete neighborhood job on a cluster

Use `embed_neighborhoods.py` for the prepared **center + +/-7 CDS** union.
It runs full-length ESM2-3B with memory-efficient SDPA attention and validates
that implementation against the original attention on the allocated GPU.
Transferred cache vectors are also checked against the loaded checkpoint on a
small calibration panel. No local temperature throttling is applied by this
cluster entry point. Do not restart the Windows `run_full_length.py` job.

The local handoff contains 320,078 unique sequences, of which 92,979 vectors
were saved before stopping; 227,099 remain. These counts describe that frozen
handoff, not subsequent cluster progress. It includes 5,211 proteins >1022 aa;
they retain full context and an `extended_context` flag. Long-context embeddings
need a separate sensitivity analysis; numerical checks are not biological validation.

1. Transfer and extract `esm_neighborhood_cluster_handoff.tar.gz` into a new
   output directory on shared cluster storage. It contains `all_proteins.faa`,
   `embedding_cache.sqlite`, `slots.tsv`, `sequences.tsv`, and checksums.
   Do not run `prepare_inputs.py` again or overwrite this cache.
2. Update this repository and activate a dedicated **CUDA-enabled** Python 3.11
   environment with the dependencies in `requirements.txt`. This runner needs
   PyTorch >=2.4 (the Windows implementation was tested with 2.10.0) and a GPU
   supporting PyTorch Flash/Efficient SDPA. The sampling environment alone does
   not establish GPU readiness. A GPU with >=24 GB VRAM is recommended; batch
   OOMs are retried at a smaller batch size, single-sequence OOMs are reported.
3. Stage the same `facebook/esmfold_v1` checkpoint and set model format
   `esmfold-stem`, or supply the corresponding full `esm2_t36_3B_UR50D` checkpoint
   and model format `esm2`. Model downloads are not included in the archive.
   Compute jobs load only local files. The existing ESMFold checkpoint revision
   was `75a3841ee059df2bf4d56688166c8fb459ddd97a`.
4. Verify the **freshly extracted** handoff on a CPU node, before any cache writes:

   ```bash
   cd /groups/rubin/projects/kuang/out/DL_entry_exclusion/DL_entry_exclusion_gene/modules/protein_embeddings
   python verify_handoff.py /groups/rubin/projects/kuang/out/DL_entry_exclusion/results/esm_neighborhoods
   ```

5. Set `ESM_PYTHON` to the actual GPU environment Python, `ESM_MODEL_DIR` to the
   staged checkpoint directory, and `ESM_OUT` to the extracted output directory.
   Submit from this module directory with the site's actual account/partition:

   ```bash
   export ESM_OUT=/groups/rubin/projects/kuang/out/DL_entry_exclusion/results/esm_neighborhoods
   export ESM_MODEL_FORMAT=esmfold-stem
   # Also set ESM_PYTHON and ESM_MODEL_DIR to existing cluster paths.
   sbatch --account=<account> --partition=<gpu-partition> slurm_neighborhoods.sbatch
   ```

The job commits each completed batch, skips existing full-sequence hashes,
and holds an OS-owned single-writer lock. Retry the same submission after a
time limit or interruption; do not rerun the original archive checksum check
after the cache has changed. Only one job may write a given output directory.
Adjust the Slurm resource request to site policy. No account, partition, or
GPU environment is assumed to exist merely because it appears in an example.

Final output is memory-mappable `embeddings.npy` (320078 x 2560 float32), with
`sequence_manifest.tsv` mapping row numbers to full sequence hashes and sources.
`slots.tsv` maps each CDS occurrence to its sequence hash and preserves absent or
excluded positions. `embedding_validation.json` confirms every exported vector
against the cache and verifies valid-slot coverage. `progress.json` states
whether the run is actually complete. The cluster bundle intentionally excludes
the old local live HTML dashboard; use progress/logs and regenerate a viewer from
the finished files if needed.
