# DL_entry_exclusion_gene

Components for exploratory analysis of entry-exclusion protein datasets.

## Modules

| Module | Purpose | Status |
|---|---|---|
| [gtdb_background](modules/gtdb_background/README.md) | Sample short proteins and their ±7 CDS genomic neighborhoods from GTDB-selected NCBI assemblies | CPU pipeline; Linux setup and Slurm template included |
| [protein_embeddings](modules/protein_embeddings/README.md) | Reuse sequence-matched vectors, verify/download Atlas entries, infer missing ESM2-3B representations | Resumable GPU pipeline and Slurm template |

The background sampler is an independent directory inside this repository, not a
Git submodule. Its records remain **unlabeled backgrounds**, not experimentally
confirmed negatives. Classifier training is not included.

Existing GTDB R214 installations can use `--source local-r214` to read local
metadata and Prodigal protein FASTAs without downloading anything. The original
`ncbi-r232` mode remains the default. From the module directory, with the Python
environment activated:

```bash
python run_pipeline.py --source local-r214 \
  --metadata /path/to/bac120_metadata_r214.tsv \
  --faa-root /path/to/protein_faa_reps \
  --out results/r214_pilot --cache cache/r214
```

For Slurm, pass the same source arguments after `slurm_analyze.sbatch` and set
`GTDB_OUT`/`GTDB_CACHE` to the preparation paths. A conda Python 3.11 environment
is also supported: set `GTDB_PYTHON` to that environment's `bin/python`.

## Linux quick start

```bash
git clone https://github.com/KuangHu/DL_entry_exclusion_gene.git
cd DL_entry_exclusion_gene/modules/gtdb_background
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s . -p 'test_*.py' -v
python run_pipeline.py --stage prepare --out results/pilot --cache cache/gtdb_r232
sbatch slurm_analyze.sbatch
```

Without Slurm, replace the last command with:

```bash
python run_pipeline.py --stage analyze --out results/pilot --cache cache/gtdb_r232
```

The preparation node needs internet access. Analysis runs offline using shared
cache paths. See the module README for resource settings, optional exact-positive
exclusion and custom sampling. Data, local environments and model weights are
excluded from Git.
