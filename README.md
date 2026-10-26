# DL_entry_exclusion_gene

Components for exploratory analysis of entry-exclusion protein datasets.

## Modules

| Module | Purpose | Status |
|---|---|---|
| [gtdb_background](modules/gtdb_background/README.md) | Sample short proteins and their ±7 CDS genomic neighborhoods from GTDB-selected NCBI assemblies | CPU pipeline; Linux setup and Slurm template included |

The background sampler is an independent directory inside this repository, not a
Git submodule. Its records remain **unlabeled backgrounds**, not experimentally
confirmed negatives. ESM inference and classifier training are not included.

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
