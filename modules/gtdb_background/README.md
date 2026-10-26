# GTDB short-protein background windows

This pipeline collects **unlabeled background cases**, not experimentally confirmed EEX negatives. It does not infer new EEX functions or change existing project annotations.

The completed pilot uses GTDB **R232 / 232.0** to select exact NCBI assembly accessions and assign taxonomy. Protein calls, coordinates and product descriptions come from the corresponding **NCBI RefSeq GBFF**, not from GTDB's Prodigal protein files. Raw downloads, assembly reports and checksums are cached for reproducibility.

## Linux cluster setup

This is the independent `modules/gtdb_background` component of the main repository,
not a Git submodule. Sampling needs Python 3.11 and CPU only; no GPU, ESM weights,
SignalP, DeepTMHMM or API key is required.

```bash
git clone https://github.com/KuangHu/DL_entry_exclusion_gene.git
cd DL_entry_exclusion_gene/modules/gtdb_background
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s . -p 'test_*.py' -v
python run_pipeline.py --out results/pilot --cache cache/gtdb_r232
```

The final atlas is `results/pilot/index.html`; copy the entire output directory
when viewing it elsewhere, including its JavaScript and image files. Defaults
are 20 genera, 20 centers per genome, 40–300 aa, ±7 CDS. A larger run can use
`--genomes 100 --per-genome 20 --min-aa 40 --max-aa 300 --flank 7`.
Use a different `--out` directory for each run.

Pinned requirements match the original pilot environment. Each invocation writes
arguments, platform and dependency versions to `run_<stage>.json`. Do not use
`python -O` or `PYTHONOPTIMIZE`, since source validation uses assertions.

### Prepare online, analyze offline

Run preparation on a node where your cluster permits outbound HTTPS to GTDB and
NCBI. Then submit from this module directory:

```bash
python run_pipeline.py --stage prepare --out results/pilot --cache cache/gtdb_r232
sbatch slurm_analyze.sbatch
```

The Slurm template requests 2 CPUs, 8 GiB and 4 hours as starting settings,
not measured minimums or runtime guarantees. Adjust for scale and cluster policy.
Add your site's account/partition using `sbatch --account=... --partition=...`.
The analysis stage never downloads data. Keep the Python environment and cache on
paths accessible from compute nodes. Without Slurm, run directly:

```bash
python run_pipeline.py --stage analyze --out results/pilot --cache cache/gtdb_r232
```

Custom shared paths and sampling settings:

```bash
export GTDB_OUT=/shared/your_project/background_run1
export GTDB_CACHE=/shared/your_project/gtdb_r232
python run_pipeline.py --stage prepare --out "$GTDB_OUT" --cache "$GTDB_CACHE" --genomes 100
sbatch slurm_analyze.sbatch --per-genome 20 --min-aa 40 --max-aa 300 --flank 7
```

`GTDB_PYTHON` overrides the job's default `.venv/bin/python`. Genome count affects
preparation; per-genome count, length bounds and flank affect analysis. Pass the
same custom `--seed` to both stages. Do not run concurrent preparation jobs using
the same cache, or concurrent jobs with the same output directory. Downloads
resume at file level; interrupted files restart. Absolute GBFF paths are stored
in the manifest: after relocating a cache, rerun `prepare` with the new cache
path to regenerate the manifest from cached ZIP packages.

### Optional exact-positive exclusion

Only code is included in Git. The earlier project's positive panel is not bundled.
Without it, the runner prints a notice and records zero known positives checked;
source annotation review still applies. This differs from the original local
pilot, which checked 1,047 exact positive sequences. Supply the panel to `all` or
`analyze`, or after the Slurm script name:

```bash
python run_pipeline.py --stage analyze --out results/pilot --cache cache/gtdb_r232 --known-positive-panel /shared/your_project/panel.json
```

The panel is a JSON array; records with `"comparison_group": "EEX_existing_label"`
need `"sequence_sha256"`, the lowercase SHA256 of the full amino-acid sequence.
Other fields are ignored; raw positive sequences are unnecessary for this exact
match filter. No homology screening is performed. All cases remain unlabeled.

The portable Python entry point also works on Windows. Individual tools expose `--help`.

## Selection and extraction

1. Select GTDB bacterial species representatives with RefSeq accessions, Complete Genome assembly level, CheckM2 completeness ≥95 and contamination ≤5, and no MAG/SAG/environmental-assembly category.
2. Sample classes deterministically, with one assembly per genus. This is a diversity pilot, not a population-weighted or EEX-taxonomy-matched sample.
3. Download exact accession versions via NCBI Datasets. Verify the returned assembly accession, ZIP CRC and GBFF SHA256. Do not silently substitute a newer assembly.
4. Order **all** CDS by genomic start, end and protein accession. Pseudogenes, partial CDS and untranslated CDS retain their index.
5. Eligible centers have full standard amino-acid translations within the chosen length interval. Exact project EEX matches and source exclusion-name annotations go to the audit rather than the background sample. Name checks are conservative review flags, not functional evidence.
6. Select centers by deterministic SHA256 rank, round-robin over length intervals. Exact centers are deduplicated globally; selected windows do not share CDS within one source contig.
7. Extract by CDS index and orient positions to the center strand. Preserve neighbor strands; never cross contigs. Wrap around an origin only on source-declared circular records with enough CDS to avoid repeating neighbors.
8. Keep unavailable protein translations in their original CDS slots. They are not removed before counting ±7.

The standard length strata are bounded by 80, 120 and 200 aa, clipped to the configured range. These strata are sampling choices, not estimates of biological abundance.

## Outputs

| File | Meaning |
|---|---|
| index.html | Main report and searchable neighborhood diagrams |
| centers.tsv | Center ID, sequence hash, accession, taxonomy, source evidence and coordinates |
| center_proteins.faa | Center proteins keyed by case_id (BG...) |
| neighbors.tsv | Every CDS position, including excluded translations and original gene/product fields |
| window_proteins_unique.faa | Exact-deduplicated valid proteins keyed by sequence_id (BP...) |
| windows.jsonl | Complete mapping of centers and neighboring positions |
| architectures.json | Architecture identity, copies, center sequence counts and frequencies |
| architecture_position_diversity.tsv | Distinct aa/valid-copy and missing/excluded counts at every position |
| center_sampling_audit.tsv | Short centers and reasons for sampling, exclusion or nonselection |
| selected_assemblies.tsv | Exact assembly selection and taxonomy |
| download_manifest.json | Source locations, exact accessions and checksums |
| validation.json | Verification against original GBFF coordinates and translations |

**negative_label remains null** and **experimental_negative remains false** for every case. Function-annotated and uncharacterized centers remain separate evidence groups. Conjugation/exclusion annotation text in a window is retained as a review flag; its absence does not establish biological absence.

A copy is one sampled center occurrence, not an independent biological replicate. Exact deduplication does not remove near homologs. The genus split_group is provided for later grouping; it is **not** a completed train/test split.

## Validation and next use

The edge-case tests cover reverse-strand indexing, linear edges, circular wrapping, tiny circular contigs, invalid CDS positions and partial annotations. Report generation independently checks every window against the original GBFF, verifies valid FASTA sequences, and checks nonoverlap and architecture counts.

Before comparing EEX embeddings, match taxonomy, length, signal peptides, topology and window visibility. This pilot's full windows differ from the truncated windows in the earlier EEX dataset. Gene-calling and annotation differences are also potential confounders. Near-homolog grouping, phenotype-negative validation and center ESM/SignalP/topology calculations are outside this extraction stage.

## Sources

- [Pinned GTDB R232 directory](https://data.gtdb.ecogenomic.org/releases/release232/232.0/)
- [GTDB genome-download guidance](https://gtdb.ecogenomic.org/faq)
- [NCBI genome data package](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/)
- [NCBI Datasets API specification](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/openapi3/openapi3.docs.yaml)
