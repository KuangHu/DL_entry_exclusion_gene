"""Portable entry point for GTDB background sampling (Python 3.11+)."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

HERE = Path(__file__).resolve().parent
METADATA_URL = 'https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_metadata_r232.tsv.gz'
METADATA_MD5 = '0ba4237077b65cfc5556e1d8be797485'


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--cache', type=Path, required=True)
    p.add_argument('--stage', choices=('all', 'prepare', 'analyze', 'report'), default='all',
                   help='prepare downloads; analyze extracts and reports without network')
    p.add_argument('--genomes', type=int, default=20)
    p.add_argument('--per-genome', type=int, default=20)
    p.add_argument('--min-aa', type=int, default=40)
    p.add_argument('--max-aa', type=int, default=300)
    p.add_argument('--flank', type=int, default=7)
    p.add_argument('--seed', default='gtdb-background-v1')
    p.add_argument('--known-positive-panel', type=Path,
                   help='Optional existing JSON panel; omission disables exact-positive filtering')
    return p


def main(argv=None):
    p = parser()
    a = p.parse_args(argv)
    if a.genomes < 1 or a.per_genome < 1 or a.min_aa < 1 or a.max_aa < a.min_aa or a.flank < 0:
        p.error('Counts and lengths must be positive, max-aa >= min-aa, flank >= 0')
    if sys.version_info < (3, 11):
        p.error('Python 3.11 or newer is required')
    if sys.flags.optimize:
        p.error('Do not use python -O or PYTHONOPTIMIZE: source validation uses assertions')
    if a.known_positive_panel and not a.known_positive_panel.is_file():
        p.error('Known-positive panel does not exist')
    a.out = a.out.resolve()
    a.cache = a.cache.resolve()
    a.out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PYTHONUTF8='1', MPLBACKEND='Agg')

    def run(script, *args):
        print('Running ' + script, flush=True)
        subprocess.run([sys.executable, '-X', 'utf8', str(HERE / script), *map(str, args)],
                       check=True, env=env)

    versions = {name: importlib.metadata.version(name) for name in ('biopython', 'numpy', 'matplotlib')}
    invocation = dict(python=sys.version, platform=platform.platform(), dependencies=versions,
                      arguments={k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()},
                      status='started')
    record_path = a.out / ('run_' + a.stage + '.json')

    def save():
        record_path.write_text(json.dumps(invocation, indent=2), encoding='utf8')

    save()
    try:
        if a.stage in ('all', 'analyze'):
            run('test_extraction.py')
        if a.stage in ('all', 'prepare'):
            metadata = a.cache / 'bac120_metadata_r232.tsv.gz'
            run('download.py', METADATA_URL, metadata, '--md5', METADATA_MD5, '--max-mib', 350)
            run('select_genomes.py', '--metadata', metadata, '--out', a.out,
                '--genomes', a.genomes, '--seed', a.seed)
            run('fetch_assemblies.py', '--selection', a.out / 'selected_assemblies.json',
                '--cache', a.cache / 'assemblies', '--out', a.out / 'download_manifest.json')
        if a.stage in ('all', 'analyze'):
            manifest = a.out / 'download_manifest.json'
            rows = json.loads(manifest.read_text(encoding='utf8'))
            if not rows or any(r.get('status') != 'downloaded' for r in rows):
                raise ValueError('Incomplete downloads: rerun prepare before analysis')
            for row in rows:
                if not Path(row['gbff']).is_file():
                    raise FileNotFoundError('GBFF path unavailable; rerun prepare with the current cache location: ' + row['gbff'])
            if not a.known_positive_panel:
                print('NOTICE: no positive panel supplied; exact-positive filtering is disabled. All cases remain unlabeled.', flush=True)
            optional = ['--known-positive-panel', a.known_positive_panel.resolve()] if a.known_positive_panel else []
            run('extract_windows.py', '--downloads', manifest, '--out', a.out,
                '--min-aa', a.min_aa, '--max-aa', a.max_aa, '--per-genome', a.per_genome,
                '--flank', a.flank, '--seed', a.seed, *optional)
        if a.stage in ('all', 'analyze', 'report'):
            run('build_report.py', '--out', a.out)
        invocation['status'] = 'complete'
    except Exception as exc:
        invocation.update(status='failed', error=str(exc))
        raise
    finally:
        save()
    print(a.out / ('download_manifest.json' if a.stage == 'prepare' else 'index.html'))


if __name__ == '__main__':
    main()
