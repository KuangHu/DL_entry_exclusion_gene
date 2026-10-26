"""Regression coverage for source selection and an offline R214 report round trip."""
import contextlib
import csv
import hashlib
import io
import json
from pathlib import Path
import random
import tempfile
import unittest

from run_pipeline import main
from select_genomes import select


def metadata(path, source):
    row = dict(accession='RS_GCF_000000001.1', gtdb_representative='t',
               gtdb_taxonomy='d__Bacteria;p__Test;c__Test;o__Test;f__Test;g__Test;s__Test',
               ncbi_assembly_level='Complete Genome', ncbi_genome_category='none')
    prefix = 'checkm' if source == 'local-r214' else 'checkm2'
    row.update({prefix + '_completeness': '99', prefix + '_contamination': '0'})
    with path.open('w', encoding='utf8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t')
        writer.writeheader()
        writer.writerow(row)


class SourceModeTests(unittest.TestCase):
    def test_ncbi_sampling_rank_matches_original_pipeline(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            metadata(root / 'metadata.tsv', 'ncbi-r232')
            with contextlib.redirect_stdout(io.StringIO()):
                select(root / 'metadata.tsv', root / 'out', count=1)
            row = json.loads((root / 'out/selected_assemblies.json').read_text())[0]
            # Before local-r214 was introduced, the seed used the bare assembly ID.
            expected = hashlib.sha256(b'gtdb-background-v1|GCF_000000001.1').hexdigest()
            self.assertEqual(row['sampling_rank'], expected)

    def test_local_pipeline_generates_source_accurate_report(self):
        with tempfile.TemporaryDirectory(prefix='r214 round trip ') as folder:
            root = Path(folder)
            metadata(root / 'metadata.tsv', 'local-r214')
            proteins = root / 'proteins/bacteria'
            proteins.mkdir(parents=True)
            rng = random.Random(7)
            with (proteins / 'RS_GCF_000000001.1_protein.faa').open('w', encoding='utf8') as handle:
                for i in range(1, 41):
                    seq = 'M' + ''.join(rng.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(79))
                    start = 1 + (i - 1) * 300
                    strand = -1 if i % 2 else 1
                    handle.write(f'>ctg_{i} # {start} # {start+242} # {strand} # ID=1_{i};partial=00\n{seq}*\n')
            out = root / 'out'
            main(['--source', 'local-r214', '--metadata', str(root / 'metadata.tsv'),
                  '--faa-root', str(root / 'proteins'), '--cache', str(root / 'cache'),
                  '--out', str(out), '--genomes', '1', '--per-genome', '2'])
            validation = json.loads((out / 'validation.json').read_text())
            self.assertEqual(validation['source_CDS_positions_verified'], 30)
            self.assertEqual(validation['source_annotation_review'], 'unavailable')
            self.assertEqual(validation['cases_wrapping_circular_origin'], 0)
            page = (out / 'index.html').read_text(encoding='utf8')
            self.assertNotIn('@@', page)
            self.assertIn('未执行已知 EEX 精确匹配排除', page)
            self.assertIn('源功能注释审查不可用', page)
            self.assertIn('没有 GTDB 官方 MD5 校验', page)
            self.assertNotIn('<option value="false">未见标记', page)
            self.assertNotIn('0 个抽样窗口含接合', page)


if __name__ == '__main__':
    unittest.main()
