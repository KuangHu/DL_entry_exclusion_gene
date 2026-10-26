"""R214 local-source edge cases: Prodigal headers, shard-prefix dedup / conflicts, wrap-off."""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from extract_windows import _parse_prodigal_gene, _iter_prodigal_records, decode_feature, neighborhood_indices
from fetch_assemblies import _parse_prodigal_file, _biological_sha256, _build_local_r214_index, resolve_local_r214


def sha(s): return hashlib.sha256(s.encode()).hexdigest()


HEADER_WITH_PREFIX = '>GB_GCA_000007325.1_protein__AE009951.2_1 # 3 # 77 # -1 # ID=1_1;partial=10;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;gc_cont=0.173\n'
HEADER_PLAIN = '>AE009951.2_1 # 3 # 77 # -1 # ID=1_1;partial=10;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;gc_cont=0.173\n'
SEQ_LINES = ['MKKFIILLFILVQGLIFSATKTLSD\n']


def write_prodigal(path, records):
    """records: list of (contig, ordinal, start_1, end_1, strand, partial, seq)."""
    with Path(path).open('w', encoding='utf8') as f:
        for contig, ordinal, start, end, strand, partial, seq in records:
            f.write(f'>GB_GCA_X.1_protein__{contig}_{ordinal} # {start} # {end} # {strand} # ID=1_{ordinal};partial={partial};start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;gc_cont=0.5\n')
            f.write(seq + '\n')


class ProdigalHeaderTests(unittest.TestCase):
    def test_prefix_and_plain_headers_decode_to_same_gene(self):
        c1, o1, g1 = _parse_prodigal_gene(HEADER_WITH_PREFIX, SEQ_LINES, set())
        c2, o2, g2 = _parse_prodigal_gene(HEADER_PLAIN, SEQ_LINES, set())
        self.assertEqual((c1, o1), (c2, o2))
        for k in ('start0', 'end0', 'strand', 'prodigal_partial', 'sequence_sha256', 'sequence_status'):
            self.assertEqual(g1[k], g2[k], k)

    def test_coord_conversion_is_1_based_to_0_based(self):
        _, _, g = _parse_prodigal_gene('>c_1 # 3 # 77 # -1 # ID=1_1;partial=00\n', ['MA\n'], set())
        self.assertEqual(g['start0'], 2)
        self.assertEqual(g['end0'], 77)
        self.assertEqual(g['strand'], -1)

    def test_partial_marks_partial_CDS(self):
        _, _, g = _parse_prodigal_gene('>c_1 # 1 # 30 # 1 # ID=1_1;partial=10\n', ['MA\n'], set())
        self.assertEqual(g['sequence_status'], 'partial_CDS')
        _, _, g2 = _parse_prodigal_gene('>c_2 # 50 # 80 # 1 # ID=1_2;partial=00\n', ['MA\n'], set())
        self.assertEqual(g2['sequence_status'], 'valid')

    def test_prodigal_and_decode_feature_share_required_gene_fields(self):
        _, _, g = _parse_prodigal_gene('>c_1 # 1 # 30 # 1 # ID=1_1;partial=00\n', ['MA\n'], set())
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        f = SeqFeature(FeatureLocation(0, 30, strand=1), type='CDS', qualifiers={'translation': ['MA']})
        n = decode_feature(f, 0, set())
        shared = {'start0', 'end0', 'strand', 'sequence', 'sequence_sha256', 'sequence_status',
                  'aa_length', 'source_unknown', 'exclusion_annotation_flag',
                  'conjugation_annotation_flag', 'exact_project_EEX', 'annotation_review',
                  'gene', 'product', 'locus_tag', 'protein_accession', 'location'}
        self.assertTrue(shared <= set(g), shared - set(g))
        self.assertTrue(shared <= set(n), shared - set(n))
        self.assertEqual(g['annotation_review'], 'unavailable')
        self.assertEqual(n['annotation_review'], 'available')
        self.assertTrue(g['source_unknown'])
        self.assertFalse(g['exclusion_annotation_flag'])
        self.assertFalse(g['conjugation_annotation_flag'])

    def test_unparseable_header_raises(self):
        with self.assertRaises(ValueError):
            _parse_prodigal_gene('>no_coords here\n', ['MA\n'], set())


class ProdigalIterTests(unittest.TestCase):
    def test_wrap_off_and_ordinal_consistency(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'g.faa'
            # 3 CDS on one contig, in-order ordinals, strand +1 then -1 then +1
            write_prodigal(p, [('ctg', 1, 1, 30, 1, '00', 'MASDEFGHIK'),
                               ('ctg', 2, 50, 110, -1, '00', 'MMMMMMMMMMMMMMMMMMMM'),
                               ('ctg', 3, 200, 260, 1, '00', 'MASDEFGHIKASDEFGHIKA')])
            records = list(_iter_prodigal_records(p, set()))
        self.assertEqual(len(records), 1)
        rec = records[0]
        self.assertEqual(rec['id'], 'ctg')
        self.assertFalse(rec['circular'])
        self.assertEqual(len(rec['genes']), 3)
        self.assertEqual([g['cds_index'] for g in rec['genes']], [0, 1, 2])
        self.assertTrue(rec['ordinal_consistency']['gap_free'])
        self.assertTrue(rec['ordinal_consistency']['monotonic'])

    def test_ordinal_gap_is_reported_but_not_fatal(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'g.faa'
            write_prodigal(p, [('ctg', 1, 1, 30, 1, '00', 'MASDEFGHIK'),
                               ('ctg', 5, 200, 260, 1, '00', 'MMMMMMMMMMMMMMMM')])
            records = list(_iter_prodigal_records(p, set()))
        self.assertFalse(records[0]['ordinal_consistency']['gap_free'])
        self.assertTrue(records[0]['ordinal_consistency']['monotonic'])

    def test_wrap_off_edge_window_marked_incomplete(self):
        res = neighborhood_indices(5, 0, 1, flank=7, circular=False)
        self.assertTrue(any(i is None for _, i, _ in res))


class ResolveLocalR214Tests(unittest.TestCase):
    def _write_shard_and_plain(self, folder, bio_identical):
        bact = folder / 'bacteria'
        bact.mkdir()
        shard = bact / '2-GB_GCA_000012345.1_protein.faa'
        plain = bact / 'GB_GCA_000012345.1_protein.faa'
        # shard version has long prefixed id; plain uses bare id
        shard.write_text('>GB_GCA_000012345.1_protein__ctg_1 # 1 # 30 # 1 # ID=1_1;partial=00\nMASDEFGHIK\n', encoding='utf8')
        if bio_identical:
            plain.write_text('>ctg_1 # 1 # 30 # 1 # ID=1_1;partial=00\nMASDEFGHIK\n', encoding='utf8')
        else:
            plain.write_text('>ctg_1 # 1 # 30 # 1 # ID=1_1;partial=00\nMMMMMMMMMMMMMMMM\n', encoding='utf8')
        return bact

    def test_shard_prefix_duplicate_same_bio_content_dedups(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)
            self._write_shard_and_plain(folder, bio_identical=True)
            cache = folder / 'cache'
            cache.mkdir()
            index = _build_local_r214_index(folder)
            self.assertEqual(len(index['GB_GCA_000012345.1']), 2)
            row = dict(assembly_accession='GCA_000012345.1', accession_prefix='GB_')
            result = resolve_local_r214(row, index, cache)
        self.assertEqual(result['status'], 'downloaded')
        self.assertTrue(result['protein_faa'].endswith('2-GB_GCA_000012345.1_protein.faa'))  # shard preferred as canonical
        self.assertEqual(len(result['duplicate_paths']), 1)
        self.assertIn('biological_sha256', result)
        self.assertEqual(result['record_count'], 1)

    def test_shard_prefix_duplicate_different_bio_content_reports_conflict(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)
            self._write_shard_and_plain(folder, bio_identical=False)
            cache = folder / 'cache'
            cache.mkdir()
            index = _build_local_r214_index(folder)
            row = dict(assembly_accession='GCA_000012345.1', accession_prefix='GB_')
            result = resolve_local_r214(row, index, cache)
        self.assertEqual(result['status'], 'content_conflict')
        self.assertEqual(len(result['conflicting_paths']), 2)
        self.assertEqual(len(result['conflicting_biological_sha256']), 2)
        self.assertNotIn('protein_faa', result)

    def test_missing_file_reports_download_failed(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)
            (folder / 'bacteria').mkdir()
            cache = folder / 'cache'
            cache.mkdir()
            index = _build_local_r214_index(folder)
            row = dict(assembly_accession='GCA_999999999.1', accession_prefix='GB_GCA_')
            result = resolve_local_r214(row, index, cache)
        self.assertEqual(result['status'], 'download_failed')
        self.assertIn('FileNotFoundError', result['error'])

    def test_biological_sha256_invariant_to_header_prefix(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)
            self._write_shard_and_plain(folder, bio_identical=True)
            shard_recs = _parse_prodigal_file(folder / 'bacteria' / '2-GB_GCA_000012345.1_protein.faa')
            plain_recs = _parse_prodigal_file(folder / 'bacteria' / 'GB_GCA_000012345.1_protein.faa')
        self.assertEqual(_biological_sha256(shard_recs), _biological_sha256(plain_recs))

    def test_index_accepts_RS_GCF_and_GB_GCA(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)
            bact = folder / 'bacteria'
            bact.mkdir()
            (bact / 'GB_GCA_000000001.1_protein.faa').write_text('>c_1 # 1 # 30 # 1 # ID=1_1;partial=00\nMA\n', encoding='utf8')
            (bact / 'RS_GCF_000000002.1_protein.faa').write_text('>c_1 # 1 # 30 # 1 # ID=1_1;partial=00\nMA\n', encoding='utf8')
            index = _build_local_r214_index(folder)
        self.assertIn('GB_GCA_000000001.1', index)
        self.assertIn('RS_GCF_000000002.1', index)


if __name__ == '__main__':
    unittest.main()
