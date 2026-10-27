"""Check exact-sequence reuse, provenance and pooling without downloading a model."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from embed import Cache, DIM, fetch_atlas, import_legacy, main, read_fasta, sha


class EmbeddingTests(unittest.TestCase):
    def test_reuse_by_sequence_hash_and_export_original_ids(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);seq='MAKDEFGHIK'
            (root/'input.faa').write_text('>new_id\n'+seq+'\n>alias_id\n'+seq+'\n')
            (root/'panel.json').write_text(json.dumps([dict(id='old_id',sequence=seq,sequence_sha256=sha(seq))]))
            np.savez(root/'old.npz',ids=np.array(['old_id']),embeddings=np.ones((1,DIM),dtype=np.float32))
            (root/'manifest.json').write_text(json.dumps(dict(model='ESM2 t36 3B UR50D',pooling='Final layer mean; full sequences; no truncation')))
            main(['--stage','prepare','--fasta',str(root/'input.faa'),'--out',str(root/'out'),
                  '--reuse',str(root/'panel.json'),str(root/'old.npz'),str(root/'manifest.json')])
            main(['--stage','export','--fasta',str(root/'input.faa'),'--out',str(root/'out')])
            with np.load(root/'out/embeddings.npz',allow_pickle=False) as data:
                self.assertEqual(list(data['ids']),['new_id','alias_id'])
                self.assertEqual(data['embeddings'].shape,(2,DIM))
            report=json.loads((root/'out/run_prepare.json').read_text())
            self.assertEqual(report['imported_legacy'],1);self.assertEqual(report['remaining_unique'],0)

    def test_incompatible_legacy_model_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);cache=Cache(root/'cache.sqlite')
            try:
                (root/'manifest.json').write_text(json.dumps(dict(model='ESM C',pooling='mean')))
                with self.assertRaisesRegex(ValueError,'not compatible'):
                    import_legacy(cache,root/'unused',root/'unused',root/'manifest.json',set())
            finally:cache.close()

    def test_atlas_similar_sequence_cannot_replace_query(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);cache=Cache(root/'cache.sqlite');digest=sha('MAKDEFGHIK')
            try:
                (root/'map.tsv').write_text('sequence_sha256\tatlas_id\n'+digest+'\tMGYP123\n')
                with patch('embed.fetch_text',return_value='>MGYP123\nMAKDEFGHIL\n') as request,patch('embed.time.sleep'):
                    audit=fetch_atlas(cache,root/'map.tsv',{digest},root)
                self.assertEqual(request.call_count,1)
                self.assertEqual(audit[0]['status'],'not_reused');self.assertIsNone(cache.get(digest))
            finally:cache.close()

    def test_atlas_exact_sequence_is_downloaded_with_provenance(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);cache=Cache(root/'cache.sqlite');digest=sha('MAKDEFGHIK')
            try:
                (root/'map.tsv').write_text('sequence_sha256\tatlas_id\n'+digest+'\tMGYP123\n')
                responses=['>MGYP123\nMAKDEFGHIK\n',json.dumps(dict(id='MGYP123',embedding=[1.0]*DIM))]
                with patch('embed.fetch_text',side_effect=responses),patch('embed.time.sleep'):
                    fetch_atlas(cache,root/'map.tsv',{digest},root)
                value=cache.get(digest);self.assertEqual(value[1],'atlas');self.assertTrue(value[2]['sequence_verified'])
            finally:cache.close()

    def test_invalid_sequences_not_silently_trimmed(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'input.faa';path.write_text('>x\nMAX*\n')
            with self.assertRaises(ValueError):read_fasta(path)


if __name__=='__main__':unittest.main()
