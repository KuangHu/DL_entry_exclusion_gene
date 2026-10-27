import hashlib,json,sqlite3,tempfile,unittest
from pathlib import Path
import numpy as np
from embed import Cache,SIGNATURE,sha
from neighborhood_io import output_lock,unique_rows
from verify_handoff import file_sha,verify

class NeighborhoodIOTests(unittest.TestCase):
    def test_duplicate_full_sequences_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'in.faa';p.write_text('>a\nMAKL\n>b\nMAKL\n')
            with self.assertRaises(ValueError):unique_rows(p)

    def test_writer_exclusion_and_release(self):
        with tempfile.TemporaryDirectory() as d:
            with output_lock(d):
                with self.assertRaises(RuntimeError):
                    with output_lock(d):pass
            with output_lock(d):pass

    def test_resume_archive_integrity_and_corruption(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'all_proteins.faa').write_text('>a\nMAKL\n>b\nMMKL\n')
            c=Cache(p/'embedding_cache.sqlite');c.put(sha('MAKL'),np.ones(2560),'prior_local',{});c.commit();c.close()
            files={name:{'sha256':file_sha(p/name)} for name in ['all_proteins.faa','embedding_cache.sqlite']}
            (p/'handoff_manifest.json').write_text(json.dumps(dict(files=files,unique_sequences=2,saved_vectors=1)))
            self.assertEqual(verify(p)['remaining'],1)
            (p/'all_proteins.faa').write_text('>a\nMAKQ\n>b\nMMKL\n')
            with self.assertRaisesRegex(ValueError,'Checksum mismatch'):verify(p)

if __name__=='__main__':unittest.main()
