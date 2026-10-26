"""Verify stage isolation, argument forwarding and failures without network access."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from run_pipeline import main


class PipelineTests(unittest.TestCase):
    def setUp(self):
        platform_patch = patch('run_pipeline.platform.platform', return_value='test-platform')
        platform_patch.start()
        self.addCleanup(platform_patch.stop)

    def test_offline_analysis_and_paths_with_spaces(self):
        with tempfile.TemporaryDirectory(prefix='background test ') as folder:
            out = Path(folder)
            gbff = out / 'input with spaces.gbff'
            gbff.touch()
            (out / 'download_manifest.json').write_text(json.dumps([
                {'status': 'downloaded', 'gbff': str(gbff)}]), encoding='utf8')
            with patch('run_pipeline.subprocess.run') as run:
                main(['--stage', 'analyze', '--out', str(out), '--cache', str(out / 'cache'), '--flank', '3'])
            calls = [c.args[0] for c in run.call_args_list]
            self.assertEqual([Path(c[3]).name for c in calls],
                             ['test_extraction.py', 'extract_windows.py', 'build_report.py'])
            self.assertEqual(calls[1][calls[1].index('--flank') + 1], '3')
            self.assertEqual(json.loads((out / 'run_analyze.json').read_text())['status'], 'complete')

    def test_failed_step_stops_pipeline(self):
        with tempfile.TemporaryDirectory() as folder:
            with patch('run_pipeline.subprocess.run', side_effect=subprocess.CalledProcessError(1, 'download')) as run:
                with self.assertRaises(subprocess.CalledProcessError):
                    main(['--stage', 'prepare', '--out', folder, '--cache', folder])
            self.assertEqual(run.call_count, 1)
            self.assertEqual(json.loads((Path(folder) / 'run_prepare.json').read_text())['status'], 'failed')

    def test_incomplete_downloads_rejected_before_extraction(self):
        with tempfile.TemporaryDirectory() as folder:
            (Path(folder) / 'download_manifest.json').write_text('[{"status":"download_failed"}]')
            with patch('run_pipeline.subprocess.run') as run:
                with self.assertRaisesRegex(ValueError, 'Incomplete downloads'):
                    main(['--stage', 'analyze', '--out', folder, '--cache', folder])
            self.assertEqual(run.call_count, 1)


if __name__ == '__main__':
    unittest.main()
