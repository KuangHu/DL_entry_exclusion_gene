"""Portable input and single-writer safeguards for neighborhood inference."""
from contextlib import contextmanager
import os
from pathlib import Path
from embed import read_fasta

def unique_rows(path):
    rows=read_fasta(path)
    if len({r['sequence_sha256'] for r in rows})!=len(rows):
        raise ValueError('Neighborhood FASTA must be deduplicated by full sequence')
    return rows

@contextmanager
def output_lock(folder):
    """OS-owned locks are released even after process death; no stale PID locks."""
    path=Path(folder)/'.embedding_writer.lock'
    with path.open('a+b') as handle:
        if os.name=='nt':
            import msvcrt
            handle.seek(0)
            try:msvcrt.locking(handle.fileno(),msvcrt.LK_NBLCK,1)
            except OSError as exc:raise RuntimeError('Another writer holds this output directory') from exc
            try:yield
            finally:
                handle.seek(0);msvcrt.locking(handle.fileno(),msvcrt.LK_UNLCK,1)
        else:
            import fcntl
            try:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except OSError as exc:raise RuntimeError('Another writer holds this output directory') from exc
            try:yield
            finally:fcntl.flock(handle,fcntl.LOCK_UN)
