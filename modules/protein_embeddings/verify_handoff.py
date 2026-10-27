"""Validate an extracted handoff archive before resuming on the cluster; CPU only."""
import argparse,hashlib,json,sqlite3
from contextlib import closing
from pathlib import Path
import numpy as np
from embed import SIGNATURE,valid_vector
from neighborhood_io import unique_rows,output_lock

def file_sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while chunk:=f.read(8*1024*1024):h.update(chunk)
    return h.hexdigest()

def verify(folder):
    folder=Path(folder).resolve();manifest=json.loads((folder/'handoff_manifest.json').read_text())
    for name,expected in manifest['files'].items():
        p=(folder/name).resolve()
        if not p.is_relative_to(folder):raise ValueError('Manifest path escapes handoff directory')
        if file_sha(p)!=expected['sha256']:raise ValueError('Checksum mismatch: '+name)
    rows=unique_rows(folder/'all_proteins.faa');wanted={r['sequence_sha256'] for r in rows}
    with closing(sqlite3.connect((folder/'embedding_cache.sqlite').as_uri()+'?mode=ro',uri=True)) as c:
        if c.execute('pragma quick_check').fetchone()[0]!='ok':raise ValueError('SQLite integrity check failed')
        seen=set()
        for h,sig,raw in c.execute('select sha,signature,vector from vectors'):
            if sig!=SIGNATURE or h not in wanted:raise ValueError('Unexpected cache signature or sequence')
            valid_vector(np.frombuffer(raw,dtype='<f4'));seen.add(h)
    if len(rows)!=manifest['unique_sequences'] or len(seen)!=manifest['saved_vectors']:
        raise ValueError('Handoff counts do not match manifest')
    result=dict(status='passed',unique_sequences=len(rows),saved_vectors=len(seen),remaining=len(rows)-len(seen))
    print(json.dumps(result,indent=2));return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path);a=p.parse_args()
    with output_lock(a.folder):verify(a.folder)

if __name__=='__main__':main()
