"""Verify that exported vectors cover the exact input sequences in input order."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from embed import Cache, DIM, SIGNATURE, read_fasta


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();rows=read_fasta(args.fasta)
    with np.load(args.out/'embeddings.npz',allow_pickle=False) as saved:
        ids=list(map(str,saved['ids']));hashes=list(map(str,saved['sequence_sha256']));matrix=saved['embeddings']
    if ids != [r['id'] for r in rows] or hashes != [r['sequence_sha256'] for r in rows]:raise ValueError('Output ID/order/sequence hashes do not match the full input')
    if matrix.shape != (len(rows),DIM) or matrix.dtype != np.float32 or not np.isfinite(matrix).all():raise ValueError('Vector dimensions, dtype or finite-value check failed')
    norms=np.linalg.norm(matrix,axis=1)
    if not (norms>0).all():raise ValueError('Zero-norm vectors')
    with (args.out/'sequence_manifest.tsv').open(encoding='utf8',newline='') as handle:
        manifest=list(csv.DictReader(handle,delimiter='\t'))
    if len(manifest)!=len(rows):raise ValueError('Sequence manifest count mismatch')
    cache=Cache(args.out/'embedding_cache.sqlite');sources=Counter()
    try:
        for row,vec,meta in zip(rows,matrix,manifest):
            cached=cache.get(row['sequence_sha256'])
            if cached is None or not np.array_equal(cached[0],vec):raise ValueError('Cache/output vector mismatch: '+row['id'])
            if meta['id']!=row['id'] or meta['sequence_sha256']!=row['sequence_sha256'] or int(meta['aa_length'])!=row['aa_length'] or meta['source']!=cached[1]:raise ValueError('Sequence provenance mismatch: '+row['id'])
            sources[cached[1]]+=1
    finally:cache.close()
    fingerprint=hashlib.sha256((args.out/'embeddings.npz').read_bytes()).hexdigest()
    report=dict(status='passed',signature=SIGNATURE,input_rows=len(rows),
                unique_sequence_hashes=len(set(hashes)),shape=list(matrix.shape),dtype=str(matrix.dtype),
                input_order_and_full_sequence_hashes_verified=True,all_vectors_finite_and_nonzero=True,
                per_sequence_cache_and_source_verified=True,source_counts=dict(sources),
                norm_min=float(norms.min()),norm_median=float(np.median(norms)),norm_max=float(norms.max()),
                distinct_vector_byte_hashes=len({hashlib.sha256(v.tobytes()).hexdigest() for v in matrix}),
                embeddings_npz_sha256=fingerprint,
                limitation='Representation integrity check only; no EEX function labels or classifier performance assessed.')
    (args.out/'embedding_validation.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
