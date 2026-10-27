"""Numerically compare padded-batch inference against previously computed references."""
import argparse
import json
from pathlib import Path
import numpy as np
from Bio import SeqIO
from embed import embed_batch, load_model, sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model-dir',type=Path,required=True)
    p.add_argument('--model-format',choices=('esm2','esmfold-stem'),default='esm2')
    p.add_argument('--panel',type=Path,required=True);p.add_argument('--vectors',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--atlas-fasta',type=Path);p.add_argument('--atlas-vector',type=Path)
    args=p.parse_args()
    panel=json.loads(args.panel.read_text(encoding='utf8'))
    with np.load(args.vectors,allow_pickle=False) as saved:
        lookup=dict(zip(map(str,saved['ids']),saved['embeddings']))
    ordered=sorted(panel,key=lambda r:len(r['sequence']))
    selected=[ordered[i] for i in np.linspace(0,len(ordered)-1,6,dtype=int)]
    rows=[dict(id=r['id'],sequence=r['sequence'],sequence_sha256=sha(r['sequence']),aa_length=len(r['sequence'])) for r in selected]
    targets=[lookup[r['id']] for r in selected]
    if bool(args.atlas_fasta) != bool(args.atlas_vector):p.error('Both Atlas calibration files are required')
    if args.atlas_fasta:
        records=list(SeqIO.parse(args.atlas_fasta,'fasta'));payload=json.loads(args.atlas_vector.read_text())
        if len(records)!=1 or records[0].id!=payload['id']:raise ValueError('Atlas calibration IDs disagree')
        seq=str(records[0].seq)
        rows.append(dict(id=records[0].id,sequence=seq,sequence_sha256=sha(seq),aa_length=len(seq)))
        targets.append(np.asarray(payload['embedding'],dtype=np.float32))
    model,vocab=load_model(args.model_dir,args.model_format,'cuda')
    vectors=embed_batch(model,vocab,rows,'cuda')
    results=[]
    for row,expected,got in zip(rows,targets,vectors):
        x=np.asarray(expected,dtype=np.float64);y=got.astype(np.float64)
        cosine=float(x@y/(np.linalg.norm(x)*np.linalg.norm(y)))
        rel=float(np.linalg.norm(x-y)/np.linalg.norm(x))
        results.append(dict(id=row['id'],aa_length=row['aa_length'],cosine_similarity=cosine,
                            relative_l2_error=rel,passed=cosine>=.99999 and rel<=.01))
    report=dict(status='passed' if all(r['passed'] for r in results) else 'failed',
                comparison='Mixed-length padded batch versus legacy single-sequence vectors and optional downloaded Atlas example.',
                model_format=args.model_format,results=results,
                limitation='Numerical compatibility check; not functional or biological validation.')
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2),encoding='utf8');print(json.dumps(report,indent=2))
    if report['status']!='passed':raise SystemExit(1)


if __name__=='__main__':main()
