"""Hash-based reuse, verified Atlas downloads, and resumable frozen ESM2-3B inference."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import sys
import time
from urllib.request import urlopen

import numpy as np
from Bio import SeqIO

SIGNATURE = 'esm2_t36_3B_UR50D:layer36:residue_mean:full_sequence:v1'
DIM = 2560
AA = set('ACDEFGHIKLMNPQRSTVWY')


def sha(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def read_fasta(path, limit=None):
    rows = []; ids = set()
    with Path(path).open(encoding='utf8') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            seq = str(record.seq)
            if record.id in ids:
                raise ValueError('Repeated FASTA identifier: ' + record.id)
            if not seq or not set(seq) <= AA:
                raise ValueError('Expected full standard-amino-acid sequence: ' + record.id)
            ids.add(record.id)
            rows.append(dict(id=record.id, sequence=seq, sequence_sha256=sha(seq), aa_length=len(seq)))
            if limit and len(rows) == limit:break
    if not rows:raise ValueError('Empty FASTA')
    return rows


def valid_vector(vector):
    vector = np.asarray(vector, dtype='<f4')
    if vector.shape != (DIM,) or not np.isfinite(vector).all() or np.linalg.norm(vector) == 0:
        raise ValueError('Invalid ESM2-3B vector (shape, finite values, or zero norm)')
    return vector


class Cache:
    def __init__(self, path):
        self.db = sqlite3.connect(path)
        self.db.execute('CREATE TABLE IF NOT EXISTS vectors (sha TEXT, signature TEXT, vector BLOB NOT NULL, source TEXT NOT NULL, detail TEXT NOT NULL, PRIMARY KEY(sha,signature))')

    def get(self, digest):
        row = self.db.execute('SELECT vector,source,detail FROM vectors WHERE sha=? AND signature=?', (digest,SIGNATURE)).fetchone()
        if row is None:return None
        return valid_vector(np.frombuffer(row[0],dtype='<f4')),row[1],json.loads(row[2])

    def put(self, digest, vector, source, detail):
        vector = valid_vector(vector)
        self.db.execute('INSERT OR IGNORE INTO vectors VALUES (?,?,?,?,?)',
                        (digest,SIGNATURE,vector.tobytes(),source,json.dumps(detail,sort_keys=True)))

    def commit(self):self.db.commit()
    def close(self):self.db.close()


def import_legacy(cache, panel_path, vectors_path, manifest_path, wanted):
    manifest = json.loads(Path(manifest_path).read_text(encoding='utf8'))
    model = manifest.get('model','').lower()
    pooling = manifest.get('pooling','').lower()
    if not ('esm2' in model and '3b' in model and 'mean' in pooling and 'full' in pooling
            and ('no truncation' in pooling or 'full original' in pooling)):
        raise ValueError('Legacy model/pooling manifest is not compatible: ' + str(manifest_path))
    panel = json.loads(Path(panel_path).read_text(encoding='utf8'))
    if not isinstance(panel,list):raise ValueError('Legacy panel must be a JSON array')
    with np.load(vectors_path,allow_pickle=False) as saved:
        ids = list(map(str,saved['ids'])); matrix = saved['embeddings']
        if matrix.shape != (len(ids),DIM) or len(set(ids)) != len(ids):raise ValueError('Legacy vector shape or ID mismatch')
        lookup = {r['id']:r for r in panel}
        if len(lookup) != len(panel) or set(ids) != set(lookup):raise ValueError('Legacy panel IDs do not match vector IDs')
        imported = 0
        for i,identifier in enumerate(ids):
            row = lookup[identifier]; digest = row.get('sequence_sha256')
            if row.get('sequence') and sha(row['sequence']) != digest:raise ValueError('Legacy sequence checksum mismatch: '+identifier)
            if digest in wanted and cache.get(digest) is None:
                cache.put(digest,matrix[i],'prior_local',dict(panel=str(Path(panel_path).resolve()),
                          vectors=str(Path(vectors_path).resolve()),manifest=str(Path(manifest_path).resolve()),id=identifier))
                imported += 1
        cache.commit()
    return imported


def fetch_text(url):
    with urlopen(url,timeout=30) as response:
        data = response.read(2*1024*1024+1)
        if len(data) > 2*1024*1024:raise ValueError('Atlas response exceeded 2 MiB cap')
        return data.decode('utf8')


def fetch_atlas(cache, mapping, wanted, out):
    audit = []
    with Path(mapping).open(encoding='utf8',newline='') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            digest = row['sequence_sha256']; identifier = row['atlas_id']
            if digest not in wanted or cache.get(digest) is not None:continue
            result = dict(sequence_sha256=digest,atlas_id=identifier)
            try:
                if not re.fullmatch(r'MGYP\d+',identifier):raise ValueError('Invalid Atlas identifier')
                raw = fetch_text('https://api.esmatlas.com/fetchSequence/'+identifier+'.fasta')
                lines = raw.splitlines()
                if not lines or not lines[0].startswith('>') or sum(l.startswith('>') for l in lines) != 1:
                    raise ValueError('Expected one Atlas FASTA record')
                seq = ''.join(l.strip() for l in lines[1:])
                if sha(seq) != digest:raise ValueError('Atlas sequence does not exactly match the input')
                payload = json.loads(fetch_text('https://api.esmatlas.com/fetchEmbedding/ESM2/'+identifier+'.json'))
                if payload.get('id') != identifier:raise ValueError('Atlas embedding ID mismatch')
                cache.put(digest,payload['embedding'],'atlas',dict(atlas_id=identifier,sequence_verified=True))
                cache.commit(); result['status']='downloaded_exact_match'
            except Exception as exc:
                result.update(status='not_reused',error=str(exc))
            audit.append(result)
            time.sleep(.25)
    (out/'atlas_download_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
    return audit


def load_model(path, model_format, device):
    # torch/transformers are optional until the inference stage.
    import torch
    from transformers import EsmConfig, EsmModel
    torch.set_num_threads(8)
    if device == 'cuda' and not torch.cuda.is_available():raise RuntimeError('CUDA GPU unavailable')
    torch.backends.cuda.matmul.allow_tf32=False
    config = EsmConfig.from_pretrained(path,local_files_only=True)
    if (config.num_hidden_layers,config.hidden_size) != (36,DIM):raise ValueError('Expected ESM2 t36 3B configuration')
    if model_format == 'esmfold-stem':
        state_all = torch.load(Path(path)/'pytorch_model.bin',map_location='cpu',weights_only=True,mmap=True)
        state = {k[4:]:v for k,v in state_all.items() if k.startswith('esm.') and not k.startswith('esm.contact_head.')}
        state.pop('embeddings.position_ids',None)
        with torch.device('meta'):model=EsmModel(config,add_pooling_layer=False)
        model.contact_head=torch.nn.Identity()
        missing,unexpected=model.load_state_dict(state,strict=False,assign=True)
        if missing or unexpected:raise ValueError('Checkpoint keys mismatch: '+str((missing,unexpected)))
        model.embeddings.position_ids=torch.arange(config.max_position_embeddings).expand((1,-1))
        del state,state_all
    else:
        model=EsmModel.from_pretrained(path,add_pooling_layer=False,local_files_only=True)
    model=model.eval().to(device)
    if device == 'cuda':model=model.half()
    vocab={token:i for i,token in enumerate(config.vocab_list)}
    return model,vocab


def embed_batch(model,vocab,rows,device):
    import torch
    width=max(r['aa_length'] for r in rows)+2
    tokens=torch.full((len(rows),width),vocab['<pad>'],dtype=torch.long,device=device)
    mask=torch.zeros_like(tokens)
    residue_mask=torch.zeros_like(tokens,dtype=torch.bool)
    for j,row in enumerate(rows):
        values=[vocab['<cls>']]+[vocab[a] for a in row['sequence']]+[vocab['<eos>']]
        tokens[j,:len(values)]=torch.tensor(values,device=device)
        mask[j,:len(values)]=1;residue_mask[j,1:len(values)-1]=True
    with torch.inference_mode():
        hidden=model(tokens,attention_mask=mask).last_hidden_state.float()
        vectors=(hidden*residue_mask.unsqueeze(-1)).sum(1)/residue_mask.sum(1,keepdim=True)
    return [valid_vector(v) for v in vectors.cpu().numpy()]


def infer(cache, rows, args):
    import torch
    pending=[r for r in rows if cache.get(r['sequence_sha256']) is None]
    if not pending:return 0
    # Explicitly require a long-sequence policy in a future extension; never truncate.
    if any(r['aa_length'] > args.max_length for r in pending):raise ValueError('Input exceeds --max-length; no sequence was truncated')
    print('Loading model; pending unique sequences:',len(pending),flush=True)
    model,vocab=load_model(args.model_dir,args.model_format,args.device)
    pending.sort(key=lambda r:(r['aa_length'],r['sequence_sha256']))
    done=0;start=time.monotonic();batch_size=args.batch_size
    while done < len(pending):
        batch=pending[done:done+batch_size]
        while len(batch)>1 and len(batch)*(max(r['aa_length'] for r in batch)+2)>args.token_budget:batch=batch[:-1]
        try:
            vectors=embed_batch(model,vocab,batch,args.device)
        except torch.cuda.OutOfMemoryError:
            if len(batch)==1:raise
            batch_size=max(1,len(batch)//2);torch.cuda.empty_cache();continue
        for row,vector in zip(batch,vectors):
            cache.put(row['sequence_sha256'],vector,'inference',dict(model_dir=str(args.model_dir),model_format=args.model_format,
                      precision='fp16' if args.device=='cuda' else 'fp32',pooling='fp32 residue mean',device=args.device))
        cache.commit();done+=len(batch)
        if done == len(pending) or done%256 < len(batch):
            print(f'Embedded {done}/{len(pending)}; elapsed {time.monotonic()-start:.1f}s',flush=True)
    return done


def export(cache,rows,out):
    cached=[cache.get(r['sequence_sha256']) for r in rows]
    if any(v is None for v in cached):raise ValueError('Not all requested sequences have embeddings')
    matrix=np.stack([v[0] for v in cached])
    temp=out/'embeddings.npz.tmp'
    with temp.open('wb') as handle:
        np.savez_compressed(handle,ids=np.array([r['id'] for r in rows]),
                            sequence_sha256=np.array([r['sequence_sha256'] for r in rows]),embeddings=matrix)
    temp.replace(out/'embeddings.npz')
    with (out/'sequence_manifest.tsv').open('w',encoding='utf8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=['id','sequence_sha256','aa_length','source'],delimiter='\t');writer.writeheader()
        for row,value in zip(rows,cached):writer.writerow({k:row[k] for k in ('id','sequence_sha256','aa_length')}|dict(source=value[1]))
    return Counter(v[1] for v in cached)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--stage',choices=('prepare','infer','export','all'),default='all')
    p.add_argument('--reuse',nargs=3,action='append',default=[],metavar=('PANEL_JSON','VECTORS_NPZ','MANIFEST_JSON'))
    p.add_argument('--atlas-map',type=Path,help='Optional TSV: sequence_sha256, atlas_id. Exact sequence verified before reuse.')
    p.add_argument('--model-dir',type=Path,help='Local checkpoint directory, already downloaded')
    p.add_argument('--model-format',choices=('esm2','esmfold-stem'),default='esm2')
    p.add_argument('--device',choices=('cuda','cpu'),default='cuda')
    p.add_argument('--batch-size',type=int,default=16);p.add_argument('--token-budget',type=int,default=4096)
    p.add_argument('--max-length',type=int,default=1022);p.add_argument('--limit',type=int)
    args=p.parse_args(argv)
    if min(args.batch_size,args.token_budget,args.max_length)<1 or (args.limit is not None and args.limit<1):p.error('Limits must be positive')
    if args.stage in ('all','infer') and not args.model_dir:p.error('--model-dir is required for inference')
    args.out.mkdir(parents=True,exist_ok=True)
    rows=read_fasta(args.fasta,args.limit);unique={r['sequence_sha256']:r for r in rows}
    cache=Cache(args.out/'embedding_cache.sqlite')
    record=dict(status='started',signature=SIGNATURE,model='ESM2 t36 3B UR50D',layer=36,dimensions=DIM,
                pooling='Final-layer mean over all residues; BOS/EOS/padding excluded; no truncation',
                input_fasta=str(args.fasta.resolve()),input_sha256=hashlib.sha256(args.fasta.read_bytes()).hexdigest(),
                requested_rows=len(rows),unique_sequences=len(unique),limit=args.limit,
                started_utc=datetime.now(timezone.utc).isoformat(),arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()})
    manifest_path=args.out/('run_'+args.stage+'.json')
    try:
        record['initially_cached']=sum(cache.get(h) is not None for h in unique)
        if args.stage in ('prepare','all'):
            record['imported_legacy']=sum(import_legacy(cache,*bundle,set(unique)) for bundle in args.reuse)
            if args.atlas_map:record['atlas_download_attempts']=len(fetch_atlas(cache,args.atlas_map,set(unique),args.out))
        if args.stage in ('infer','all'):record['computed_this_run']=infer(cache,list(unique.values()),args)
        missing=[r for h,r in unique.items() if cache.get(h) is None]
        with (args.out/'pending.faa').open('w',encoding='ascii') as handle:
            for row in missing:handle.write('>'+row['id']+'\n'+row['sequence']+'\n')
        record['remaining_unique']=len(missing)
        if args.stage in ('infer','export','all'):record['source_counts']=dict(export(cache,rows,args.out))
        record['status']='prepared' if args.stage=='prepare' else 'complete'
    except Exception as exc:
        record.update(status='failed',error=type(exc).__name__+': '+str(exc));raise
    finally:
        record['finished_utc']=datetime.now(timezone.utc).isoformat()
        manifest_path.write_text(json.dumps(record,indent=2),encoding='utf8');cache.close()
    print(json.dumps(record,indent=2),flush=True)


if __name__=='__main__':main()
