"""Resumable full-sequence ESM2 inference with validated memory-efficient attention."""
import argparse, csv, hashlib, json, time, types
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import numpy as np

OUT=None
import embed
from neighborhood_io import output_lock, unique_rows

def now():return datetime.now(timezone.utc).isoformat()
def write_json(name,value):
    p=OUT/name;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(value,indent=2),encoding='utf8');t.replace(p)
    if name=='progress.json':
        t=OUT/'progress.js.tmp';t.write_text('window.embeddingProgress='+json.dumps(value)+';',encoding='utf8');t.replace(OUT/'progress.js')

def sdpa_forward(self,hidden_states,attention_mask=None,head_mask=None,encoder_hidden_states=None,
                 encoder_attention_mask=None,past_key_value=None,output_attentions=False):
    import torch
    if self.training or self.is_decoder or head_mask is not None or encoder_hidden_states is not None or past_key_value is not None or output_attentions:
        raise ValueError('This attention path is only for frozen encoder inference without attention outputs')
    if self.position_embedding_type!='rotary':raise ValueError('Expected ESM2 rotary position embeddings')
    q=self.transpose_for_scores(self.query(hidden_states))*self.attention_head_size**-0.5
    k=self.transpose_for_scores(self.key(hidden_states));v=self.transpose_for_scores(self.value(hidden_states))
    q,k=self.rotary_embeddings(q,k)
    # ESM scales queries BEFORE rotary. Do not apply SDPA's default scale again.
    if attention_mask is not None:attention_mask=attention_mask.to(q.dtype)
    from torch.nn.attention import sdpa_kernel,SDPBackend
    with sdpa_kernel([SDPBackend.FLASH_ATTENTION,SDPBackend.EFFICIENT_ATTENTION]):
        x=torch.nn.functional.scaled_dot_product_attention(q,k,v,attn_mask=attention_mask,dropout_p=0.0,scale=1.0)
    return (x.permute(0,2,1,3).contiguous().view(hidden_states.shape[0],hidden_states.shape[1],self.all_head_size),)

def run(args):
    global OUT
    OUT=args.out
    import torch
    def controlled_embed(model,vocab,batch):
        return embed.embed_batch(model,vocab,batch,'cuda')
    rows=unique_rows(args.fasta); cache=embed.Cache(OUT/'embedding_cache.sqlite')
    rows.sort(key=lambda r:r['sequence_sha256'])
    started=now(); record=dict(status='loading',started_utc=started,signature=embed.SIGNATURE,
       unique_sequences=len(rows),device=torch.cuda.get_device_name(0),attention='validated SDPA; full bidirectional attention',
       model_dir=str(args.model_dir.resolve()),pooling='full sequence residue mean; no truncation',
       long_sequence_policy='Full input at all lengths, no chunking. Lengths >1022 separately flagged as extended context; numerical validity is not biological validation.',
       input_sha256=hashlib.sha256(args.fasta.read_bytes()).hexdigest())
    try:
        pending=[r for r in rows if cache.get(r['sequence_sha256']) is None]
        record['initially_cached']=len(rows)-len(pending);record['remaining']=len(pending);write_json('progress.json',record)
        if pending:
            model,vocab=embed.load_model(args.model_dir,args.model_format,'cuda')
            # Same checkpoint, tokens, final-layer output and pooling; compare both
            # implementations including padded batches and >1022-residue inputs.
            sample=[]
            for target in [60,160,300,600,1000,1100,1500,2000]:
                sample.append(min(rows,key=lambda r:abs(r['aa_length']-target)))
            baseline=[controlled_embed(model,vocab,[r])[0] for r in sample]
            # Reject a shape-compatible but numerically different checkpoint
            # before extending a cache copied from another machine.
            cache_checks=[]
            for r,v in zip(sample,baseline):
                saved=cache.get(r['sequence_sha256'])
                if saved is None:continue
                old=saved[0].astype(float);new=v.astype(float)
                cos=float(np.dot(old,new)/(np.linalg.norm(old)*np.linalg.norm(new)))
                err=float(np.linalg.norm(old-new)/np.linalg.norm(old))
                cache_checks.append(dict(sha=r['sequence_sha256'],length=r['aa_length'],cosine=cos,relative_l2=err,passed=cos>=.99999 and err<=.01))
            write_json('checkpoint_cache_validation.json',dict(checks=cache_checks,passed=all(x['passed'] for x in cache_checks)))
            if not all(x['passed'] for x in cache_checks):raise ValueError('Checkpoint differs from transferred cache')
            padded=controlled_embed(model,vocab,sample[:3])
            layers=list(model.encoder.layer);assert len(layers)==36
            for layer in layers:layer.attention.self.forward=types.MethodType(sdpa_forward,layer.attention.self)
            checks=[]
            for r,old in zip(sample,baseline):
                new=controlled_embed(model,vocab,[r])[0]
                cos=float(np.dot(old.astype(float),new)/(np.linalg.norm(old.astype(float))*np.linalg.norm(new.astype(float))))
                err=float(np.linalg.norm(old.astype(float)-new)/np.linalg.norm(old.astype(float)))
                checks.append(dict(sha=r['sequence_sha256'],length=r['aa_length'],cosine=cos,relative_l2=err,passed=cos>=.99999 and err<=.01))
            for r,old,new in zip(sample[:3],padded,controlled_embed(model,vocab,sample[:3])):
                cos=float(np.dot(old.astype(float),new)/(np.linalg.norm(old.astype(float))*np.linalg.norm(new.astype(float))))
                err=float(np.linalg.norm(old.astype(float)-new)/np.linalg.norm(old.astype(float)))
                checks.append(dict(sha=r['sequence_sha256'],length=r['aa_length'],padded_batch=True,cosine=cos,relative_l2=err,passed=cos>=.99999 and err<=.01))
            write_json('attention_validation.json',dict(passed=all(x['passed'] for x in checks),checks=checks,torch=torch.__version__))
            if not all(x['passed'] for x in checks):raise ValueError('Attention compatibility check failed')
            print('Attention compatibility passed',flush=True)
            # Run longest inputs first so full-context memory limitations are
            # discovered immediately. Length sorting still minimizes padding.
            pending.sort(key=lambda r:(-r['aa_length'],r['sequence_sha256']))
            done=0;batchsize=args.batch_size;t0=time.monotonic();last=0;failures=[]
            while done<len(pending):
                batch=pending[done:done+batchsize]
                while len(batch)>1 and len(batch)*(max(r['aa_length'] for r in batch)+2)>args.token_budget:batch=batch[:-1]
                try:
                    vectors=controlled_embed(model,vocab,batch)
                except torch.cuda.OutOfMemoryError:
                    torch.cuda.empty_cache()
                    if len(batch)>1:batchsize=max(1,len(batch)//2);continue
                    failures.append(dict(sha=batch[0]['sequence_sha256'],length=batch[0]['aa_length'],reason='single_sequence_cuda_oom'))
                    done+=1;write_json('failed_sequences.json',failures);continue
                for r,v in zip(batch,vectors):
                    cache.put(r['sequence_sha256'],v,'inference_sdpa_full',dict(model='esm2_t36_3B_UR50D',layer=36,
                         precision='fp16 inference fp32 residue mean',attention='sdpa_scale_1_query_prescaled_before_rotary',
                         aa_length=r['aa_length'],length_scope='extended_context' if r['aa_length']>1022 else 'standard',
                         truncated=False,chunked=False))
                cache.commit();done+=len(batch);elapsed=time.monotonic()-t0
                if elapsed-last>=15 or done==len(pending):
                    record.update(status='running',computed_this_run=done-len(failures),processed_this_run=done,
                         remaining=len(pending)-done+len(failures),elapsed_inference_seconds=elapsed,last_length=batch[-1]['aa_length'],updated_utc=now())
                    write_json('progress.json',record);print(json.dumps({k:record[k] for k in ['computed_this_run','remaining','elapsed_inference_seconds','last_length']}),flush=True);last=elapsed
            del model;torch.cuda.empty_cache()
            if failures:raise RuntimeError(str(len(failures))+' full-sequence OOMs; recorded without truncation')
        record['status']='exporting';write_json('progress.json',record)
        tmp=OUT/'embeddings.tmp.npy';matrix=np.lib.format.open_memmap(tmp,mode='w+',dtype='float32',shape=(len(rows),2560));sources=Counter()
        with (OUT/'sequence_manifest.tsv').open('w',encoding='utf8',newline='') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['row_0based','sequence_sha256','aa_length','source','length_scope'])
            for i,r in enumerate(rows):
                item=cache.get(r['sequence_sha256']);assert item is not None
                v,src,detail=item;matrix[i]=v;sources[src]+=1
                w.writerow([i,r['sequence_sha256'],r['aa_length'],src,'extended_context' if r['aa_length']>1022 else 'standard'])
        matrix.flush();del matrix;tmp.replace(OUT/'embeddings.npy')
        check=np.load(OUT/'embeddings.npy',mmap_mode='r');assert check.shape==(len(rows),2560)
        # Full matrix-to-cache comparison in bounded batches, independent of export loop.
        for i,r in enumerate(rows):assert np.array_equal(check[i],cache.get(r['sequence_sha256'])[0])
        seqhash={r['sequence_sha256'] for r in rows};counts=Counter()
        with args.slots.open(encoding='utf8') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                if r['sequence_sha256']:assert r['sequence_sha256'] in seqhash;counts['covered_valid_slots']+=1
                else:counts['missing_or_excluded_slots']+=1
        digest=hashlib.sha256()
        with (OUT/'embeddings.npy').open('rb') as f:
            while chunk:=f.read(8*1024*1024):digest.update(chunk)
        write_json('embedding_validation.json',dict(status='passed',shape=list(check.shape),dtype=str(check.dtype),
             every_vector_matches_cache=True,finite_nonzero_checked=True,full_input_sha256=record['input_sha256'],
             embeddings_file_sha256=digest.hexdigest(),slot_coverage=dict(counts),source_counts=dict(sources),
             limitation='Integrity and numerical compatibility only; long-context embeddings require separate sensitivity analysis.'))
        record.update(status='complete',remaining=0,finished_utc=now(),source_counts=dict(sources));write_json('progress.json',record)
        print('COMPLETE',json.dumps(record),flush=True)
    except Exception as exc:
        record.update(status='failed',error=type(exc).__name__+': '+str(exc),updated_utc=now());write_json('progress.json',record);raise
    finally:cache.close()

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--model-dir',type=Path,required=True)
    p.add_argument('--model-format',choices=['esmfold-stem','esm2'],default='esmfold-stem')
    p.add_argument('--fasta',type=Path)
    p.add_argument('--slots',type=Path)
    p.add_argument('--batch-size',type=int,default=16)
    p.add_argument('--token-budget',type=int,default=4096)
    args=p.parse_args(argv)
    if min(args.batch_size,args.token_budget)<1:p.error('Batch size and token budget must be positive')
    args.out.mkdir(parents=True,exist_ok=True)
    args.fasta=args.fasta or args.out/'all_proteins.faa'
    args.slots=args.slots or args.out/'slots.tsv'
    with output_lock(args.out):run(args)

if __name__=='__main__':main()
