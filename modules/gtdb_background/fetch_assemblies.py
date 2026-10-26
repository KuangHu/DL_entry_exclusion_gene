"""Fetch exact assembly versions selected by GTDB; retain source annotations unchanged.

Two modes:
- ncbi-r232: downloads NCBI RefSeq GBFF via Datasets API; verifies ZIP CRC + GBFF SHA256.
- local-r214: resolves local GTDB R214 Prodigal protein FASTA files from a user-specified
  faa-root. For each selected accession, scans bacteria/ and archaea/ for files matching
  <prefix><acc>_protein.faa with an optional leading "<n>-" shard tag. Records raw file
  SHA256; computes a canonical biological SHA256 over parsed Prodigal records so that
  shard-prefix duplicates with identical biological content dedup to one, while files with
  genuinely different biological content are reported as content_conflict without picking.
"""
from pathlib import Path
import argparse,concurrent.futures,hashlib,json,re,zipfile
from urllib.parse import urlencode
from download import download,digest

PRODIGAL_HEADER=re.compile(r'^>(?P<id>\S+)\s*#\s*(?P<start>-?\d+)\s*#\s*(?P<end>-?\d+)\s*#\s*(?P<strand>-?1)\s*#\s*(?P<attrs>.*)$')
FAA_NAME=re.compile(r'^(?P<shard>\d+-)?(?P<base>(?:RS_GCF|GB_GCA)_\d+\.\d+)_protein\.faa$')

def _parse_prodigal_file(path):
    """Return list of (contig, ordinal, start, end, strand, partial, seq_sha256) sorted by (contig,ordinal)."""
    records=[]
    with Path(path).open('r',encoding='utf8') as f:
        current=None;seq_parts=[]
        for line in f:
            if line.startswith('>'):
                if current is not None:
                    seq=''.join(seq_parts).rstrip('*').upper()
                    records.append(current+(hashlib.sha256(seq.encode()).hexdigest(),))
                m=PRODIGAL_HEADER.match(line.rstrip())
                if not m:raise ValueError('Unparseable Prodigal header in '+str(path)+': '+line.rstrip()[:200])
                raw_id=m.group('id')
                suffix=raw_id.rsplit('_protein__',1)[-1] if '_protein__' in raw_id else raw_id
                last_us=suffix.rfind('_')
                if last_us<=0:raise ValueError('Prodigal id lacks contig_ordinal split: '+raw_id)
                contig=suffix[:last_us];ordinal_str=suffix[last_us+1:]
                if not ordinal_str.isdigit():raise ValueError('Non-numeric Prodigal ordinal in '+raw_id)
                ordinal=int(ordinal_str)
                attrs=dict(kv.split('=',1) for kv in m.group('attrs').split(';') if '=' in kv)
                partial=attrs.get('partial','00')
                current=(contig,ordinal,int(m.group('start')),int(m.group('end')),int(m.group('strand')),partial)
                seq_parts=[]
            else:
                seq_parts.append(line.strip())
        if current is not None:
            seq=''.join(seq_parts).rstrip('*').upper()
            records.append(current+(hashlib.sha256(seq.encode()).hexdigest(),))
    records.sort(key=lambda r:(r[0],r[1]))
    return records

def _biological_sha256(records):
    h=hashlib.sha256()
    for r in records:h.update(('|'.join(map(str,r))+'\n').encode())
    return h.hexdigest()

def fetch_one_ncbi(r,cache):
    acc=r['assembly_accession'];assert re.fullmatch(r'GCF_\d+\.\d+',acc)
    folder=Path(cache)/acc;folder.mkdir(parents=True,exist_ok=True)
    url='https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'+acc+'/download?'+urlencode({'include_annotation_type':'GENOME_GBFF','filename':acc+'.zip'})
    try:
        d=download(url,folder/'ncbi_dataset.zip',max_bytes=40*1024*1024)
        with zipfile.ZipFile(folder/'ncbi_dataset.zip') as z:
            gb='ncbi_dataset/data/'+acc+'/genomic.gbff';report='ncbi_dataset/data/assembly_data_report.jsonl'
            if gb not in z.namelist():raise ValueError('Exact assembly GBFF not supplied')
            if z.getinfo(gb).file_size>150*1024*1024:raise ValueError('Uncompressed GBFF exceeds cap')
            reports=[json.loads(line) for line in z.read(report).decode().splitlines() if line.strip()]
            if len(reports)!=1 or reports[0].get('accession')!=acc:raise ValueError('Assembly version mismatch')
            content=z.read(gb)
            (folder/'genomic.gbff').write_bytes(content)
            (folder/'assembly_data_report.json').write_text(json.dumps(reports[0],indent=2),encoding='utf8')
            record={**r,'status':'downloaded','source_mode':'ncbi-r232','gbff':str((folder/'genomic.gbff').resolve()),'gbff_sha256':hashlib.sha256(content).hexdigest(),'download_url':url,'zip_sha256':d['sha256'],'annotation_source':'NCBI GenBank/RefSeq annotation; not GTDB Prodigal gene calls'}
            (folder/'provenance.json').write_text(json.dumps(record,indent=2),encoding='utf8')
            print(acc,'verified',len(content),'GBFF bytes',flush=True);return record
    except Exception as e:
        print(acc,'FAILED',type(e).__name__,str(e),flush=True)
        return {**r,'status':'download_failed','source_mode':'ncbi-r232','error':type(e).__name__+': '+str(e),'download_url':url}

def _build_local_r214_index(faa_root):
    """Scan faa-root/{bacteria,archaea} once; return {<prefix><acc>: [resolved_path, ...]}."""
    index={}
    for sub in ('bacteria','archaea'):
        d=Path(faa_root)/sub
        if not d.is_dir():continue
        for p in d.iterdir():
            m=FAA_NAME.match(p.name)
            if not m:continue
            key=m.group('base')
            index.setdefault(key,[]).append(p.resolve())
    return index

def resolve_local_r214(r,index,cache):
    """Look up precomputed index; verify biological content; report conflicts without picking."""
    acc=r['assembly_accession'];prefix=r.get('accession_prefix','')
    key=prefix+acc
    matches=sorted(index.get(key,[]))
    if not matches:
        return {**r,'status':'download_failed','source_mode':'local-r214','error':'FileNotFoundError: no local R214 FAA for '+key+'_protein.faa (checked bacteria/, archaea/)'}
    try:
        parsed=[(p,_parse_prodigal_file(p)) for p in matches]
    except Exception as e:
        return {**r,'status':'download_failed','source_mode':'local-r214','error':type(e).__name__+': '+str(e)}
    bio_hashes={_biological_sha256(recs) for _,recs in parsed}
    file_hashes=[(str(p),digest(p)) for p,_ in parsed]
    if len(bio_hashes)>1:
        return {**r,'status':'content_conflict','source_mode':'local-r214',
                'conflicting_paths':[str(p) for p,_ in parsed],
                'conflicting_biological_sha256':sorted(bio_hashes),
                'error':'Multiple R214 FAA files for '+acc+' disagree on biological content; must be resolved manually'}
    bio_sha=bio_hashes.pop();chosen=max(parsed,key=lambda pr:len(pr[0].name)) # prefer shard-prefixed (longer) name as canonical
    folder=Path(cache)/acc;folder.mkdir(parents=True,exist_ok=True)
    record_count=len(chosen[1])
    record={**r,'status':'downloaded','source_mode':'local-r214',
            'protein_faa':str(chosen[0]),
            'protein_faa_sha256':digest(chosen[0]),
            'biological_sha256':bio_sha,
            'record_count':record_count,
            'duplicate_paths':[p for p,_ in file_hashes if p!=str(chosen[0])],
            'duplicate_sha256':{p:h for p,h in file_hashes if p!=str(chosen[0])},
            'annotation_source':'GTDB R214 Prodigal gene calls; no NCBI functional annotations'}
    (folder/'provenance.json').write_text(json.dumps(record,indent=2),encoding='utf8')
    print(acc,'resolved',record_count,'Prodigal CDS in',chosen[0].name,'(duplicates deduplicated:',len(record['duplicate_paths']),')',flush=True)
    return record

def resolve_many(rows,cache,source,faa_root=None):
    if source=='ncbi-r232':
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:return list(pool.map(lambda r:fetch_one_ncbi(r,cache),rows))
    if source=='local-r214':
        if not faa_root:raise ValueError('local-r214 requires --faa-root')
        index=_build_local_r214_index(faa_root)
        print('Local R214 index built:',len(index),'accessions indexed from',faa_root,flush=True)
        return [resolve_local_r214(r,index,cache) for r in rows]
    raise ValueError('Unknown source: '+source)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--selection',required=True)
    p.add_argument('--cache',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--source',choices=('ncbi-r232','local-r214'),default='ncbi-r232')
    p.add_argument('--faa-root',help='Root for local-r214 containing bacteria/ and archaea/ subfolders')
    a=p.parse_args()
    rows=json.loads(Path(a.selection).read_text(encoding='utf8'))
    results=resolve_many(rows,a.cache,a.source,a.faa_root)
    Path(a.out).write_text(json.dumps(results,indent=2),encoding='utf8')
    print('Verified',sum(r['status']=='downloaded' for r in results),'/',len(results),flush=True)
    if any(r['status']!='downloaded' for r in results):raise SystemExit(1)
