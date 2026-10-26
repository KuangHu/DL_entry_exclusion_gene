"""Fetch exact assembly versions selected by GTDB; retain NCBI annotations unchanged."""
from pathlib import Path
import argparse,concurrent.futures,hashlib,json,re,zipfile
from urllib.parse import urlencode
from download import download,digest

def fetch_one(r,cache):
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
            content=z.read(gb) # ZIP CRC is checked during read. No archive path is used for extraction.
            (folder/'genomic.gbff').write_bytes(content)
            (folder/'assembly_data_report.json').write_text(json.dumps(reports[0],indent=2),encoding='utf8')
            record=dict(**r,status='downloaded',gbff=str((folder/'genomic.gbff').resolve()),gbff_sha256=hashlib.sha256(content).hexdigest(),download_url=url,zip_sha256=d['sha256'],annotation_source='NCBI GenBank/RefSeq annotation; not GTDB Prodigal gene calls')
            (folder/'provenance.json').write_text(json.dumps(record,indent=2),encoding='utf8')
            print(acc,'verified',len(content),'GBFF bytes',flush=True);return record
    except Exception as e:
        print(acc,'FAILED',type(e).__name__,str(e),flush=True)
        return dict(**r,status='download_failed',error=type(e).__name__+': '+str(e),download_url=url)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--selection',required=True);p.add_argument('--cache',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    rows=json.loads(Path(a.selection).read_text(encoding='utf8'))
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(lambda r:fetch_one(r,a.cache),rows))
    Path(a.out).write_text(json.dumps(results,indent=2),encoding='utf8')
    print('Verified',sum(r['status']=='downloaded' for r in results),'/',len(results),flush=True)
    if any(r['status']!='downloaded' for r in results):raise SystemExit(1)
