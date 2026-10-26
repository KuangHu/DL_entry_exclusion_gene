"""Bounded downloads of public GTDB/NCBI files with atomic completion and checksums."""
from pathlib import Path
from urllib.request import Request,urlopen
import hashlib,json,time

def digest(path,algorithm='sha256'):
    h=hashlib.new(algorithm)
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def download(url,path,max_bytes=400*1024*1024,expected_md5=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        if expected_md5 and digest(path,'md5')!=expected_md5:raise ValueError('Cached checksum mismatch: '+str(path))
        return dict(url=url,path=str(path),bytes=path.stat().st_size,sha256=digest(path),cached=True)
    part=path.with_suffix(path.suffix+'.part')
    last=None
    for attempt in range(3):
        try:
            req=Request(url,headers={'User-Agent':'GTDB-background-audit/0.1','Accept':'*/*'})
            with urlopen(req,timeout=45) as response,part.open('wb') as f:
                total=int(response.headers.get('Content-Length',0))
                if total>max_bytes:raise ValueError('Download exceeds byte limit')
                n=0;t=time.monotonic();start=t
                while True:
                    block=response.read(1024*1024)
                    if not block:break
                    n+=len(block)
                    if n>max_bytes:raise ValueError('Download exceeds byte limit')
                    f.write(block)
                    if time.monotonic()-t>=15:
                        print(f'{path.name}: {n/1024**2:.1f} MiB'+(f' / {total/1024**2:.1f}' if total else ''),flush=True);t=time.monotonic()
                if total and total!=n:raise IOError('Truncated download')
            if expected_md5 and digest(part,'md5')!=expected_md5:raise ValueError('Downloaded MD5 mismatch')
            result=dict(url=url,path=str(path),bytes=n,sha256=digest(part),md5=digest(part,'md5'),seconds=time.monotonic()-start,cached=False)
            part.replace(path);path.with_suffix(path.suffix+'.download.json').write_text(json.dumps(result,indent=2),encoding='utf8')
            print('Downloaded',path.name,n,'bytes',flush=True);return result
        except (OSError,TimeoutError) as e:
            last=e
            if attempt<2:time.sleep(2**attempt)
    raise last

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('url');p.add_argument('path');p.add_argument('--md5');p.add_argument('--max-mib',type=int,default=400)
    a=p.parse_args();print(json.dumps(download(a.url,a.path,a.max_mib*1024*1024,a.md5),indent=2))
