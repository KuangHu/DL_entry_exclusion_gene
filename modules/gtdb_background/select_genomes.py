"""Select a reproducible, genus-deduplicated pilot from pinned GTDB metadata."""
from pathlib import Path
from collections import defaultdict,Counter
import argparse,csv,gzip,hashlib,json

def number(r,key):
    try:return float(r[key])
    except (KeyError,ValueError):return None

def select(metadata,out,count=20,seed='gtdb-background-v1'):
    if count<1:raise ValueError('genomes must be positive')
    out=Path(out);out.mkdir(parents=True,exist_ok=True);eligible=[];counts=Counter();total=0
    with gzip.open(metadata,'rt',encoding='utf8',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        fields=reader.fieldnames
        required={'accession','gtdb_taxonomy','gtdb_representative','ncbi_assembly_level','checkm2_completeness','checkm2_contamination'}
        assert required<=set(fields),required-set(fields)
        for r in reader:
            total+=1;reason=''
            if r['gtdb_representative'].lower() not in ('t','true'):reason='not_species_representative'
            elif not r['accession'].startswith('RS_GCF_'):reason='not_RefSeq'
            elif r['ncbi_assembly_level'].lower()!='complete genome':reason='not_complete_assembly'
            elif number(r,'checkm2_completeness') is None or number(r,'checkm2_contamination') is None:reason='missing_quality'
            elif number(r,'checkm2_completeness')<95 or number(r,'checkm2_contamination')>5:reason='quality_threshold'
            elif r.get('ncbi_genome_category','').lower() not in ('none','na','n/a',''):reason='nonisolate_category'
            if reason:counts[reason]+=1;continue
            taxa=dict(t.split('__',1) for t in r['gtdb_taxonomy'].split(';') if '__' in t)
            if not taxa.get('g') or not taxa.get('c'):counts['missing_genus_or_class']+=1;continue
            counts['eligible']+=1
            accession=r['accession'][3:]
            rank=hashlib.sha256((seed+'|'+accession).encode()).hexdigest()
            eligible.append(dict(gtdb_accession=r['accession'],assembly_accession=accession,gtdb_taxonomy=r['gtdb_taxonomy'],phylum=taxa.get('p',''),taxonomic_class=taxa['c'],genus=taxa['g'],species=taxa.get('s',''),checkm2_completeness=r['checkm2_completeness'],checkm2_contamination=r['checkm2_contamination'],ncbi_assembly_level=r['ncbi_assembly_level'],ncbi_organism_name=r.get('ncbi_organism_name',''),ncbi_isolation_source=r.get('ncbi_isolation_source',''),genome_size=r.get('genome_size',''),sampling_rank=rank))
    classes=defaultdict(list)
    for r in eligible:classes[r['taxonomic_class']].append(r)
    for rs in classes.values():rs.sort(key=lambda r:r['sampling_rank'])
    class_order=sorted(classes,key=lambda c:hashlib.sha256((seed+'|class|'+c).encode()).hexdigest())
    chosen=[];genera=set()
    while len(chosen)<count:
        progressed=False
        for cls in class_order:
            while classes[cls] and classes[cls][0]['genus'] in genera:classes[cls].pop(0)
            if not classes[cls]:continue
            r=classes[cls].pop(0);chosen.append(r);genera.add(r['genus']);progressed=True
            if len(chosen)==count:break
        if not progressed:break
    assert len(chosen)==count,('Insufficient eligible genera',len(chosen))
    for name,rows in [('eligible_assemblies.tsv',eligible),('selected_assemblies.tsv',chosen)]:
        with (out/name).open('w',encoding='utf8',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    (out/'selected_assemblies.json').write_text(json.dumps(chosen,indent=2),encoding='utf8')
    manifest=dict(gtdb_release='R232 / 232.0',metadata_file=str(Path(metadata).resolve()),rows_scanned=total,filter_counts=dict(counts),selection='GTDB bacterial species representatives; RefSeq; Complete Genome; CheckM2 completeness >=95, contamination <=5; genome category none/NA; one genome per genus; deterministic round-robin across taxonomic classes.',seed=seed,selected=len(chosen),classes=len({r['taxonomic_class'] for r in chosen}),phyla=len({r['phylum'] for r in chosen}),note='Diversity pilot, not taxonomically matched to EEX positives; sampling weights do not estimate population prevalence.')
    (out/'genome_selection_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
    print(json.dumps(manifest,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--metadata',required=True);p.add_argument('--out',required=True);p.add_argument('--genomes',type=int,default=20);p.add_argument('--seed',default='gtdb-background-v1')
    a=p.parse_args();select(a.metadata,a.out,a.genomes,a.seed)
