"""Select a reproducible, genus-deduplicated pilot from pinned GTDB metadata."""
from pathlib import Path
from collections import defaultdict,Counter
import argparse,csv,gzip,hashlib,json

SOURCES={
    'ncbi-r232':dict(release='R232 / 232.0',completeness_field='checkm2_completeness',contamination_field='checkm2_contamination',quality_label='CheckM2',accept_prefix=('RS_',),annotation_source='NCBI exact-version RefSeq assemblies, not GTDB Prodigal protein calls'),
    'local-r214':dict(release='R214 / 214.0',completeness_field='checkm_completeness',contamination_field='checkm_contamination',quality_label='CheckM',accept_prefix=('RS_','GB_'),annotation_source='GTDB Prodigal protein calls on R214 reps; no NCBI functional annotations'),
}

def number(r,key):
    try:return float(r[key])
    except (KeyError,ValueError):return None

def _open(metadata):
    p=Path(metadata)
    return gzip.open(p,'rt',encoding='utf8',newline='') if p.suffix=='.gz' else p.open('r',encoding='utf8',newline='')

def select(metadata,out,count=20,seed='gtdb-background-v1',source='ncbi-r232'):
    if count<1:raise ValueError('genomes must be positive')
    if source not in SOURCES:raise ValueError('Unknown source: '+source)
    spec=SOURCES[source]
    out=Path(out);out.mkdir(parents=True,exist_ok=True);eligible=[];counts=Counter();total=0
    with _open(metadata) as f:
        reader=csv.DictReader(f,delimiter='\t')
        fields=reader.fieldnames
        required={'accession','gtdb_taxonomy','gtdb_representative','ncbi_assembly_level',spec['completeness_field'],spec['contamination_field']}
        assert required<=set(fields),required-set(fields)
        for r in reader:
            total+=1;reason=''
            if r['gtdb_representative'].lower() not in ('t','true'):reason='not_species_representative'
            elif not any(r['accession'].startswith(p) for p in spec['accept_prefix']):reason='accession_prefix_not_accepted'
            elif r['ncbi_assembly_level'].lower()!='complete genome':reason='not_complete_assembly'
            elif number(r,spec['completeness_field']) is None or number(r,spec['contamination_field']) is None:reason='missing_quality'
            elif number(r,spec['completeness_field'])<95 or number(r,spec['contamination_field'])>5:reason='quality_threshold'
            elif r.get('ncbi_genome_category','').lower() not in ('none','na','n/a',''):reason='nonisolate_category'
            if reason:counts[reason]+=1;continue
            taxa=dict(t.split('__',1) for t in r['gtdb_taxonomy'].split(';') if '__' in t)
            if not taxa.get('g') or not taxa.get('c'):counts['missing_genus_or_class']+=1;continue
            counts['eligible']+=1
            accession_prefix=next(p for p in spec['accept_prefix'] if r['accession'].startswith(p))
            assembly_accession=r['accession'][len(accession_prefix):]
            # Preserve the original R232 sampling order; R214 uses its source prefix.
            rank_accession=assembly_accession if source=='ncbi-r232' else r['accession']
            rank=hashlib.sha256((seed+'|'+rank_accession).encode()).hexdigest()
            row=dict(gtdb_accession=r['accession'],assembly_accession=assembly_accession,accession_prefix=accession_prefix,
                     gtdb_taxonomy=r['gtdb_taxonomy'],phylum=taxa.get('p',''),taxonomic_class=taxa['c'],
                     genus=taxa['g'],species=taxa.get('s',''),
                     quality_label=spec['quality_label'],
                     completeness=r[spec['completeness_field']],contamination=r[spec['contamination_field']],
                     ncbi_assembly_level=r['ncbi_assembly_level'],
                     ncbi_organism_name=r.get('ncbi_organism_name',''),ncbi_isolation_source=r.get('ncbi_isolation_source',''),
                     genome_size=r.get('genome_size',''),sampling_rank=rank,source_mode=source,
                     ncbi_genbank_assembly_accession=r.get('ncbi_genbank_assembly_accession',''))
            eligible.append(row)
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
    manifest=dict(gtdb_release=spec['release'],source_mode=source,quality_label=spec['quality_label'],
                  metadata_file=str(Path(metadata).resolve()),rows_scanned=total,filter_counts=dict(counts),
                  accepted_accession_prefixes=list(spec['accept_prefix']),annotation_source=spec['annotation_source'],
                  selection=f"GTDB bacterial species representatives; accession prefix in {list(spec['accept_prefix'])}; Complete Genome; {spec['quality_label']} completeness >=95, contamination <=5; genome category none/NA; one genome per genus; deterministic round-robin across taxonomic classes.",
                  seed=seed,selected=len(chosen),classes=len({r['taxonomic_class'] for r in chosen}),phyla=len({r['phylum'] for r in chosen}),
                  note='Diversity pilot, not taxonomically matched to EEX positives; sampling weights do not estimate population prevalence.')
    (out/'genome_selection_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
    print(json.dumps(manifest,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--metadata',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--genomes',type=int,default=20)
    p.add_argument('--seed',default='gtdb-background-v1')
    p.add_argument('--source',choices=sorted(SOURCES),default='ncbi-r232')
    a=p.parse_args();select(a.metadata,a.out,a.genomes,a.seed,a.source)
