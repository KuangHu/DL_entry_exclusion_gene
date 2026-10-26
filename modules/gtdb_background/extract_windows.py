"""Sample short-protein genomic windows as unlabeled backgrounds, never confirmed negatives.

Two source modes:
- ncbi-r232: GBFF parsed via Biopython; replicon, topology, and product/gene names available.
- local-r214: GTDB R214 Prodigal protein FASTA headers parsed directly. No product/gene names,
  no replicon, no topology; wrap-around is forced OFF (we never fabricate a circular record
  without authoritative source evidence). Prodigal ordinals are validated for monotonic,
  gap-free order per contig as a consistency check only; the sorted (contig, start, end)
  ordering drives the neighbor index exactly like the NCBI path, so partial/invalid CDS
  retain their position and are never collapsed.
"""
from pathlib import Path
from collections import Counter,defaultdict
import argparse,csv,hashlib,json,re
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition
from fetch_assemblies import PRODIGAL_HEADER

AA=set('ACDEFGHIKLMNPQRSTVWY')
UNCERTAIN=re.compile(r'hypothetical|uncharacteri[sz]ed|unknown|\bDUF\d*|\bUPF\d*',re.I)
EEX_PRODUCT=re.compile(r'(entry|surface|superinfection)\s+exclusion',re.I)
EEX_GENE=re.compile(r'trbK|traS|excA|eex[A-Za-z0-9]*|yddJ',re.I)
CONJ_PRODUCT=re.compile(r'conjugat|conjugal|type[\s-]*(IV|4)\s+secretion|relaxase|type[\s-]*IV\s+coupling',re.I)
CONJ_GENE=re.compile(r'virB\d+|virD4|tra[A-Z]|trb[A-Z]|mob[A-Z]')

def sha(s):return hashlib.sha256(s.encode()).hexdigest()
def save(path,rows,fields):
    with Path(path).open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)

def ordered_cds(record):
    """Keep invalid/missing CDS in the index; never close gaps after filtering."""
    fs=[f for f in record.features if f.type=='CDS']
    if any(f.location is None for f in fs):raise ValueError('Unlocated CDS cannot be ordered: '+record.id)
    return sorted(fs,key=lambda f:(int(f.location.start),int(f.location.end),str(f.qualifiers.get('protein_id',[''])[0])))

def neighborhood_indices(n,center,strand,flank=7,circular=False):
    if strand not in (-1,1):raise ValueError('Center strand must be +1 or -1')
    result=[]
    for rel in range(-flank,flank+1):
        raw=center+rel*strand
        if 0<=raw<n:idx=raw;wrapped=False
        elif circular and n>=2*flank+1:idx=raw%n;wrapped=True
        else:idx=None;wrapped=False
        result.append((rel,idx,wrapped))
    found=[i for _,i,_ in result if i is not None]
    assert len(found)==len(set(found)) # Never fabricate repeated neighbors in tiny circular contigs.
    return result

def decode_feature(feature,index,known):
    q=feature.qualifiers;gene=q.get('gene',[''])[0];product=q.get('product',[''])[0];seq=q.get('translation',[''])[0].replace(' ','').replace('\n','').rstrip('*')
    partial=bool(q.get('partial')) or any(not isinstance(x,ExactPosition) for p in feature.location.parts for x in (p.start,p.end)) or q.get('codon_start',['1'])[0]!='1'
    status='pseudo' if 'pseudo' in q or 'pseudogene' in q else 'missing_translation' if not seq else 'nonstandard_translation' if not set(seq)<=AA else 'partial_CDS' if partial else 'valid'
    unknown=not product.strip() or bool(UNCERTAIN.search(product))
    exact=sha(seq) in known if seq else False
    eex=bool(EEX_PRODUCT.search(product) or EEX_GENE.fullmatch(gene))
    conj=bool(CONJ_PRODUCT.search(product) or CONJ_GENE.fullmatch(gene))
    return dict(cds_index=index,start0=int(feature.location.start),end0=int(feature.location.end),strand=feature.location.strand or 0,
        location=str(feature.location),locus_tag=q.get('locus_tag',[''])[0],protein_accession=q.get('protein_id',[''])[0],
        gene=gene,product=product,aa_length=len(seq),sequence=seq,sequence_sha256=sha(seq) if seq else '',sequence_status=status,
        source_unknown=unknown,exact_project_EEX=exact,exclusion_annotation_flag=eex,conjugation_annotation_flag=conj,
        annotation_review='available')

def _parse_prodigal_gene(header_line,seq_lines,known):
    """Decode a single Prodigal record. Returns (contig, ordinal, gene_dict) where gene_dict
    has the same keys the NCBI decode_feature produces plus Prodigal-specific metadata; the
    contig and ordinal are returned alongside so the caller can group/order without having
    them collide with case-level keys of the same name."""
    m=PRODIGAL_HEADER.match(header_line.rstrip())
    if not m:raise ValueError('Unparseable Prodigal header: '+header_line.rstrip()[:200])
    raw_id=m.group('id')
    suffix=raw_id.rsplit('_protein__',1)[-1] if '_protein__' in raw_id else raw_id
    last_us=suffix.rfind('_')
    if last_us<=0:raise ValueError('Prodigal id lacks contig_ordinal split: '+raw_id)
    contig=suffix[:last_us];ordinal=int(suffix[last_us+1:])
    start_1=int(m.group('start'));end_1=int(m.group('end'));strand=int(m.group('strand'))
    attrs=dict(kv.split('=',1) for kv in m.group('attrs').split(';') if '=' in kv)
    partial_code=attrs.get('partial','00')
    seq=''.join(l.strip() for l in seq_lines).rstrip('*').upper()
    is_partial=partial_code!='00'
    status='missing_translation' if not seq else 'nonstandard_translation' if not set(seq)<=AA else 'partial_CDS' if is_partial else 'valid'
    gene=dict(start0=start_1-1,end0=end_1,strand=strand,
              location=f'{contig}:{start_1-1}-{end_1}({"+"if strand>0 else"-"})',
              locus_tag=raw_id,protein_accession='',
              gene='',product='',
              aa_length=len(seq),sequence=seq,sequence_sha256=sha(seq) if seq else '',sequence_status=status,
              source_unknown=True, # Prodigal provides no product name
              exact_project_EEX=sha(seq) in known if seq else False,
              exclusion_annotation_flag=False,conjugation_annotation_flag=False,
              annotation_review='unavailable',
              prodigal_partial=partial_code,prodigal_start_type=attrs.get('start_type',''),
              prodigal_rbs_motif=attrs.get('rbs_motif',''),prodigal_ordinal=ordinal)
    return contig,ordinal,gene

def _iter_prodigal_records(path,known):
    """Group genes by contig, assign cds_index per sorted order, and verify Prodigal ordinal
    consistency. Yields one record per contig; wrap-around is off; nt_length is approximated
    from the max end0 for display only."""
    by_contig=defaultdict(list)
    ordinals_by_contig=defaultdict(list)
    with Path(path).open('r',encoding='utf8') as f:
        header=None;seq_lines=[]
        for line in f:
            if line.startswith('>'):
                if header is not None:
                    c,o,g=_parse_prodigal_gene(header,seq_lines,known)
                    by_contig[c].append(g);ordinals_by_contig[c].append(o)
                header=line;seq_lines=[]
            else:seq_lines.append(line)
        if header is not None:
            c,o,g=_parse_prodigal_gene(header,seq_lines,known)
            by_contig[c].append(g);ordinals_by_contig[c].append(o)
    ordinal_consistency={}
    for contig,gs in by_contig.items():
        order=sorted(range(len(gs)),key=lambda j:(gs[j]['start0'],gs[j]['end0'],ordinals_by_contig[contig][j]))
        gs_sorted=[gs[j] for j in order]
        ordinals_sorted=[ordinals_by_contig[contig][j] for j in order]
        monotonic=all(b>a for a,b in zip(ordinals_sorted,ordinals_sorted[1:]))
        expected=list(range(ordinals_sorted[0],ordinals_sorted[0]+len(ordinals_sorted))) if ordinals_sorted else []
        gap_free=ordinals_sorted==expected
        ordinal_consistency[contig]=dict(count=len(ordinals_sorted),min=min(ordinals_sorted) if ordinals_sorted else None,max=max(ordinals_sorted) if ordinals_sorted else None,monotonic=monotonic,gap_free=gap_free)
        for i,g in enumerate(gs_sorted):g['cds_index']=i
        by_contig[contig]=gs_sorted
    for contig,gs in by_contig.items():
        yield dict(id=contig,genes=gs,circular=False,replicon='unspecified',description='',
                   nt_length=max((g['end0'] for g in gs),default=0),
                   ordinal_consistency=ordinal_consistency[contig])

def _iter_ncbi_records(path,known):
    """Decode an NCBI GBFF via Biopython using the original sorted-features logic."""
    for record in SeqIO.parse(path,'genbank'):
        fs=ordered_cds(record);genes=[decode_feature(f,i,known) for i,f in enumerate(fs)]
        circular=record.annotations.get('topology','').lower()=='circular'
        source=next((f for f in record.features if f.type=='source'),None)
        replicon='plasmid' if source is not None and 'plasmid' in source.qualifiers else 'chromosome' if source is not None and 'chromosome' in source.qualifiers else 'unspecified'
        yield dict(id=record.id,genes=genes,circular=circular,replicon=replicon,description=record.description,
                   nt_length=len(record),ordinal_consistency=None)

def length_bin(n,minimum=40,maximum=300):
    edges=[minimum]+[x for x in (80,120,200) if minimum<x<=maximum]+[maximum+1]
    return next(f'{a:04d}-{b-1:04d}' for a,b in zip(edges,edges[1:]) if a<=n<b)
def clean_gene(g):
    return {k:v for k,v in g.items() if k!='sequence'}
def arch_token(g):
    if g['source_unknown']:return 'unknown'
    return g['gene'] or g['product']

def extract(download_manifest,out,positive_panel=None,minimum=40,maximum=300,per_genome=20,flank=7,seed='gtdb-background-v1'):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    records=json.loads(Path(download_manifest).read_text(encoding='utf8'))
    if not any(r['status']=='downloaded' for r in records):raise ValueError('No verified assemblies in download manifest')
    known=set()
    if positive_panel:
        panel=json.loads(Path(positive_panel).read_text(encoding='utf8'))
        known={r['sequence_sha256'] for r in panel if r.get('comparison_group')=='EEX_existing_label'}
    cases=[];audit=[];counts=Counter();genome_summary=[];sequences={};seen_centers=set();groups={}
    source_modes=Counter();ordinal_warnings=[]
    for assembly in records:
        if assembly['status']!='downloaded':continue
        source_mode=assembly.get('source_mode','ncbi-r232');source_modes[source_mode]+=1
        acc=assembly['assembly_accession'];candidates=[];contigs={};all_cds=0;local_audit=[]
        if source_mode=='local-r214':
            iterator=_iter_prodigal_records(assembly['protein_faa'],known)
        else:
            iterator=_iter_ncbi_records(assembly['gbff'],known)
        for rec in iterator:
            genes=rec['genes'];all_cds+=len(genes)
            contigs[rec['id']]=dict(genes=genes,circular=rec['circular'],replicon=rec['replicon'],description=rec['description'],nt_length=rec['nt_length'],ordinal_consistency=rec['ordinal_consistency'])
            if rec['ordinal_consistency'] and not rec['ordinal_consistency']['gap_free']:
                ordinal_warnings.append(dict(assembly_accession=acc,contig=rec['id'],**rec['ordinal_consistency']))
            for g in genes:
                if g['sequence_status']!='valid':counts['CDS_'+g['sequence_status']]+=1;continue
                if not minimum<=g['aa_length']<=maximum:counts['CDS_outside_center_length']+=1;continue
                reason='eligible'
                if g['strand'] not in (-1,1):reason='unknown_strand'
                elif g['exact_project_EEX']:reason='exact_project_EEX_review'
                elif g['exclusion_annotation_flag']:reason='source_exclusion_annotation_review'
                else:
                    indexes=neighborhood_indices(len(genes),g['cds_index'],g['strand'],flank,rec['circular'])
                    if any(i is None for _,i,_ in indexes):reason='incomplete_CDS_window'
                ar=dict(assembly_accession=acc,contig=rec['id'],center_cds=g['cds_index'],
                        protein_accession=g.get('protein_accession',''),locus_tag=g.get('locus_tag',''),
                        aa_length=g['aa_length'],sequence_sha256=g['sequence_sha256'],
                        gene=g['gene'],product=g['product'],source_unknown=g['source_unknown'],
                        annotation_review=g['annotation_review'],reason=reason)
                local_audit.append(ar);counts['short_'+reason]+=1
                if reason!='eligible':continue
                rank=sha(seed+'|'+acc+'|'+rec['id']+'|'+str(g['cds_index'])+'|'+g['sequence_sha256'])
                candidates.append(dict(contig=rec['id'],gene=g,rank=rank,audit=ar,indexes=indexes))
        bins=defaultdict(list)
        for c in sorted(candidates,key=lambda c:c['rank']):bins[length_bin(c['gene']['aa_length'],minimum,maximum)].append(c)
        used=defaultdict(set);selected=[]
        while len(selected)<per_genome:
            progressed=False
            for key in sorted(bins):
                while bins[key]:
                    c=bins[key].pop(0);g=c['gene'];occupied={i for _,i,_ in c['indexes']}
                    if g['sequence_sha256'] in seen_centers:c['audit']['reason']='exact_duplicate_selected_center';continue
                    if occupied&used[c['contig']]:c['audit']['reason']='overlaps_selected_CDS_window';continue
                    c['audit']['reason']='sampled';selected.append(c);seen_centers.add(g['sequence_sha256']);used[c['contig']].update(occupied);progressed=True;break
                if len(selected)==per_genome:break
            if not progressed:break
        for c in candidates:
            if c['audit']['reason']=='eligible':c['audit']['reason']='not_selected_per_genome_cap'
        for c in selected:
            cid=f'BG{len(cases)+1:06d}';g=c['gene'];ct=contigs[c['contig']];genes=[]
            for rel,i,wrapped in c['indexes']:
                raw=ct['genes'][i];item=clean_gene(raw)
                item.update(relative_cds=rel,relative_strand=raw['strand']*g['strand'],wrapped_origin=wrapped)
                if raw['sequence_status']=='valid':
                    sequences[raw['sequence_sha256']]=raw['sequence']
                    item['sequence_id']='BP'+raw['sequence_sha256'][:20]
                else:item['sequence_id']=''
                genes.append(item)
            context_flag=any(x['exclusion_annotation_flag'] or x['exact_project_EEX'] or x['conjugation_annotation_flag'] for x in genes)
            evidence='source_function_annotated' if not g['source_unknown'] else 'source_function_uncharacterized'
            identity=[(arch_token(x),x['relative_strand']) for x in genes]
            aid='GBGA'+sha(json.dumps(identity,ensure_ascii=True))[:14]
            case=dict(case_id=cid,architecture_id=aid,assembly_accession=acc,contig=c['contig'],center_cds=g['cds_index'],center_strand=g['strand'],
                center_sequence_id='BP'+g['sequence_sha256'][:20],center_sha256=g['sequence_sha256'],center_aa_length=g['aa_length'],
                center_gene=g['gene'],center_product=g['product'],center_evidence=evidence,
                source_annotation_review=g['annotation_review'],
                context_annotation_review=context_flag,negative_label=None,experimental_negative=False,
                gtdb_taxonomy=assembly['gtdb_taxonomy'],genus=assembly['genus'],taxonomic_class=assembly['taxonomic_class'],phylum=assembly['phylum'],
                split_group=assembly['genus'],replicon=ct['replicon'],contig_circular=ct['circular'],contig_nt_length=ct['nt_length'],
                full_CDS_window=True,valid_neighbor_sequences=sum(x['sequence_status']=='valid' for x in genes if x['relative_cds']!=0),
                wrapped_origin=any(x['wrapped_origin'] for x in genes),
                hydrophobic_fraction=sum(g['sequence'].count(a) for a in 'AVILMFWY')/g['aa_length'],
                charge_proxy=(g['sequence'].count('K')+g['sequence'].count('R')-g['sequence'].count('D')-g['sequence'].count('E'))/g['aa_length'],
                annotation_source=assembly['annotation_source'],
                source_mode=source_mode,
                source_file_sha256=assembly.get('gbff_sha256') or assembly.get('protein_faa_sha256',''),
                genes=genes)
            if source_mode=='ncbi-r232':
                case['source_gbff_sha256']=assembly['gbff_sha256'] # Preserve the original provenance column.
            cases.append(case)
            if aid not in groups:groups[aid]=dict(architecture_id=aid,cases=[],center_frequencies=Counter(),positions=[dict(relative_cds=x['relative_cds'],identity_label=arch_token(x),relative_strand=x['relative_strand'],valid=0,excluded=0,missing=0,seqs=set()) for x in genes])
            ag=groups[aid];ag['cases'].append(cid);ag['center_frequencies'][g['sequence_sha256']]+=1
            for x,p in zip(genes,ag['positions']):
                if x['sequence_status']=='valid':p['valid']+=1;p['seqs'].add(x['sequence_sha256'])
                else:p['excluded']+=1
        audit.extend(local_audit)
        genome_summary.append(dict(assembly_accession=acc,source_mode=source_mode,genus=assembly['genus'],taxonomic_class=assembly['taxonomic_class'],phylum=assembly['phylum'],contigs=len(contigs),CDS_count=all_cds,eligible_short_centers=len(candidates),selected_centers=len(selected),requested_centers=per_genome))
        print(acc,source_mode,'CDS',all_cds,'eligible centers',len(candidates),'sampled',len(selected),flush=True)
    arch=[]
    for ag in groups.values():
        freq=ag.pop('center_frequencies');ag['copies']=len(ag['cases']);ag['distinct_center_aa']=len(freq)
        ag['center_sequence_frequencies']=dict(freq);ag['center_frequency_distribution']=dict(Counter(freq.values()))
        ag['architecture_identity_note']='Positions aggregate ordered (identity_label, relative_strand) tuples. When identity_label is "unknown" the arrow pattern reflects strand topology only; identical all-unknown patterns do not imply a shared functional operon.'
        for p in ag['positions']:p['distinct_aa']=len(p.pop('seqs'));assert p['valid']+p['excluded']+p['missing']==ag['copies']
        arch.append(ag)
    if not cases:raise ValueError('No eligible sampled windows; inspect source annotations and thresholds')
    assert len(cases)==len({c['center_sha256'] for c in cases})
    assert len({'BP'+s[:20] for s in sequences})==len(sequences)
    (out/'windows.jsonl').write_text(''.join(json.dumps(c,ensure_ascii=False)+'\n' for c in cases),encoding='utf8')
    (out/'architectures.json').write_text(json.dumps(arch,ensure_ascii=False,indent=2),encoding='utf8')
    fields=[k for k in cases[0] if k!='genes'] if cases else ['case_id']
    save(out/'centers.tsv',[{k:c[k] for k in fields} for c in cases],fields)
    auditfields=list(audit[0]) if audit else ['reason']
    save(out/'center_sampling_audit.tsv',audit,auditfields)
    save(out/'genome_summary.tsv',genome_summary,list(genome_summary[0]))
    neighbors=[]
    for c in cases:
        for g in c['genes']:neighbors.append(dict(case_id=c['case_id'],architecture_id=c['architecture_id'],assembly_accession=c['assembly_accession'],contig=c['contig'],**g))
    save(out/'neighbors.tsv',neighbors,list(neighbors[0]))
    save(out/'architecture_position_diversity.tsv',[dict(architecture_id=a['architecture_id'],copies=a['copies'],**p) for a in arch for p in a['positions']],['architecture_id','copies','relative_cds','identity_label','relative_strand','valid','excluded','missing','distinct_aa'])
    (out/'center_proteins.faa').write_text(''.join('>'+c['case_id']+'\n'+sequences[c['center_sha256']]+'\n' for c in cases),encoding='ascii')
    (out/'window_proteins_unique.faa').write_text(''.join('>BP'+s[:20]+'\n'+seq+'\n' for s,seq in sorted(sequences.items())),encoding='ascii')
    source_mode=next(iter(source_modes)) if len(source_modes)==1 else 'mixed'
    GTDB_release='R214 / 214.0' if source_mode=='local-r214' else 'R232 / 232.0'
    annotation_source=('GTDB R214 Prodigal protein predictions; no NCBI functional annotations.'
                       if source_mode=='local-r214'
                       else 'NCBI exact-version RefSeq assemblies, not GTDB Prodigal protein calls')
    source_review_summary='unavailable' if source_mode=='local-r214' else 'available'
    manifest=dict(status='complete',GTDB_release=GTDB_release,source_mode=source_mode,
                  annotation_source=annotation_source,
                  source_annotation_review=source_review_summary,
                  assemblies_requested=len(records),assemblies_downloaded=sum(r['status']=='downloaded' for r in records),
                  assemblies_with_cases=len({c['assembly_accession'] for c in cases}),case_count=len(cases),
                  unique_center_sequences=len(seen_centers),unique_window_sequences=len(sequences),
                  architecture_count=len(arch),flank_CDS=flank,min_aa=minimum,max_aa=maximum,per_genome_cap=per_genome,seed=seed,
                  center_evidence_counts=dict(Counter(c['center_evidence'] for c in cases)),
                  context_review_cases=sum(c['context_annotation_review'] for c in cases),
                  cases_with_unavailable_neighbor_sequences=sum(c['valid_neighbor_sequences']<2*flank for c in cases),
                  source_annotation_audit_counts=dict(Counter(r['reason'] for r in audit)),
                  CDS_filter_counts=dict(counts),
                  known_EEX_exact_sequences_checked=len(known),homology_screen='not_run',
                  negative_label='unlabeled; none are experimentally established negatives',
                  sampling='Deterministic length-bin round robin within each genome; selected center sequences exact-deduplicated globally; selected windows do not share CDS within a source contig.',
                  coordinates=f"All original CDS retained in ordered index, including pseudo/partial/untranslated CDS. Center oriented forward. Circular wrapping only on source-declared circular records with >={2*flank+1} CDS; no crossing contigs.",
                  architecture_identity='Ordered original gene names/products, unknown identity token for uncharacterized source annotations, and relative strands; locus tags and full products retained for each observed case. Evidence categories do not enter identity.',
                  copy_definition='One sampled center occurrence; not an independent biological replicate. Selected center exact duplicates are excluded. Neighbor exact duplicates retained in occurrence counts.',
                  limitations=['Diversity pilot is not taxonomically matched to EEX positives.','No EEX-phenotype negatives are asserted.','Absence of a source annotation keyword is not evidence of biological absence.','Near-homolog clusters and final train/test separation remain pending.','Center SignalP, topology and ESM features have not been computed in this sampling stage.'],
                  annotations_changed=False)
    if source_mode=='local-r214':
        manifest['local_r214_notes']=[
            'Centers and neighbors are Prodigal predictions, not NCBI annotated CDS; set membership differs.',
            'No functional product names: every center is source_function_uncharacterized; architectures reduce to strand patterns of "unknown" tokens.',
            'Wrap-around is disabled; linear-contig edges produce incomplete_CDS_window and no eligibility.',
            'Prodigal ordinals were validated per contig for monotonic gap-free order (see ordinal_warnings). They are a consistency check only; neighbor indexing uses the sorted (start, end, ordinal) order.',
            'source_annotation_review=unavailable; absence of a product name is NOT a signal of EEX exclusion or presence.',
            ]
        manifest['ordinal_warnings']=ordinal_warnings
    (out/'extraction_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8');print(json.dumps(manifest,indent=2),flush=True)
    return cases

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--downloads',required=True);p.add_argument('--out',required=True);p.add_argument('--known-positive-panel');p.add_argument('--min-aa',type=int,default=40);p.add_argument('--max-aa',type=int,default=300);p.add_argument('--per-genome',type=int,default=20);p.add_argument('--flank',type=int,default=7);p.add_argument('--seed',default='gtdb-background-v1')
    a=p.parse_args()
    if not 1<=a.min_aa<=a.max_aa or a.flank<1 or a.per_genome<1:raise SystemExit('Require 1 <= min-aa <= max-aa, flank >=1 and per-genome >=1.')
    extract(a.downloads,a.out,a.known_positive_panel,a.min_aa,a.max_aa,a.per_genome,a.flank,a.seed)
