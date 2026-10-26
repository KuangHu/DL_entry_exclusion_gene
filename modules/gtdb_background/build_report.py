"""Validate source windows and render a directly viewable background-sampling atlas."""
from pathlib import Path
from collections import Counter,defaultdict
import argparse,csv,hashlib,html,json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
from extract_windows import ordered_cds,neighborhood_indices,decode_feature

def build(out):
    out=Path(out);manifest=json.loads((out/'extraction_manifest.json').read_text())
    selection=json.loads((out/'genome_selection_manifest.json').read_text())
    windows=[json.loads(l) for l in (out/'windows.jsonl').read_text(encoding='utf8').splitlines()]
    arch={a['architecture_id']:a for a in json.loads((out/'architectures.json').read_text())}
    downloads=json.loads((out/'download_manifest.json').read_text())
    fasta={r.id:str(r.seq) for r in SeqIO.parse(out/'window_proteins_unique.faa','fasta')}
    central={r.id:str(r.seq) for r in SeqIO.parse(out/'center_proteins.faa','fasta')}
    byassembly=defaultdict(list)
    for c in windows:byassembly[c['assembly_accession']].append(c)
    used=set();seqstatuses=Counter();positions=0;wrapcases=0
    for r in downloads:
        if r['status']!='downloaded':continue
        assert hashlib.sha256(Path(r['gbff']).read_bytes()).hexdigest()==r['gbff_sha256']
        sources={rec.id:rec for rec in SeqIO.parse(r['gbff'],'genbank')}
        for c in byassembly[r['assembly_accession']]:
            rec=sources[c['contig']];fs=ordered_cds(rec);ori=c['center_strand']
            expected=neighborhood_indices(len(fs),c['center_cds'],ori,manifest['flank_CDS'],rec.annotations.get('topology','').lower()=='circular')
            assert len(c['genes'])==2*manifest['flank_CDS']+1
            assert manifest['min_aa']<=len(central[c['case_id']])<=manifest['max_aa']
            assert hashlib.sha256(central[c['case_id']].encode()).hexdigest()==c['center_sha256']
            for g,(rel,i,wrapped) in zip(c['genes'],expected):
                assert i is not None and g['cds_index']==i and g['relative_cds']==rel and g['wrapped_origin']==wrapped
                key=(c['assembly_accession'],c['contig'],i);assert key not in used;used.add(key)
                raw=decode_feature(fs[i],i,set())
                for field in ['start0','end0','strand','gene','product','sequence_status','sequence_sha256']:assert g[field]==raw[field],(key,field)
                if g['sequence_status']=='valid':assert fasta[g['sequence_id']]==raw['sequence']
                else:assert g['sequence_id']==''
                assert g['relative_strand']==g['strand']*ori
                seqstatuses[g['sequence_status']]+=1;positions+=1
            wrapcases+=c['wrapped_origin']
    assert len(central)==len(windows)==len({c['center_sha256'] for c in windows})
    assert not any(c['negative_label'] is not None or c['experimental_negative'] for c in windows)
    for a in arch.values():
        assert a['copies']==len(a['cases'])
        for p in a['positions']:assert p['valid']+p['excluded']+p['missing']==a['copies']
    validation=dict(status='passed',source_assemblies_verified=len(byassembly),source_CDS_positions_verified=positions,source_sequence_status_counts=dict(seqstatuses),cases_wrapping_circular_origin=wrapcases,exact_distinct_centers=len(central),no_selected_CDS_window_overlap=True,source_coordinates_and_translations_verified=True,architecture_diversity_counts_verified=True,no_experimental_negative_labels=True,unit_tests='Run test_extraction.py separately; this report verifies the actual source records.',browser_preview='Not performed by this report generator.')
    (out/'validation.json').write_text(json.dumps(validation,indent=2),encoding='utf8')

    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    fig,axes=plt.subplots(1,3,figsize=(14,4.6))
    lens=[c['center_aa_length'] for c in windows]
    bins=list(np.arange(manifest['min_aa'],manifest['max_aa'],20))+[manifest['max_aa']+1]
    axes[0].hist(lens,bins=bins,color='#2c7b91',edgecolor='white')
    axes[0].set_xlabel('Center protein length (aa)');axes[0].set_ylabel('Sampled centers');axes[0].set_title('A | Length-stratified pilot',loc='left')
    evidence=manifest['center_evidence_counts'];values=[evidence.get('source_function_annotated',0),evidence.get('source_function_uncharacterized',0)]
    axes[1].bar(['Source function\nannotated','Source function\nuncharacterized'],values,color=['#2c7b91','#a8b4be'])
    for j,n in enumerate(values):axes[1].text(j,n+5,str(n),ha='center')
    axes[1].set_ylim(0,max(values)*1.16);axes[1].set_title('B | Annotation evidence',loc='left')
    available=Counter(c['valid_neighbor_sequences'] for c in windows)
    axes[2].bar(sorted(available),[available[k] for k in sorted(available)],color='#cc963d')
    for k,n in available.items():axes[2].text(k,n+5,str(n),ha='center',fontsize=9)
    axes[2].set_xticks(sorted(available));axes[2].set_xlabel('Available neighbor proteins / '+str(2*manifest['flank_CDS']))
    axes[2].set_ylim(0,max(available.values())*1.16);axes[2].set_title('C | CDS positions are never dropped',loc='left')
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    fig.suptitle(f"GTDB R232 pilot | {len(byassembly)} genomes | {len(windows)} unlabeled genomic backgrounds",x=.05,ha='left',fontweight='bold',fontsize=15)
    fig.tight_layout(rect=[0,.09,1,.91],w_pad=3)
    fig.text(.05,.015,'Taxonomic diversity and length stratification are sampling choices, not natural prevalence estimates.\nAll windows have '+str(2*manifest['flank_CDS']+1)+' observed CDS positions. Pseudogenes or unavailable translations remain at their original positions.',fontsize=9,color='#5d7485')
    for ext in ('png','svg'):fig.savefig(out/('sampling_overview.'+ext),dpi=170,bbox_inches='tight')
    plt.close(fig)

    # Source values are JSON-serialized into an external script; local file pages need no fetch().
    (out/'atlas_data.js').write_text('window.backgroundCases='+json.dumps(windows,ensure_ascii=True,separators=(',',':'))+';\nwindow.backgroundArchitectures='+json.dumps(arch,ensure_ascii=True,separators=(',',':'))+';',encoding='utf8')
    (out/'atlas.js').write_text(ATLAS_JS,encoding='utf8')
    rows=''
    for r in downloads:
        rows+='<tr><td>'+html.escape(r['assembly_accession'])+'</td><td>'+html.escape(r['genus'])+'</td><td>'+html.escape(r['taxonomic_class'])+'</td><td>'+str(len(byassembly[r['assembly_accession']]))+'</td></tr>'
    excluded=positions-seqstatuses['valid']
    page=HTML.replace('@@GENOME_ROWS@@',rows).replace('@@CASES@@',str(len(windows))).replace('@@GENOMES@@',str(len(byassembly))).replace('@@FASTA@@',str(len(fasta))).replace('@@EXCLUDED@@',str(excluded)).replace('@@UNCERTAIN@@',str(evidence.get('source_function_uncharacterized',0))).replace('@@ANNOTATED@@',str(evidence.get('source_function_annotated',0))).replace('@@PARTIAL@@',str(manifest['cases_with_unavailable_neighbor_sequences'])).replace('@@REVIEW@@',str(manifest['context_review_cases'])).replace('@@PHYLA@@',str(selection['phyla'])).replace('@@ARCH@@',str(len(arch)))
    page=page.replace('±7 CDS','±'+str(manifest['flank_CDS'])+' CDS').replace('40–300 aa',str(manifest['min_aa'])+'–'+str(manifest['max_aa'])+' aa').replace('20 个实际使用的组装',str(len(byassembly))+' 个实际使用的组装').replace('每个基因组最多取 20 个中心','每个基因组最多取 '+str(manifest['per_genome_cap'])+' 个中心').replace('1,047 条 EEX',str(manifest['known_EEX_exact_sequences_checked'])+' 条 EEX').replace('完整的 15 个 CDS','完整的 '+str(2*manifest['flank_CDS']+1)+' 个 CDS').replace('400 个中心蛋白',str(len(windows))+' 个中心蛋白')
    # Bundle a code snapshot so a report remains portable outside this workspace.
    snapshot=out/'pipeline';snapshot.mkdir(exist_ok=True)
    source=Path(__file__).resolve().parent
    for item in sorted(source.glob('*.py')) + [source/'README.md',source/'requirements.txt',source/'slurm_analyze.sbatch']:
        (snapshot/item.name).write_bytes(item.read_bytes())
    (out/'index.html').write_text(page,encoding='utf8')
    for url in re.findall(r'(?:href|src)="([^"#]+)"',page):
        if not re.match(r'https?://',url):assert (out/url).exists(),url
    print(json.dumps(validation,indent=2))

ATLAS_JS=r"""'use strict';
const cases=window.backgroundCases, architectures=window.backgroundArchitectures;
const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
let page=0;
const genus=document.getElementById('genus'), evidence=document.getElementById('evidence'), review=document.getElementById('review'), query=document.getElementById('query');
[...new Set(cases.map(c=>c.genus))].sort().forEach(g=>{const o=document.createElement('option');o.value=g;o.textContent=g;genus.appendChild(o);});
function draw(c,a){
  let s='<svg viewBox="0 0 '+(40+c.genes.length*98)+' 210" class="track" role="img" aria-label="'+esc(c.case_id+' observed CDS neighborhood')+'">';
  c.genes.forEach((g,i)=>{
    const x=20+i*98,cx=x+42,y=76,p=a.positions[i],invalid=g.sequence_status!=='valid';
    const col=g.relative_cds===0?'#247d91':g.source_unknown?'#bcc6ce':'#718fab';
    const flag=g.conjugation_annotation_flag||g.exclusion_annotation_flag||g.exact_project_EEX;
    const lab=g.gene||(g.source_unknown?'unknown':g.product);
    const short=lab.length>14?lab.slice(0,13)+'…':lab;
    const pts=g.relative_strand<0?[[x+82,y-12],[x+14,y-12],[x,y],[x+14,y+12],[x+82,y+12]]:[[x,y-12],[x+68,y-12],[x+82,y],[x+68,y+12],[x,y+12]];
    s+='<g><title>'+esc(g.product+' | gene='+g.gene+' | '+g.locus_tag+' | '+g.protein_accession+' | '+g.sequence_status+' | '+g.location)+'</title>';
    s+='<text x="'+cx+'" y="35" text-anchor="middle" font-size="11" fill="#213b4e">'+esc(short)+'</text>';
    s+='<text x="'+cx+'" y="52" text-anchor="middle" font-size="9" fill="#607886">'+g.aa_length+' aa'+(flag?' ‡':'')+'</text>';
    s+='<polygon points="'+pts.map(p=>p.join(',')).join(' ')+'" fill="'+col+'" stroke="'+(flag?'#cb8f2d':invalid?'#ad5961':'white')+'" stroke-width="'+(flag||invalid?'2':'1')+'"'+(invalid?' stroke-dasharray="4 2"':'')+'/>';
    s+='<text x="'+cx+'" y="117" text-anchor="middle" font-size="11">D/V '+p.distinct_aa+'/'+p.valid+'</text><text x="'+cx+'" y="138" text-anchor="middle" font-size="10" fill="#607886">M/X '+p.missing+'/'+p.excluded+'</text>';
    s+='<text x="'+cx+'" y="159" text-anchor="middle" font-size="10" fill="#607886">'+(g.relative_cds>0?'+':'')+g.relative_cds+(g.wrapped_origin?' ↻':'')+'</text>';
    s+='<text x="'+cx+'" y="181" text-anchor="middle" font-size="8" fill="#607886">'+esc(g.locus_tag.length>16?g.locus_tag.slice(-16):g.locus_tag)+'</text></g>';
  });return s+'</svg>';
}
function card(c){
  const a=architectures[c.architecture_id],freq=Object.entries(a.center_frequency_distribution).map(([n,k])=>k+' × '+n).join('；');
  return '<article class="card"><h3>'+esc(c.case_id)+' · '+esc(c.genus)+' <small>'+esc(c.architecture_id)+'</small></h3><p><b>中心 '+c.center_aa_length+' aa · '+esc(c.center_gene||c.center_product||'unknown')+'</b></p><p>源注释：'+(c.center_evidence==='source_function_annotated'?'有功能描述':'功能未明')+'；窗口文字标记复核：'+(c.context_annotation_review?'有':'未见')+'；负标签：未赋予。</p><p>架构 copies '+a.copies+'；中心完整 aa '+a.distinct_center_aa+' 种；频率（序列数 × 出现次数）：'+freq+'</p><div class="scroll">'+draw(c,a)+'</div><p class="muted">'+esc(c.assembly_accession)+' · '+esc(c.contig)+' · CDS index '+c.center_cds+'（0-based） · '+(c.contig_circular?'源记录标注为环状':'线性记录')+' · '+esc(c.replicon)+'</p><details><summary>完整源注释与分类</summary><p>'+esc(c.gtdb_taxonomy)+'</p><table><thead><tr><th>相对 CDS</th><th>Gene / locus</th><th>原始 product</th><th>序列状态</th></tr></thead><tbody>'+c.genes.map(g=>'<tr><td>'+g.relative_cds+'</td><td>'+esc(g.gene||g.locus_tag)+'</td><td>'+esc(g.product||'未提供')+'</td><td>'+esc(g.sequence_status)+'</td></tr>').join('')+'</tbody></table></details></article>';
}
function filtered(){const q=query.value.trim().toLowerCase();return cases.filter(c=>(!genus.value||c.genus===genus.value)&&(!evidence.value||c.center_evidence===evidence.value)&&(!review.value||String(c.context_annotation_review)===review.value)&&(!q||[c.case_id,c.contig,c.assembly_accession,c.architecture_id,...c.genes.map(g=>g.gene+' '+g.product+' '+g.locus_tag)].join(' ').toLowerCase().includes(q)));}
function render(){const a=filtered(),n=Math.max(1,Math.ceil(a.length/10));page=Math.min(page,n-1);document.getElementById('cards').innerHTML=a.slice(page*10,(page+1)*10).map(card).join('')||'<p>没有匹配的抽样记录。</p>';document.getElementById('count').textContent=a.length+' 个窗口 · '+(page+1)+' / '+n+' 页';document.getElementById('prev').disabled=page===0;document.getElementById('next').disabled=page>=n-1;}
[genus,evidence,review].forEach(e=>e.addEventListener('change',()=>{page=0;render();}));query.addEventListener('input',()=>{page=0;render();});
document.getElementById('prev').addEventListener('click',()=>{page--;render();document.getElementById('atlas').scrollIntoView();});document.getElementById('next').addEventListener('click',()=>{page++;render();document.getElementById('atlas').scrollIntoView();});render();
"""

HTML=r'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>GTDB 短蛋白背景窗口 · Pilot</title><style>
*{box-sizing:border-box}body{margin:0;background:#f2f6f8;color:#20394b;font:16px/1.7 "Segoe UI","Microsoft YaHei",sans-serif}main{max-width:1420px;margin:auto;padding:38px 30px 70px}header{border-top:5px solid #247d91;padding-top:25px}h1{font-size:35px;line-height:1.3;margin:12px 0}h2{font-size:25px;margin-top:36px}h3{font-size:18px}small,.muted{color:#627c8e}a{color:#1c7088}p{max-width:1150px}.eyebrow{letter-spacing:2px;font-size:12px;color:#247d91;font-weight:700}.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:15px;margin:24px 0}.metric,.card,figure{background:white;border:1px solid #d8e3e9;border-radius:9px;padding:18px}.metric b{display:block;font-size:30px}.metric span{color:#627c8e;font-size:14px}.note{background:#fff7e8;border-left:4px solid #c49443;padding:16px 20px}figure{margin:22px 0}figure img{display:block;width:100%;height:auto}figcaption{font-size:14px;color:#627c8e;padding-top:10px}table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:8px 12px;border-bottom:1px solid #dde6ec;text-align:left}th{background:#e9f0f4}.scroll{overflow:auto}.track{width:100%;min-width:1250px;height:auto;display:block}.card{margin:16px 0}.card p{font-size:13px;margin:5px 0;max-width:none}.card small{font-size:12px;font-weight:normal}.toolbar{display:flex;gap:10px;flex-wrap:wrap;align-items:center;background:#e4eef3;padding:14px;border-radius:8px}select,input,button{font:inherit;border:1px solid #adc1ce;padding:7px 10px;border-radius:5px;background:white;color:#20394b}button{cursor:pointer}button:disabled{opacity:.4;cursor:default}summary{cursor:pointer;font-weight:600}details{margin:15px 0}.files{columns:2}.files a{display:block;margin:7px 0;break-inside:avoid}.foot{border-top:1px solid #d8e3e9;margin-top:30px;padding-top:18px;font-size:13px;color:#627c8e}@media(max-width:760px){main{padding:20px 14px}.metrics{grid-template-columns:repeat(2,1fr)}h1{font-size:27px}.files{columns:1}}
</style></head><body><main><header><div class="eyebrow">GTDB R232 · GENERAL GENOMIC BACKGROUNDS · 2026-10-26</div><h1>走出接合操纵子：短蛋白及 ±7 CDS 背景样本</h1><p>已将抽样与邻域提取写成可重复运行的代码，并完成 @@GENOMES@@ 个不同属的基因组试跑。这里收集的是<strong>未标记的背景样本</strong>；序列短、没有 exclusion 注释或周围未见接合关键词，都不等于经过验证的 EEX 阴性。</p></header>
<div class="metrics"><div class="metric"><b>@@GENOMES@@</b><span>GTDB 分类的完整 RefSeq 组装</span></div><div class="metric"><b>@@CASES@@</b><span>完整短蛋白中心 · 40–300 aa</span></div><div class="metric"><b>@@FASTA@@</b><span>全部窗口内去重后的有效蛋白序列</span></div><div class="metric"><b>±7 CDS</b><span>保留顺序、方向、坐标和无效 CDS</span></div></div>
<nav><a href="#design">抽样设计</a> · <a href="#atlas">全部窗口图集</a> · <a href="#files">代码与数据</a></nav>
<section id="design"><h2>这批数据回答“能否构建外部背景”，尚未回答识别效果</h2>
<p>固定 GTDB R232 的分类与组装编号，筛选物种代表、完整 RefSeq 组装、CheckM2 完整度 ≥95、污染度 ≤5、非 MAG/SAG 类别。按分类阶元轮流抽样并限制每属一份，覆盖 @@PHYLA@@ 个门；再从相应 NCBI 注释中取短蛋白及邻域。这些蛋白是 <b>NCBI 的基因调用与注释</b>，并非 GTDB 批量 Prodigal 蛋白集。</p>
<p>每个基因组最多取 20 个中心，按长度区间分层、固定种子和 SHA256 顺序选择。全局排除完全重复的中心序列；同一 contig 内，选定窗口不共享 CDS。中心已与项目中 1,047 条 EEX 相关完整序列做精确相同检查；近同源筛查尚未运行。</p>
<figure><img src="sampling_overview.png" alt="短蛋白长度、源注释证据和邻居序列可用性"><figcaption>这是分层抽样后的分布，不能解释为自然界的出现频率。@@ANNOTATED@@ 个中心有源功能描述，@@UNCERTAIN@@ 个功能未明；二者都没有被赋予实验阴性标签。</figcaption></figure>
<div class="note"><b>保持比较公平：</b>这批基因组分类跨度大，且所有窗口都有完整的 15 个 CDS 位置。现有 EEX 数据的许多窗口被截断。后续需要匹配分类、长度、膜特征和窗口可见性，防止模型仅凭分类差异或裁剪方式区分两组。</div>
<p>@@PARTIAL@@ 个窗口含不可用的邻居蛋白序列，共 @@EXCLUDED@@ 个位置被单独记为序列排除。它们仍占据原来的 CDS 位置，未被删除后重新编号。只有源记录明确标注环状且 CDS 足够多时才跨原点取邻居；不跨 contig 拼接。</p>
<p>@@REVIEW@@ 个抽样窗口含接合 / exclusion 相关文字标记，单独保留为复核层；这只是对源注释文字的检查，未做新的功能推断。其余窗口的“未见文字标记”也不能证明相关功能不存在。中心的 SignalP、膜拓扑及 ESM 尚未在本次抽样阶段计算。</p>
<details><summary>20 个实际使用的组装</summary><div class="scroll"><table><thead><tr><th>组装编号</th><th>GTDB 属</th><th>GTDB 类</th><th>窗口数</th></tr></thead><tbody>@@GENOME_ROWS@@</tbody></table></div></details>
</section><section id="atlas"><h2>全部 @@CASES@@ 个观察窗口</h2>
<p>绿色中心；蓝灰色为有源功能描述的邻居；浅灰色为源功能未明。虚线边框表示蛋白序列不可用，金色边框与 ‡ 表示源注释文字需要复核。鼠标停在基因上可看完整 product、坐标和序列状态；箭头等宽表示 CDS 顺序，不代表物理长度。</p>
<p class="muted">D/V = 对应架构位置的不同完整 aa 条数 / 有效出现次数；M/X = 不可见位置 / 序列排除次数。Copy = 一次抽样中心出现，不等于独立生物学样本。原始未知描述与 locus 标签保留。当前有 @@ARCH@@ 种源注释架构；本批尚未建立逐中心已知功能参考匹配。</p>
<div class="toolbar"><label>属 <select id="genus"><option value="">全部</option></select></label><label>中心证据 <select id="evidence"><option value="">全部</option><option value="source_function_annotated">有功能描述</option><option value="source_function_uncharacterized">功能未明</option></select></label><label>文字复核 <select id="review"><option value="">全部</option><option value="true">有标记</option><option value="false">未见标记</option></select></label><input id="query" placeholder="窗口 / 基因 / product / 记录编号"><button id="prev">上一页</button><button id="next">下一页</button><span id="count" aria-live="polite"></span></div>
<div id="cards"></div><noscript>图集需要 JavaScript；上方静态图与下方数据下载仍可用。</noscript></section>
<section id="files"><h2>可重复运行的代码与数据</h2><p>入口脚本支持设置基因组数、每基因组中心数、长度范围及邻域宽度；下载缓存保留固定版本、原始注释、SHA256 和 GTDB 官方 MD5 校验。样本提取可从已下载的本地 GBFF 独立重跑。</p>
<div class="files"><a href="pipeline/run_pipeline.py">运行入口 · Python</a><a href="pipeline/README.md">使用说明与数据口径</a><a href="pipeline/extract_windows.py">邻域提取代码 · Python</a><a href="pipeline/test_extraction.py">坐标与注释边界测试</a>
<a href="centers.tsv">中心序列、分类与证据表</a><a href="center_proteins.faa">400 个中心蛋白 · FASTA</a><a href="neighbors.tsv">全部邻居 CDS、原始标签及坐标</a><a href="window_proteins_unique.faa">窗口内去重的完整蛋白 · FASTA</a>
<a href="windows.jsonl">完整窗口结构与序列 ID 映射</a><a href="architectures.json">架构 copies 与中心序列频率</a><a href="architecture_position_diversity.tsv">每个架构、每个位置的 D/V 与 M/X</a><a href="center_sampling_audit.tsv">中心筛选与抽样审计</a>
<a href="selected_assemblies.tsv">GTDB 组装抽样清单</a><a href="download_manifest.json">实际下载来源与文件校验</a><a href="genome_selection_manifest.json">基因组筛选规则</a><a href="extraction_manifest.json">提取参数与限制</a><a href="validation.json">源坐标 / 序列验证</a><a href="sampling_overview.svg">抽样概览 · SVG</a>
</div>
<p>来源：<a href="https://data.gtdb.ecogenomic.org/releases/release232/232.0/">GTDB 固定版本目录</a> · <a href="https://gtdb.ecogenomic.org/faq">GTDB 的 NCBI 下载说明</a> · <a href="https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/">NCBI 组装包与注释文件说明</a>。</p></section>
<p class="foot">本阶段完成数据采样与邻域提取。所有背景的 negative_label 均为空；没有更新项目已有的 EEX 注释、候选分类或 novelty shortlist。</p></main><script src="atlas_data.js"></script><script src="atlas.js"></script></body></html>'''

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args();build(a.out)
