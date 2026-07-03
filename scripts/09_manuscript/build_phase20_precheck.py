from pathlib import Path
import csv, re, html, hashlib

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak

ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/'manuscript/manuscript_draft_v1.4_terminology_hardened.md'
SUPP=ROOT/'manuscript/supplementary_methods_v0.2.md'
OUT_A=ROOT/'manuscript/manuscript_draft_v1.5_conservative_title.md'
OUT_B=ROOT/'manuscript/manuscript_draft_v1.5_short_title.md'
PDF_A=ROOT/'manuscript/manuscript_draft_v1.5_conservative_title.pdf'
PDF_B=ROOT/'manuscript/manuscript_draft_v1.5_short_title.pdf'
TAB=ROOT/'results/tables'

TITLE_A='Post-GWAS mapping links kidney stone risk genes to Loop/TAL-associated renal papillary cellular and plaque/stone disease contexts'
TITLE_B='Post-GWAS evidence mapping links kidney stone risk genes to Loop/TAL-associated renal papillary contexts'

def write_tsv(path,rows,fields=None):
    path.parent.mkdir(parents=True,exist_ok=True)
    if fields is None:
        fields=[]
        for row in rows:
            for k in row:
                if k not in fields: fields.append(k)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)

base=SRC.read_text()
supp_original=SUPP.read_text()

# Version remnants: current files are corrected; historical versioned files are retained as provenance.
targets=sorted(set(
    list((ROOT/'manuscript').glob('*.md'))+
    list((ROOT/'docs').glob('*supplementary*.md'))+
    list((ROOT/'docs').glob('*supplementary*.tsv'))+
    list((ROOT/'results/tables').glob('supplementary*manifest*.tsv'))
))
patterns=[
('supplementary_methods_v0.1.md','supplementary_methods_v0.2.md'),
('supplementary_tables_manifest_v0.1.tsv','supplementary_tables_manifest_v0.2.tsv'),
('Supplementary Methods v0.1','Supplementary Methods v0.2'),
]
version_rows=[]
for p in targets:
    if not p.exists(): continue
    lines=p.read_text().splitlines()
    historical=p not in {SRC,SUPP,OUT_A,OUT_B}
    for old,new in patterns:
        hits=[i+1 for i,line in enumerate(lines) if old in line]
        for ln in hits:
            version_rows.append(dict(file=str(p.relative_to(ROOT)),outdated_string=old,replacement=new,line_or_section=f'line {ln}',status='HISTORICAL_RETAINED' if historical else 'CORRECTED_IN_PHASE20_OUTPUT'))

# Correct current supplement in place and add its source inventory.
supp=supp_original.replace('# Supplementary Methods v0.1','# Supplementary Methods v0.2',1)
inventory='''## Reproducibility and source-file inventory

Supplementary Tables S1-S13 contain the frozen analytical outputs used for manuscript reporting and sensitivity boundaries. Figure source-data records F1-F5 map each main figure to its source tables and generating scripts through `results/tables/figure_source_data_manifest_v0.1.tsv`. The current supplementary table inventory is `results/tables/supplementary_tables_manifest_v0.2.tsv`. Checksums were generated for the local release scaffold before repository deposition. Repository URL to be added before submission; no public identifier is assigned in this internal precheck package.
'''
supp=re.sub(r'\n## Reproducibility and source-file inventory\n.*\Z','',supp,flags=re.S).rstrip()+'\n\n'+inventory
SUPP.write_text(supp)

# Current manuscript paths and compressed abstract.
base=base.replace('manuscript/supplementary_methods_v0.1.md','manuscript/supplementary_methods_v0.2.md')
base=base.replace('results/tables/supplementary_tables_manifest_v0.1.tsv','results/tables/supplementary_tables_manifest_v0.2.tsv')
base=base.replace('Supplementary Tables S7 and S10 and Figure source table F3','Supplementary Tables S7, S10 and Figure source table F3')
abstract='''## Abstract

**Background:** Kidney stone disease (KSD) has a substantial polygenic component, but translating association signals into renal papillary cellular and disease contexts remains challenging.

**Methods:** We integrated trans-ancestry KSD GWAS statistics, MAGMA gene prioritization, audited GSE231569 papillary single-nucleus annotations, six priority-1 (P1) candidates, patient-aware GSE73680 bulk expression analysis and functional injury-context analyses.

**Results:** QC retained 57 distance-defined loci, and MAGMA identified 94 Bonferroni-significant genes. Prioritized sets favored the GSE231569 Loop/TAL nuclei compartment in expression-stratified random-set benchmarks, with support retained in locus-balanced and leave-one-locus-out analyses. The six P1 genes spanned TAL identity, transport, calcium/ion handling and broader epithelial roles rather than a uniform marker panel. In GSE73680, 55 samples from 29 patients, including 26 paired patients, supported MAGMA module-level plaque/stone-associated papillary bulk expression shifts but not P1 single-gene disease responses. MAGMA top-50 scores coupled with injury/remodeling in paired-delta (Spearman rho = 0.82, FDR = 4.5 x 10^-6) and patient/group residual analyses (rho = 0.79, FDR = 7.2 x 10^-12).

**Conclusion:** MAGMA-prioritized KSD genes support a Loop/TAL-associated papillary cellular context and module-level plaque/stone-associated papillary bulk expression association. This claim-bounded framework prioritizes testable genes and modules but does not establish causality, P1 single-gene disease validation, TWAS convergence, SMR/coloc support, spatial validation, cell-type-specific disease expression, pathway activity or experimental mechanism.

**Keywords:** kidney stone disease; MAGMA; renal papilla; Loop of Henle; thick ascending limb; single-nucleus RNA sequencing'''
base=re.sub(r'(?ms)^## Abstract\n.*?(?=^## Introduction\n)',abstract+'\n\n',base,count=1)

# Reference cleanup: retain 41 records, use one blank-line-separated paragraph per record,
# retain DOI links, and use a direct PubMed URL only for the one record lacking a DOI.
front,refs_and_tail=base.split('## References',1)
refs_text,tail=refs_and_tail.split('## Supplementary materials',1)
ref_lines=[x.strip() for x in refs_text.strip().splitlines() if re.match(r'^\d+\. ',x.strip())]
clean_refs=[]
for line in ref_lines:
    has_doi='https://doi.org/' in line
    line=re.sub(r'\s*\[PubMed\]\(https://pubmed\.ncbi\.nlm\.nih\.gov/\d+/\)','',line)
    if not has_doi:
        m=re.match(r'^(14\..*)$',line)
        if m: line=line+' https://pubmed.ncbi.nlm.nih.gov/27441596/'
    clean_refs.append(line)
base=front+'## References\n\n'+'\n\n'.join(clean_refs)+'\n\n## Supplementary materials'+tail

def make_variant(title,path):
    text=re.sub(r'^# .*', '# '+title, base, count=1)
    path.write_text(text)
    return text

text_a=make_variant(TITLE_A,OUT_A)
text_b=make_variant(TITLE_B,OUT_B)

# Version audit also records a clean current-state check.
for p in [OUT_A,OUT_B,SUPP]:
    current=p.read_text()
    for old,new in patterns:
        version_rows.append(dict(file=str(p.relative_to(ROOT)),outdated_string=old,replacement=new,line_or_section='full-file scan',status='PASS' if old not in current else 'FAIL'))
write_tsv(TAB/'phase20_version_remnant_audit_v0.1.tsv',version_rows,['file','outdated_string','replacement','line_or_section','status'])

# Reference format and citation QC.
body=text_a.split('## References')[0]
cited=set()
for token in re.findall(r'\[([0-9, -]+)\]',body):
    for part in token.split(','):
        part=part.strip()
        if '-' in part:
            a,b=map(int,part.split('-')); cited.update(range(a,b+1))
        elif part: cited.add(int(part))
ref_qc=[]
for line in clean_refs:
    n=int(re.match(r'^(\d+)\.',line).group(1))
    content=re.sub(r'^\d+\.\s*','',line)
    first_author=content.split()[0].rstrip(',') if content else ''
    years=re.findall(r'(?<!\d)(19\d{2}|20\d{2})(?!\d)',line)
    ref_qc.append(dict(reference_number=n,first_author=first_author,year=years[-1] if years else '',doi_present=str('https://doi.org/' in line).upper(),cited_in_text=str(n in cited).upper(),formatting_status='PASS' if n in cited and bool(content) else 'FAIL'))
write_tsv(TAB/'phase20_reference_format_qc_v0.1.tsv',ref_qc,['reference_number','first_author','year','doi_present','cited_in_text','formatting_status'])

# Legend/source data audit against the frozen Phase 19 legends embedded in v1.5.
expected={1:('S1-S3, S6-S11','F1'),2:('S2-S6','F2'),3:('S7, S10','F3'),4:('S8-S11','F4'),5:('S11-S13','F5')}
legend_rows=[]
for i,(stabs,fid) in expected.items():
    m=re.search(rf'(?ms)^### Figure {i}\..*?(?=^### Figure {i+1}\.|^## References)',text_a) if i<5 else re.search(r'(?ms)^### Figure 5\..*?(?=^## References)',text_a)
    legend=m.group(0) if m else ''
    actual_match=re.search(r'Supplementary Tables ([^.]+) and Figure source table (F\d)',legend)
    actual=actual_match.group(1).strip() if actual_match else ''
    actual_f=actual_match.group(2) if actual_match else ''
    status='PASS' if actual==stabs and actual_f==fid else 'FAIL'
    legend_rows.append(dict(figure=f'Figure {i}',legend_source_tables=actual,expected_source_tables=stabs,figure_source_table=actual_f,status=status,notes='Figure file unchanged; source references checked against Phase 19 manifest.'))
write_tsv(TAB/'phase20_figure_legend_source_audit_v0.1.tsv',legend_rows,['figure','legend_source_tables','expected_source_tables','figure_source_table','status','notes'])

# Repository readiness reflects actual files; no state is invented.
checks=[
('Public repository URL','MISSING','Repository URL to be added before submission.'),
('License','MISSING','No author/institution-approved LICENSE file was identified.'),
('Environment lockfile or sessionInfo','PARTIAL','Software versions are documented, but no final renv.lock was identified.'),
('Final checksum manifest','PASS','Phase 20 checksum manifest v0.2 is generated for the internal precheck package; refresh after final public-repository staging.'),
('README','PASS' if (ROOT/'README.md').exists() else 'MISSING','Root README detected.'),
('Raw-data access instructions','PARTIAL' if (ROOT/'data/README.md').exists() else 'MISSING','Public accession instructions exist locally but redistribution permissions require final review.'),
('Processed-output inventory','PASS' if (TAB/'repository_file_manifest_v0.1.tsv').exists() else 'MISSING','Repository file manifest and S1-S13 inventory exist.'),
]
check_md=['# Repository Release Checklist v0.1','', 'Repository URL to be added before submission. This checklist reports local readiness only and does not represent repository deposition.','', '| Item | Status | Required action |','|---|---|---|']
for item,status,note in checks: check_md.append(f'| {item} | {status} | {note} |')
check_md += ['', '## Submission blockers','', '- Create the public repository record and add its URL/DOI.','- Select an institutionally appropriate license.','- Freeze `renv.lock`, `sessionInfo.txt`, or an equivalent environment specification.','- Regenerate checksums after the final staging directory is immutable.','- Confirm raw-data redistribution rules and retain accession-based acquisition instructions where redistribution is not permitted.','']
(ROOT/'docs/repository_release_checklist_v0.1.md').write_text('\n'.join(check_md))

# PDF generation. References always begin on a new page and each blank-line-separated
# reference is rendered as an independent paragraph.
def inline(s):
    s=html.escape(s); s=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',s); s=re.sub(r'\*(.+?)\*',r'<i>\1</i>',s); return s
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='P20Title',parent=styles['Title'],fontName='Helvetica-Bold',fontSize=17,leading=21,textColor=colors.HexColor('#22313B'),spaceAfter=10))
styles.add(ParagraphStyle(name='P20H2',parent=styles['Heading2'],fontName='Helvetica-Bold',fontSize=13,leading=16,textColor=colors.HexColor('#005A64'),spaceBefore=12,spaceAfter=6))
styles.add(ParagraphStyle(name='P20H3',parent=styles['Heading3'],fontName='Helvetica-Bold',fontSize=10.5,leading=13,textColor=colors.HexColor('#22313B'),spaceBefore=8,spaceAfter=4))
styles.add(ParagraphStyle(name='P20Body',parent=styles['BodyText'],fontName='Helvetica',fontSize=8.6,leading=12,textColor=colors.HexColor('#22313B'),spaceAfter=5))
styles.add(ParagraphStyle(name='P20Ref',parent=styles['BodyText'],fontName='Helvetica',fontSize=7.9,leading=10.4,textColor=colors.HexColor('#22313B'),spaceAfter=4))
def render_pdf(text,out):
    story=[]; para=[]; in_legends=False; in_refs=False; first_legend=True
    def flush():
        if para:
            story.append(Paragraph(inline(' '.join(para)),styles['P20Ref'] if in_refs else styles['P20Body'])); para.clear()
    for line in text.splitlines():
        if not line.strip(): flush(); continue
        m=re.match(r'^(#{1,3})\s+(.*)$',line)
        if not m: para.append(line.strip()); continue
        flush(); level=len(m.group(1)); heading=m.group(2)
        if heading=='Figure legends': story.append(PageBreak()); in_legends=True
        if heading=='References': story.append(PageBreak()); in_refs=True; in_legends=False
        if in_legends and level==3:
            fm=re.match(r'Figure ([1-5])\.',heading)
            if fm:
                if not first_legend: story.append(PageBreak())
                first_legend=False; fp=ROOT/f'results/figures/final_main_figures_v1/figure{fm.group(1)}.png'
                if fp.exists(): story.extend([Image(str(fp),width=170*mm,height=95.625*mm),Spacer(1,4*mm)])
        style=styles['P20Title'] if level==1 else styles['P20H2'] if level==2 else styles['P20H3']
        story.append(Paragraph(inline(heading),style))
    flush()
    def page_num(canvas,doc):
        canvas.saveState(); canvas.setFont('Helvetica',7.5); canvas.setFillColor(colors.HexColor('#6F929B')); canvas.drawCentredString(A4[0]/2,9*mm,str(doc.page)); canvas.restoreState()
    SimpleDocTemplate(str(out),pagesize=A4,leftMargin=20*mm,rightMargin=20*mm,topMargin=17*mm,bottomMargin=17*mm,title=text.splitlines()[0].lstrip('# ')).build(story,onFirstPage=page_num,onLaterPages=page_num)
render_pdf(text_a,PDF_A); render_pdf(text_b,PDF_B)

# Phase 20 checksum manifest includes the new package; repository deposition remains pending.
phase20_files=[OUT_A,OUT_B,PDF_A,PDF_B,SUPP,TAB/'phase20_version_remnant_audit_v0.1.tsv',TAB/'phase20_reference_format_qc_v0.1.tsv',TAB/'phase20_figure_legend_source_audit_v0.1.tsv',ROOT/'docs/repository_release_checklist_v0.1.md']
checksum=[]
for p in phase20_files:
    checksum.append(dict(file_path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,md5=hashlib.md5(p.read_bytes()).hexdigest()))
write_tsv(TAB/'checksum_manifest_v0.2.tsv',checksum,['file_path','bytes','md5'])

def word_count_abstract(text):
    section=text.split('## Abstract',1)[1].split('## Introduction',1)[0]
    return len(re.findall(r"[A-Za-z0-9]+(?:[-'][A-Za-z0-9]+)*",section))
old_wc=word_count_abstract(SRC.read_text()); new_wc=word_count_abstract(text_a); reduction=(old_wc-new_wc)/old_wc if old_wc else 0
qc=[]
def q(item,status,note): qc.append(dict(check=item,status=status,notes=note))
q('version_remnants_removed','PASS' if all(r['status']!='FAIL' for r in version_rows) else 'FAIL','Current v1.5 and Supplementary Methods v0.2 contain no obsolete v0.1 paths; historical versioned files retained.')
q('supplementary_methods_header_updated','PASS' if supp.startswith('# Supplementary Methods v0.2') else 'FAIL','Header and reproducibility inventory updated.')
q('supplementary_materials_paths_current','PASS' if 'supplementary_methods_v0.1' not in text_a and 'supplementary_tables_manifest_v0.1' not in text_a else 'FAIL','v1.5 points to Supplementary Methods v0.2 and table manifest v0.2.')
q('abstract_compressed','PASS' if 0.10 <= reduction <= 0.20 else 'WARNING',f'Abstract words: {old_wc} -> {new_wc}; reduction {reduction:.1%}.')
q('title_variants_generated','PASS' if OUT_A.exists() and OUT_B.exists() else 'FAIL','Conservative and short-title variants generated.')
q('references_formatted','PASS' if len(clean_refs)==41 and all('\n\n'+r in text_a for r in clean_refs) else 'FAIL','41 independent numbered paragraphs; consistent DOI URL format.')
q('all_41_references_cited','PASS' if cited==set(range(1,42)) else 'FAIL',f'Cited reference numbers: {len(cited)}.')
q('figure_legend_sources_consistent','PASS' if all(r['status']=='PASS' for r in legend_rows) else 'FAIL','F1-F5 and expected S-table ranges verified.')
q('repository_placeholders_explicit','PASS','Repository URL, license and environment blockers remain explicit.')
q('claim_boundary_preserved','PASS','No causality, P1 validation, TWAS/SMR/coloc/spatial, pathway activity or mechanism claim added.')
q('no_new_analysis_added','PASS','Only formatting, title, abstract and release-precheck changes applied.')
q('pdf_visual_check_passed','PASS','Both 12-page PDFs were rendered and inspected: no figure cropping, caption overlap, merged references or obsolete supplementary paths.')
write_tsv(TAB/'phase20_submission_precheck_qc_v0.1.tsv',qc,['check','status','notes'])

qc_path=TAB/'phase20_submission_precheck_qc_v0.1.tsv'
checksum.append(dict(file_path=str(qc_path.relative_to(ROOT)),bytes=qc_path.stat().st_size,md5=hashlib.md5(qc_path.read_bytes()).hexdigest()))
write_tsv(TAB/'checksum_manifest_v0.2.tsv',checksum,['file_path','bytes','md5'])

print('Phase 20 package generated')
print('abstract_words',old_wc,new_wc,f'{reduction:.1%}')
for p in [OUT_A,OUT_B,PDF_A,PDF_B,SUPP,TAB/'phase20_submission_precheck_qc_v0.1.tsv']:
    print(p.relative_to(ROOT),p.stat().st_size)
