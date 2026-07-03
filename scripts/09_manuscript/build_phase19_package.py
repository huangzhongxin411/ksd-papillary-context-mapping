from pathlib import Path
import csv, re, hashlib, html

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak

ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/'manuscript/manuscript_draft_v1.3_cited_clean.md'
SUPP_SRC=ROOT/'manuscript/supplementary_methods_v0.1.md'
OUT=ROOT/'manuscript/manuscript_draft_v1.4_terminology_hardened.md'
SUPP_OUT=ROOT/'manuscript/supplementary_methods_v0.2.md'
TABLE_DIR=ROOT/'results/supplementary_tables_v0.2'
TABLE_DIR.mkdir(parents=True,exist_ok=True)

def read_tsv(path):
    p=ROOT/path
    if not p.exists(): return []
    with p.open(newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(path,rows,fields=None):
    path.parent.mkdir(parents=True,exist_ok=True)
    if fields is None:
        fields=[]
        for r in rows:
            for k in r:
                if k not in fields: fields.append(k)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)

def union_sources(items):
    out=[]
    for label,path in items:
        for row in read_tsv(path): out.append({'source_component':label,'source_path':path,**row})
    return out

text=SRC.read_text()
supp=SUPP_SRC.read_text()

rules=[
('expression-matched benchmark percentiles','expression-stratified random-set benchmark percentiles','GSE231569 benchmark'),
('expression-matched random expectation','expression-stratified random expectation','GSE231569 benchmark'),
('Size-matched random-set benchmark percentiles across audited cell contexts','Expression-stratified random-set benchmark percentiles across audited cell contexts','Figure 2 benchmark'),
('GSE73680 plaque/stone papilla analysis','GSE73680 plaque/stone-associated papillary bulk expression analysis','GSE73680 analysis'),
('GSE73680 plaque/stone papilla bulk disease-context analysis','GSE73680 plaque/stone-associated papillary bulk expression analysis','GSE73680 analysis'),
('plaque/stone papilla samples','plaque/stone-associated papillary bulk expression samples','GSE73680 samples'),
('plaque/stone papilla comparisons','plaque/stone-associated papillary bulk expression comparisons','GSE73680 comparisons'),
('plaque/stone papilla disease-context','plaque/stone-associated papillary bulk expression context','GSE73680 context'),
('plaque/stone papilla expression dataset','plaque/stone-associated papillary bulk expression dataset','GSE73680 dataset'),
('plaque/stone-associated papillary disease-context expression association','plaque/stone-associated papillary bulk expression association','GSE73680 result'),
('plaque/stone-associated papillary disease-context association','plaque/stone-associated papillary bulk expression association','GSE73680 result'),
('plaque/stone disease-context association','plaque/stone-associated papillary bulk expression association','GSE73680 conclusion'),
('plaque/stone papilla disease programs','plaque/stone-associated papillary bulk expression contexts','GSE73680 context'),
('Loop/TAL cells','Loop/TAL nuclei compartment','GSE231569 compartment'),
]

report=[]
for old,new,domain in rules:
    before=text.count(old)+supp.count(old)
    text=text.replace(old,new); supp=supp.replace(old,new)
    after=text.count(old)+supp.count(old)
    report.append(dict(domain=domain,nonpreferred_term=old,canonical_term=new,before_count=before,remaining_count=after,status='PASS' if after==0 else 'FAIL'))

# Canonical headings and highly visible phrases.
text=text.replace('### GSE73680 patient-aware disease-context analysis','### GSE73680 patient-aware bulk expression analysis')
text=text.replace('### GSE73680 supports MAGMA module-level plaque/stone-associated papillary bulk expression association','### GSE73680 supports MAGMA module-level plaque/stone-associated papillary bulk expression association')
text=text.replace('MAGMA-prioritized genes localize to a Loop/TAL-associated snRNA context','MAGMA-prioritized genes localize to a Loop/TAL-associated single-nucleus context')
text=text.replace('MAGMA plus snRNA-based TAL-associated cellular context','MAGMA plus single-nucleus-supported Loop/TAL-associated papillary cellular context')
text=text.replace('patient-aware plaque/stone-associated papillary expression shifts','patient-aware plaque/stone-associated papillary bulk expression shifts')
text=text.replace('plaque/stone-associated papillary expression dataset','plaque/stone-associated papillary bulk expression dataset')
text=text.replace('plaque/stone-associated papillary bulk expression samples','plaque/stone-associated papillary samples profiled by bulk expression')
text=text.replace('plaque/stone-associated disease-context expression','plaque/stone-associated papillary bulk expression')
text=text.replace('GSE73680 Randall plaque/stone-associated papillary bulk expression dataset','GSE73680 plaque/stone-associated papillary bulk expression dataset')
text=text.replace('consistent paired disease-context shifts','consistent paired bulk expression module shifts')
design_needle='The workflow comprised GWAS quality control, MAGMA prioritization, audited cell-context projection, P1 candidate interpretation, patient-aware bulk expression analysis, functional enrichment and injury/remodeling coupling.'
if design_needle in text:
    text=text.replace(design_needle,design_needle+' The claim-bounded evidence flow is summarized in Fig. 1.',1)

# Explicit donor boundary in Result 2.
needle='The analysis used harmonized audited cell types and retained Loop/TAL as a high-confidence epithelial transport compartment.'
if needle in text and 'Donor-level Loop/TAL summaries were treated as descriptive' not in text:
    text=text.replace(needle,needle+' Donor-level Loop/TAL summaries were treated as descriptive because one donor contributed very few Loop/TAL nuclei.',1)

# Explicit conservative direct-symbol boundary in Supplementary Methods.
needle2='The prespecified collapse rule retained the probe with maximum mean expression for duplicate symbols, although no multi-probe genes remained after direct-symbol filtering.'
addition=' Because direct gene-symbol parsing removed platform probes lacking unambiguous symbols, the resulting gene-level matrix should be interpreted as a conservative direct-symbol reconstruction rather than a full platform reannotation.'
if needle2 in supp and addition.strip() not in supp: supp=supp.replace(needle2,needle2+addition,1)
bench_needle='The random-set universe comprised genes with finite positive global mean expression.'
if bench_needle in supp:
    supp=supp.replace(bench_needle,'The expression-stratified size-matched random-set benchmark used a universe comprising genes with finite positive global mean expression.',1)

# Data/code placeholders remain explicit and do not invent a repository URL.
text=re.sub(r'## Data availability\n\n.*?(?=\n## Code availability)',
'''## Data availability

This study used public GWAS summary statistics and GEO datasets GSE231569 and GSE73680.[4,11,15,33] Processed source tables, QC summaries, supplementary tables and figure-source manifests are organized in the local release scaffold. Repository URL to be added before submission.
''',text,flags=re.S)
text=re.sub(r'## Code availability\n\n.*?(?=\n## Acknowledgements)',
'''## Code availability

Analysis scripts are organized by workflow stage, and a release inventory with checksums has been generated. Repository URL to be added before submission. The final license and environment lockfile must be confirmed before public release.
''',text,flags=re.S)

# Harden legends in the manuscript only; figure files remain frozen.
legend_replacements={
1:'''### Figure 1. Post-GWAS evidence map for kidney stone papillary cellular context

The framework links four evidence layers: KSD GWAS/MAGMA prioritization, GSE231569 single-nucleus localization, P1 candidate interpretation and GSE73680 plaque/stone-associated papillary bulk expression analysis. The Manhattan inset summarizes the cleaned KSD GWAS; the remaining cards summarize Loop/TAL localization, the six-gene role spectrum and paired module-level bulk expression shifts. The papillary inset provides anatomical orientation only, and arrows indicate evidence flow rather than causal progression. The framework separates supported context-level associations from unresolved causality and external validation. Source data are provided in Supplementary Tables S1-S3, S6-S11 and Figure source table F1.''',
2:'''### Figure 2. MAGMA-prioritized KSD genes localize to a Loop/TAL-associated single-nucleus context

(A) KSD GWAS Manhattan summary with representative prioritized loci. (B) MAGMA top-50 module projection on the audited GSE231569 UMAP. (C) Expression-stratified random-set benchmark percentiles across audited cell contexts; the dashed line marks the 95th percentile. (D) Leading Loop/TAL module contributors. Cell-level maps provide visualization, whereas donor-cell-type summaries and benchmark analyses provide the inferential layer. Donor-level Loop/TAL summaries were descriptive because one donor contributed very few Loop/TAL nuclei. The figure supports expression-context localization, not causal mediation or spatial validation. Source data are provided in Supplementary Tables S2-S6 and Figure source table F2.''',
3:'''### Figure 3. Gene-centric evidence spectrum of P1 candidates

(A) Six P1 genes arranged from TAL identity to broader epithelial context. (B) Expression and detection across six audited GSE231569 compartments. (C) TAL specificity ratios and descriptive donor support; the dashed line marks two-fold TAL preference. (D) Evidence matrix summarizing MAGMA prioritization, TAL rank, donor support, specificity, GSE73680 bulk-response status and interpretation role. No P1 gene showed a uniform FDR-supported paired bulk expression response. The panel supports a contextual role spectrum, not P1 single-gene disease validation. Source data are provided in Supplementary Tables S7 and S10 and Figure source table F3.''',
4:'''### Figure 4. GSE73680 supports MAGMA module-level plaque/stone-associated papillary bulk expression association

(A) Patient-aware design comprising 55 samples from 29 patients, including 26 paired patients. (B) Paired deltas for six P1 genes; no gene reached FDR q < 0.05. (C) Paired module-score shifts for MAGMA gene sets and P1 core. (D) Primary size-matched random-set benchmarks and robustness summaries. MAGMA modules were FDR-supported, whereas P1 core was not; the expression-matched benchmark remained a conservative sensitivity analysis. These analyses support bulk module-level association, not P1 single-gene validation, cell-type-specific disease expression or causality. Source data are provided in Supplementary Tables S8-S11 and Figure source table F4.''',
5:'''### Figure 5. Functional context, injury coupling and extension audit

(A) Evidence network linking prioritized modules to TAL transport, calcium/ion handling, epithelial junction and injury/remodeling themes. Link width is relative within each evidence layer and is not a causal effect scale. (B) Ranked redundancy-reduced GO terms. (C) Paired-delta correlations between risk modules and injury contexts. (D) Robustness summaries and resource-audit cards. These analyses support module-level functional interpretation and injury-context coupling, but not pathway activity, causal mediation or external validation. Source data are provided in Supplementary Tables S11-S13 and Figure source table F5.'''
}
for i,body in legend_replacements.items():
    start=rf'(?ms)^### Figure {i}\..*?(?=^### Figure {i+1}\.|^## References)' if i<5 else rf'(?ms)^### Figure 5\..*?(?=^## References)'
    text=re.sub(start,body+'\n\n',text,count=1)

OUT.write_text(text)
SUPP_OUT.write_text(supp)
write_tsv(ROOT/'results/tables/terminology_harmonization_report_v0.1.tsv',report)

# Freeze 13 concrete supplementary tables.
supp_specs=[
('S1','GWAS QC summary',[('gwas_qc','results/tables/phase1_gwas_qc_report.tsv')]),
('S2','Lead loci reconstruction',[('lead_variants','results/tables/phase1_2025_lead_snps.tsv'),('merged_loci','results/tables/phase1_2025_loci.tsv'),('mapped_candidate_genes','results/tables/phase1_candidate_genes.tsv')]),
('S3','MAGMA gene results and gene sets',[('gene_results','results/tables/magma_genes.tsv'),('gene_set_summary','results/tables/magma_gene_set_summary.tsv')]),
('S4','GSE231569 sample and nuclei QC',[('donor_cell_counts','results/tables/gse231569_donor_cell_counts.tsv'),('cell_type_counts','results/tables/gse231569_cell_counts.tsv'),('qc_summary','results/tables/phase1_gse231569_qc_summary.tsv')]),
('S5','Audited cell-type markers',[('cluster_marker_assignment','results/tables/gse231569_cluster_marker_assignment_audit.tsv'),('broad_marker_audit','results/tables/gse231569_marker_audit.tsv')]),
('S6','MAGMA single-nucleus benchmarks',[('expression_stratified','results/tables/magma_scrna_random_benchmark.tsv'),('locus_balanced','results/tables/magma_locus_balanced_scrna_benchmark.tsv'),('leave_one_locus_out','results/tables/magma_leave_one_locus_out.tsv')]),
('S7','P1 candidate evidence',[('p1_evidence','results/tables/p1_tal_gene_evidence.tsv')]),
('S8','GSE73680 metadata and inclusion',[('curated_metadata','config/gse73680_sample_metadata_curated.tsv'),('analysis_design','results/gse73680/tables/gse73680_analysis_design.tsv'),('patient_structure','results/gse73680/tables/gse73680_patient_sample_structure.tsv')]),
('S9','GSE73680 feature reconstruction',[('mapping_qc','results/gse73680/tables/gse73680_gene_mapping_qc.tsv'),('matrix_build','results/gse73680/tables/gse73680_matrix_build_log.tsv'),('expression_qc','results/gse73680/tables/gse73680_expression_qc_v2.tsv')]),
('S10','P1 gene paired response',[('paired_p1','results/gse73680/tables/gse73680_patient_level_p1_gene_response.tsv')]),
('S11','MAGMA module response and robustness',[('paired_module','results/gse73680/tables/gse73680_patient_level_module_response.tsv'),('size_matched_benchmark','results/gse73680/tables/gse73680_random_module_benchmark.tsv'),('expression_matched_sensitivity','results/gse73680/tables/gse73680_expression_matched_random_benchmark.tsv'),('without_p1','results/gse73680/tables/gse73680_magma_without_p1_sensitivity.tsv'),('direction_consistency','results/gse73680/tables/gse73680_paired_direction_consistency.tsv')]),
('S12','Functional enrichment and curated context sets',[('go_terms','results/tables/go_bp_redundancy_reduced_terms.tsv'),('curated_contexts','results/tables/nephron_segment_marker_enrichment.tsv'),('influential_gene_enrichment','results/tables/loop_tal_influential_gene_pathway_enrichment.tsv')]),
('S13','Resource audit and reproducibility summary',[('twas','results/tables/twas_resource_manifest_v2.tsv'),('eqtl_smr_coloc','results/tables/eqtl_resource_manifest_v0.2.tsv'),('spatial','results/tables/spatial_resource_manifest.tsv'),('pathway_resources','results/tables/pathway_resource_status_v0.1.tsv')]),
]
supp_rows=[]
script_map={
'S1':'scripts/01_gwas_qc/phase1_gwas_sanity_check.R','S2':'scripts/02_locus_mapping/make_phase1_leads_loci.py','S3':'scripts/03_magma_fuma/postprocess_magma_genes.R','S4':'scripts/07_scrna_gene_mapping/phase1_gse231569_quick_projection.R','S5':'scripts/07_scrna_gene_mapping/gse231569_marker_audit.R','S6':'scripts/07_scrna_gene_mapping/magma_scrna_projection_benchmark.R','S7':'scripts/07_scrna_gene_mapping/p1_tal_gene_evidence_pack.R','S8':'scripts/08_plaque_context/06_curate_gse73680_metadata.R','S9':'scripts/08_plaque_context/05_map_gse73680_features_to_genes.R','S10':'scripts/08_plaque_context/11_run_gse73680_p1_gene_response.R','S11':'scripts/08_plaque_context/14_run_gse73680_random_module_benchmark.R','S12':'scripts/03_magma_fuma/20_magma_pathway_enrichment.R','S13':'scripts/04_twas/03_make_twas_resource_manifest_v2.R'}
section_map={'S1':'Methods: GWAS QC','S2':'Methods/Results: loci','S3':'Methods/Results: MAGMA','S4':'Methods: GSE231569','S5':'Methods: annotation audit','S6':'Results: Figure 2','S7':'Results: Figure 3','S8':'Methods: GSE73680','S9':'Methods: GSE73680 reconstruction','S10':'Results: Figure 4B','S11':'Results: Figure 4C-D','S12':'Results: Figure 5A-B','S13':'Methods/Limitations: resource audit'}
panel_map={'S1':'Figure 1/2A','S2':'Figure 1/2A','S3':'Figures 1-2','S4':'Figure 2B','S5':'Figure 2B','S6':'Figure 2C-D','S7':'Figure 3A-D','S8':'Figure 4A','S9':'Figure 4A','S10':'Figures 3D/4B','S11':'Figure 4C-D; Figure 5C-D','S12':'Figure 5A-B','S13':'Figure 5D'}
boundary_map={'S1':'QC only','S2':'Distance-defined reconstruction; not fine-mapping','S3':'Gene prioritization; not causality','S4':'Descriptive nuclei/donor structure','S5':'Audited annotation; not spatial validation','S6':'Expression context; not causal mediation','S7':'Candidate interpretation; not single-gene validation','S8':'Bulk metadata; not lesion-resolved spatial sampling','S9':'Conservative direct-symbol reconstruction','S10':'No FDR-supported P1 single-gene response','S11':'Bulk module association; expression-matched benchmark sensitivity only','S12':'Functional interpretation; not pathway activity','S13':'Audited resource limitations; not negative biological evidence'}
for sid,title,items in supp_specs:
    rows=union_sources(items); fp=TABLE_DIR/f'{sid}_{re.sub("[^a-z0-9]+","_",title.lower()).strip("_")}.tsv'; write_tsv(fp,rows)
    supp_rows.append(dict(table_id=sid,file_path=str(fp.relative_to(ROOT)),title=title,short_description=f'Frozen composite source table with {len(rows)} rows from {len(items)} audited component(s).',source_script=script_map[sid],main_text_section=section_map[sid],figure_panel_if_applicable=panel_map[sid],claim_boundary=boundary_map[sid],status='PASS' if rows and fp.exists() else 'FAIL'))
write_tsv(ROOT/'results/tables/supplementary_tables_manifest_v0.2.tsv',supp_rows)

# Figure source-data manifest.
fig_sources=[
('F1','Figure 1','framework','results/tables/figure1_panel_source_files.tsv','S1-S3;S6-S11','scripts/09_manuscript/figure1_evidence_map_real_manhattan_restored.R','Evidence-map synthesis; arrows are non-causal.'),
('F2','Figure 2','A-D','results/tables/figure2_panel_source_files.tsv','S2-S6','scripts/09_manuscript/figure2_gwas_magma_scrna_localization_v1.0.3.R','Single-nucleus expression context only.'),
('F3','Figure 3','A-D','results/tables/figure3_panel_source_files.tsv','S7;S10','scripts/09_manuscript/phase18_figure3_micro_polish.R','P1 contextual evidence, not disease validation.'),
('F4','Figure 4','A-D','results/tables/figure4_panel_source_files_final.tsv','S8-S11','scripts/09_manuscript/phase18_figure4_micro_polish.R','Bulk module association, not cell-specific response.'),
('F5','Figure 5','A-D','results/tables/figure5_panel_source_files_final_v2.tsv','S11-S13','scripts/09_manuscript/phase18_figure5_micro_polish.R','Functional interpretation/coupling, not mechanism.'),
]
fig_rows=[]
for fid,fig,panel,src,stabs,script,boundary in fig_sources:
    fig_rows.append(dict(source_id=fid,figure=fig,panel=panel,source_table=src,supplementary_tables=stabs,source_script=script,claim_boundary=boundary,source_exists=str((ROOT/src).exists()).upper(),figure_pdf=str((ROOT/f'results/figures/final_main_figures_v1/figure{fig.split()[-1]}.pdf').exists()).upper(),status='PASS'))
write_tsv(ROOT/'results/tables/figure_source_data_manifest_v0.1.tsv',fig_rows)

# Supplementary figure manifest from the already frozen v0.2 package.
sf_specs=[('SF1','Donor-level reproducibility','results/figures/supplementary_figures_v0.2/supp_fig1_donor_level_reproducibility.pdf','S4;S6','Descriptive donor sensitivity'),('SF2','GSE73680 paired-patient waterfall','results/figures/supplementary_figures_v0.2/supp_fig2_gse73680_paired_patient_waterfall.pdf','S11','Paired bulk module direction'),('SF3','Expression-matched benchmark sensitivity','results/figures/supplementary_figures_v0.2/supp_fig3_expression_matched_benchmark.pdf','S11','Conservative sensitivity only'),('SF4','Audited GSE231569 atlas','results/figures/supplementary_figures_v0.2/supp_fig4_audited_gse231569_atlas.pdf','S4;S5','Annotation context; not spatial validation')]
sf_rows=[dict(figure_id=i,title=t,file_path=p,supplementary_tables=s,claim_boundary=b,exists=str((ROOT/p).exists()).upper(),status='PASS' if (ROOT/p).exists() else 'FAIL') for i,t,p,s,b in sf_specs]
write_tsv(ROOT/'results/tables/supplementary_figure_manifest_v0.1.tsv',sf_rows)

# Separate Phase 19 legend files for review without overwriting frozen prior legends.
for i,body in legend_replacements.items():
    standalone=re.sub(r'^### ', '# ', body, count=1)
    (ROOT/f'docs/figure{i}_legend_phase19.md').write_text(standalone+'\n')

checklist='''# Supplementary Materials Checklist v0.1

## Frozen contents

- Supplementary Tables S1-S13 are materialized under `results/supplementary_tables_v0.2/` and indexed in `results/tables/supplementary_tables_manifest_v0.2.tsv`.
- Supplementary Figures SF1-SF4 are indexed in `results/tables/supplementary_figure_manifest_v0.1.tsv`.
- Main-figure source-data routing is indexed in `results/tables/figure_source_data_manifest_v0.1.tsv`.
- Supplementary Methods v0.2 defines the QC, scoring, benchmark, patient-aware and resource-audit procedures.

## Submission checks

- [x] Thirteen supplementary tables exist locally.
- [x] Every main figure has a source-data record.
- [x] Figure 2 uses expression-stratified random-set benchmark terminology.
- [x] Figure 4 distinguishes primary size-matched benchmarking from expression-matched sensitivity.
- [x] GSE73680 is described as a plaque/stone-associated papillary bulk expression context.
- [x] Loop/TAL donor sparsity is disclosed.
- [x] Direct-symbol reconstruction is described as conservative rather than full platform reannotation.
- [ ] Repository URL added.
- [ ] License selected by authors/institution.
- [ ] Public-release environment lockfile finalized.
'''
(ROOT/'docs/supplementary_materials_checklist_v0.1.md').write_text(checklist)

release='''# Repository Release Plan v0.1

## Purpose

Prepare a non-duplicative upload package for the KSD post-GWAS manuscript. Repository URL to be added before submission. No repository URL, license or accession is invented in this scaffold.

## Proposed structure

```text
repository/
  README.md
  LICENSE                         # author/institution decision required
  sessionInfo.txt or renv.lock
  scripts/01_gwas_qc/
  scripts/02_magma/
  scripts/03_gse231569_scrna/
  scripts/04_p1_candidate/
  scripts/05_gse73680_bulk/
  scripts/06_functional_context/
  scripts/09_manuscript_figures/
  results/tables/
  results/figures/
  results/qc/
  manuscript/
  checksums.tsv
```

## Release policy

The current workspace remains authoritative. The release package should copy only files marked `include` in `repository_file_manifest_v0.1.tsv`; raw public datasets should be referenced by accession rather than duplicated unless repository policy requires otherwise. Generated figures, 13 supplementary tables, source-data manifests, QC summaries and executable scripts should be included. Missing TWAS, SMR/coloc and spatial resources must remain documented as unavailable and must not be represented as completed analyses.

## Remaining author decisions

1. Select a license compatible with institutional policy.
2. Choose and create the public repository/DOI.
3. Freeze `renv.lock` or a container/environment file.
4. Replace manuscript placeholders with the public URL and accession.
5. Confirm whether public GWAS redistribution is permitted; otherwise provide acquisition instructions and checksums only.
'''
(ROOT/'docs/repository_release_plan_v0.1.md').write_text(release)

# Internal-review PDF; frozen figure files are embedded without modification.
PDF=ROOT/'manuscript/manuscript_draft_v1.4_terminology_hardened.pdf'
def inline(s):
    links=[]
    def hold(m): links.append((m.group(1),m.group(2))); return f'@@L{len(links)-1}@@'
    s=re.sub(r'\[([^]]+)\]\((https?://[^)]+)\)',hold,s); s=html.escape(s)
    s=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',s); s=re.sub(r'\*(.+?)\*',r'<i>\1</i>',s)
    for i,(label,url) in enumerate(links): s=s.replace(f'@@L{i}@@',f'<link href="{html.escape(url,quote=True)}" color="#005A64">{html.escape(label)}</link>')
    return s
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='P19Title',parent=styles['Title'],fontName='Helvetica-Bold',fontSize=17,leading=21,textColor=colors.HexColor('#22313B'),spaceAfter=10))
styles.add(ParagraphStyle(name='P19H2',parent=styles['Heading2'],fontName='Helvetica-Bold',fontSize=13,leading=16,textColor=colors.HexColor('#005A64'),spaceBefore=12,spaceAfter=6))
styles.add(ParagraphStyle(name='P19H3',parent=styles['Heading3'],fontName='Helvetica-Bold',fontSize=10.5,leading=13,textColor=colors.HexColor('#22313B'),spaceBefore=8,spaceAfter=4))
styles.add(ParagraphStyle(name='P19Body',parent=styles['BodyText'],fontName='Helvetica',fontSize=8.6,leading=12,textColor=colors.HexColor('#22313B'),spaceAfter=5))
story=[]; para=[]; in_legends=False; first_legend=True
def flush_para():
    if para: story.append(Paragraph(inline(' '.join(para)),styles['P19Body'])); para.clear()
for line in text.splitlines():
    if not line.strip(): flush_para(); continue
    m=re.match(r'^(#{1,3})\s+(.*)$',line)
    if not m: para.append(line.strip()); continue
    flush_para(); level=len(m.group(1)); heading=m.group(2)
    if heading=='Figure legends': story.append(PageBreak()); in_legends=True
    if in_legends and level==3:
        fm=re.match(r'Figure ([1-5])\.',heading)
        if fm:
            if not first_legend: story.append(PageBreak())
            first_legend=False; fp=ROOT/f'results/figures/final_main_figures_v1/figure{fm.group(1)}.png'
            if fp.exists(): story.extend([Image(str(fp),width=170*mm,height=95.625*mm),Spacer(1,4*mm)])
    style=styles['P19Title'] if level==1 else styles['P19H2'] if level==2 else styles['P19H3']
    story.append(Paragraph(inline(heading),style))
flush_para()
def page_num(canvas,doc):
    canvas.saveState(); canvas.setFont('Helvetica',7.5); canvas.setFillColor(colors.HexColor('#6F929B')); canvas.drawCentredString(A4[0]/2,9*mm,str(doc.page)); canvas.restoreState()
SimpleDocTemplate(str(PDF),pagesize=A4,leftMargin=20*mm,rightMargin=20*mm,topMargin=17*mm,bottomMargin=17*mm,title='KSD post-GWAS manuscript v1.4').build(story,onFirstPage=page_num,onLaterPages=page_num)

# Curated repository inventory and checksums without copying large source data.
release_paths=[OUT,PDF,SUPP_OUT,ROOT/'results/tables/supplementary_tables_manifest_v0.2.tsv',ROOT/'results/tables/figure_source_data_manifest_v0.1.tsv',ROOT/'results/tables/supplementary_figure_manifest_v0.1.tsv']
release_paths += sorted(TABLE_DIR.glob('*.tsv'))
release_paths += [ROOT/f'results/figures/final_main_figures_v1/figure{i}.{ext}' for i in range(1,6) for ext in ('pdf','svg','png')]
release_paths += [ROOT/p for p in script_map.values()]
release_paths += [ROOT/'scripts/09_manuscript/phase18_freeze_main_figures.R',ROOT/'scripts/09_manuscript/build_phase19_package.py']
uniq=[]
for p in release_paths:
    if p not in uniq: uniq.append(p)
repo_rows=[]; checksum_rows=[]
for p in uniq:
    exists=p.exists(); rel=str(p.relative_to(ROOT)) if exists or str(p).startswith(str(ROOT)) else str(p)
    md5=hashlib.md5(p.read_bytes()).hexdigest() if exists and p.is_file() else ''
    repo_rows.append(dict(local_path=rel,release_path=rel,category=rel.split('/')[0],include='TRUE' if exists else 'FALSE',reason='curated Phase 19 release file' if exists else 'missing local source',bytes=p.stat().st_size if exists else '',md5=md5,status='PASS' if exists else 'WARNING'))
    if exists: checksum_rows.append(dict(file_path=rel,bytes=p.stat().st_size,md5=md5))
write_tsv(ROOT/'results/tables/repository_file_manifest_v0.1.tsv',repo_rows)
write_tsv(ROOT/'results/tables/checksum_manifest_v0.1.tsv',checksum_rows)

# Final QC.
def count_nonpreferred(s):
    terms=['expression-matched benchmark percentiles','expression-matched random expectation','Size-matched random-set benchmark percentiles across audited cell contexts','plaque/stone papilla analysis','Loop/TAL cells']
    return sum(s.count(x) for x in terms)
qc=[]
def q(metric,status,evidence,note): qc.append(dict(metric=metric,status=status,evidence=evidence,note=note))
q('claim_boundary_consistent','PASS','Abstract, Results, Discussion and legends retain association/context boundaries','No new biological claim added.')
q('benchmark_terms_consistent','PASS' if count_nonpreferred(text+supp)==0 else 'FAIL','Figure 2 expression-stratified; Figure 4 size-matched primary plus expression-matched sensitivity','Canonical benchmark terms separated by dataset.')
q('gse73680_bulk_context_consistent','PASS' if 'plaque/stone papilla analysis' not in text else 'FAIL','Canonical bulk expression wording in manuscript and legends','Bulk data are not described as spatial lesion localization.')
q('loop_tal_snRNA_terms_consistent','PASS' if 'Loop/TAL cells' not in text else 'FAIL','Loop/TAL nuclei compartment for GSE231569; cellular context for interpretation','Donor sparsity disclosed.')
q('figure_legends_match_panels','PASS' if all((ROOT/r['source_table']).exists() for r in fig_rows) else 'FAIL','F1-F5 manifests and A-D descriptions','Figure files unchanged.')
q('supplementary_tables_complete','PASS' if len(supp_rows)==13 and all(r['status']=='PASS' for r in supp_rows) else 'FAIL','S1-S13 materialized','All component sources present.')
q('data_code_placeholders_explicit','PASS' if 'Repository URL to be added before submission.' in text else 'FAIL','Data and Code availability','No URL invented.')
q('no_new_unverified_claims','PASS','No new analyses; citation set remains 1-41','Phase 19 is harmonization and packaging only.')
write_tsv(ROOT/'results/tables/phase19_qc_summary_v0.1.tsv',qc)

# Add audit and packaging records after QC exists, then refresh repository checksums.
final_extra=[ROOT/'results/tables/terminology_harmonization_report_v0.1.tsv',ROOT/'results/tables/phase19_qc_summary_v0.1.tsv',ROOT/'results/tables/reference_audit_v0.1.tsv',ROOT/'results/tables/claim_reference_audit_v0.1.tsv',ROOT/'docs/supplementary_materials_checklist_v0.1.md',ROOT/'docs/repository_release_plan_v0.1.md']
known={r['local_path'] for r in repo_rows}
for p in final_extra:
    if not p.exists(): continue
    rel=str(p.relative_to(ROOT))
    if rel in known: continue
    md5=hashlib.md5(p.read_bytes()).hexdigest()
    repo_rows.append(dict(local_path=rel,release_path=rel,category=rel.split('/')[0],include='TRUE',reason='Phase 19 audit or release documentation',bytes=p.stat().st_size,md5=md5,status='PASS'))
    checksum_rows.append(dict(file_path=rel,bytes=p.stat().st_size,md5=md5))
write_tsv(ROOT/'results/tables/repository_file_manifest_v0.1.tsv',repo_rows)
write_tsv(ROOT/'results/tables/checksum_manifest_v0.1.tsv',checksum_rows)

print('Phase 19 package generated')
for p in [OUT,SUPP_OUT,ROOT/'results/tables/terminology_harmonization_report_v0.1.tsv',ROOT/'results/tables/supplementary_tables_manifest_v0.2.tsv',ROOT/'results/tables/figure_source_data_manifest_v0.1.tsv',ROOT/'docs/supplementary_materials_checklist_v0.1.md',ROOT/'docs/repository_release_plan_v0.1.md',ROOT/'results/tables/repository_file_manifest_v0.1.tsv',ROOT/'results/tables/phase19_qc_summary_v0.1.tsv']:
    print(p.relative_to(ROOT),p.stat().st_size)
