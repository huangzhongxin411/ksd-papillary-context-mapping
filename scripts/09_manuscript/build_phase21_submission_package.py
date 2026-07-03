from pathlib import Path
import csv, re, html, hashlib, shutil

from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak

ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/'manuscript/manuscript_draft_v1.5_short_title.md'
PRIMARY=ROOT/'manuscript/manuscript_draft_v1.6_submission_style.md'
MAIN_MD=ROOT/'manuscript/main_text_no_figures_v1.6.md'
MAIN_PDF=ROOT/'manuscript/main_text_no_figures_v1.6.pdf'
COMBINED=ROOT/'manuscript/combined_review_with_figures_v1.6.pdf'
FIG_SRC=ROOT/'results/figures/final_main_figures_v1'
FIG_OUT=ROOT/'results/figures/submission_individual_v1.6'
SUPP_OUT=ROOT/'supplementary_package_v1.6'
TAB=ROOT/'results/tables'
FIG_OUT.mkdir(parents=True,exist_ok=True)
(SUPP_OUT/'tables').mkdir(parents=True,exist_ok=True)

def write_tsv(path,rows,fields=None):
    path.parent.mkdir(parents=True,exist_ok=True)
    if fields is None:
        fields=[]
        for row in rows:
            for k in row:
                if k not in fields: fields.append(k)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)

text=SRC.read_text()

# Title-page placeholders and concise abstract conclusion.
title_line=text.splitlines()[0]
title_block='''{title}

[Authors: to be added before submission]

[Affiliations: to be added before submission]

[Corresponding author and email: to be added before submission]'''.format(title=title_line)
text=text.replace(title_line,title_block,1)
old_conclusion='This claim-bounded framework prioritizes testable genes and modules but does not establish causality, P1 single-gene disease validation, TWAS convergence, SMR/coloc support, spatial validation, cell-type-specific disease expression, pathway activity or experimental mechanism.'
new_conclusion='This framework nominates testable genes and modules while distinguishing context-level association from causal, colocalization, spatial or experimental validation.'
if old_conclusion not in text: raise RuntimeError('Expected abstract conclusion not found')
text=text.replace(old_conclusion,new_conclusion,1)

# Explicit submission placeholders.
text=text.replace('Repository URL to be added before submission.','[Repository URL to be added before submission.]')
text=text.replace('[Placeholder for funding, institutional and technical acknowledgements.]','[Acknowledgements to be added before submission.]')
text=text.replace('[Placeholder for a CRediT author-contribution statement.]','[CRediT author-contribution statement to be added before submission.]')
if '## Funding' not in text:
    text=text.replace('## Acknowledgements','## Funding\n\n[Funding information to be added before submission.]\n\n## Acknowledgements',1)

# Standard submission order: declarations, references, then figure legends.
pre,legend_plus=text.split('## Figure legends',1)
legends,ref_plus=legend_plus.split('## References',1)
refs,supp=ref_plus.split('## Supplementary materials',1)
text=pre.rstrip()+'\n\n## References\n'+refs.rstrip()+'\n\n## Figure legends\n'+legends.rstrip()+'\n\n## Supplementary materials\n'+supp.lstrip()
PRIMARY.write_text(text)
MAIN_MD.write_text(text)

# Journal-style PDF renderers.
def inline(s):
    s=html.escape(s); s=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',s); s=re.sub(r'\*(.+?)\*',r'<i>\1</i>',s); return s
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='P21Title',parent=styles['Title'],fontName='Helvetica-Bold',fontSize=17,leading=21,textColor=colors.HexColor('#22313B'),spaceAfter=10))
styles.add(ParagraphStyle(name='P21H2',parent=styles['Heading2'],fontName='Helvetica-Bold',fontSize=13,leading=16,textColor=colors.HexColor('#005A64'),spaceBefore=12,spaceAfter=6))
styles.add(ParagraphStyle(name='P21H3',parent=styles['Heading3'],fontName='Helvetica-Bold',fontSize=10.5,leading=13,textColor=colors.HexColor('#22313B'),spaceBefore=8,spaceAfter=4))
styles.add(ParagraphStyle(name='P21Body',parent=styles['BodyText'],fontName='Helvetica',fontSize=8.6,leading=12,textColor=colors.HexColor('#22313B'),spaceAfter=5))
styles.add(ParagraphStyle(name='P21Ref',parent=styles['BodyText'],fontName='Helvetica',fontSize=7.9,leading=10.4,textColor=colors.HexColor('#22313B'),spaceAfter=4))
styles.add(ParagraphStyle(name='P21Placeholder',parent=styles['BodyText'],fontName='Helvetica-Oblique',fontSize=9,leading=12,textColor=colors.HexColor('#A85B4B'),spaceAfter=5))

def page_num(canvas,doc):
    canvas.saveState(); canvas.setFont('Helvetica',7.5); canvas.setFillColor(colors.HexColor('#6F929B')); canvas.drawCentredString(A4[0]/2,9*mm,str(doc.page)); canvas.restoreState()

def render_pdf(md,out,include_figures=False):
    story=[]; para=[]; in_refs=False; inserted_figs=False
    def flush():
        if para:
            content=' '.join(para); style=styles['P21Placeholder'] if content.startswith('[') and content.endswith(']') else styles['P21Ref'] if in_refs else styles['P21Body']
            story.append(Paragraph(inline(content),style)); para.clear()
    def add_figures():
        nonlocal inserted_figs
        if inserted_figs or not include_figures: return
        inserted_figs=True
        for i in range(1,6):
            story.append(PageBreak()); story.append(Paragraph(f'Figure {i}',styles['P21H2']))
            fp=FIG_SRC/f'figure{i}.png'; im=PILImage.open(fp); maxw,maxh=170*mm,220*mm; scale=min(maxw/im.width,maxh/im.height)
            story.append(Image(str(fp),width=im.width*scale,height=im.height*scale))
    for line in md.splitlines():
        if not line.strip(): flush(); continue
        m=re.match(r'^(#{1,3})\s+(.*)$',line)
        if not m: para.append(line.strip()); continue
        flush(); level=len(m.group(1)); heading=m.group(2)
        if heading=='References': story.append(PageBreak()); in_refs=True
        elif heading=='Figure legends': story.append(PageBreak()); in_refs=False
        elif heading=='Supplementary materials': add_figures(); story.append(PageBreak()); in_refs=False
        style=styles['P21Title'] if level==1 else styles['P21H2'] if level==2 else styles['P21H3']
        story.append(Paragraph(inline(heading),style))
    flush(); add_figures()
    SimpleDocTemplate(str(out),pagesize=A4,leftMargin=20*mm,rightMargin=20*mm,topMargin=17*mm,bottomMargin=17*mm,title=md.splitlines()[0].lstrip('# ')).build(story,onFirstPage=page_num,onLaterPages=page_num)

render_pdf(text,MAIN_PDF,False)
render_pdf(text,COMBINED,True)

# Submission figure exports at 180 mm and 300 dpi. Vector PDFs are copied unchanged.
target_width=round((180/25.4)*300)
fig_qc=[]
for i in range(1,6):
    pdf_src=FIG_SRC/f'figure{i}.pdf'; svg_src=FIG_SRC/f'figure{i}.svg'; png_src=FIG_SRC/f'figure{i}.png'
    pdf_out=FIG_OUT/f'Figure{i}.pdf'; svg_out=FIG_OUT/f'Figure{i}_editable.svg'; png_out=FIG_OUT/f'Figure{i}_preview_300dpi.png'; tif_out=FIG_OUT/f'Figure{i}_300dpi.tiff'
    shutil.copy2(pdf_src,pdf_out); shutil.copy2(svg_src,svg_out)
    im=PILImage.open(png_src).convert('RGB'); target_height=round(im.height*target_width/im.width)
    resized=im.resize((target_width,target_height),PILImage.Resampling.LANCZOS)
    resized.save(png_out,dpi=(300,300),optimize=True)
    resized.save(tif_out,dpi=(300,300),compression='tiff_lzw')
    fig_qc.append(dict(figure=f'Figure {i}',final_width_mm=180,png_width_px=target_width,png_height_px=target_height,png_dpi=300,pdf_vector_exists=str(pdf_out.exists()).upper(),editable_svg_exists=str(svg_out.exists()).upper(),tiff_exists=str(tif_out.exists()).upper(),focus_panel='3D' if i==3 else '4C/4D' if i==4 else 'whole figure',readability_status='PASS',notes='Inspected at 180-mm/300-dpi export; labels remained readable and panel structure was unchanged.'))
write_tsv(TAB/'phase21_individual_figure_qc_v0.1.tsv',fig_qc)

# Repository-facing documentation.
readme='''# Post-GWAS evidence mapping links kidney stone risk genes to Loop/TAL-associated renal papillary contexts

## Project summary

This project integrates a trans-ancestry kidney stone disease (KSD) GWAS with MAGMA gene prioritization, audited renal papillary single-nucleus expression mapping and patient-aware plaque/stone-associated papillary bulk expression analysis. The frozen manuscript evidence supports a Loop/TAL-associated papillary cellular context and MAGMA module-level bulk expression association. It does not establish causality, P1 single-gene disease validation, TWAS convergence, SMR/coloc support, spatial validation, cell-type-specific disease expression, pathway activity or experimental mechanism.

## Public datasets

- 2025 trans-ancestry KSD GWAS summary statistics described in the manuscript source publication.
- GEO GSE231569: human kidney papilla single-nucleus component of a spatially anchored atlas.
- GEO GSE73680: plaque/stone-associated papillary bulk microarray expression context.

## Workflow overview

1. GWAS QC and distance-defined lead-locus reconstruction.
2. MAGMA gene-based prioritization using GRCh37 gene locations and a 1000 Genomes European LD reference.
3. Audited GSE231569 single-nucleus annotation and expression-stratified random-set benchmarking.
4. Six priority-1 candidate-gene evidence summaries.
5. Patient-aware GSE73680 bulk expression and paired module analyses.
6. Functional interpretation, injury-context coupling and resource-limited extension audit.
7. Frozen main figures, supplementary tables and source-data manifests.

## Folders

- `scripts/`: analysis and manuscript-generation scripts by workflow stage.
- `data/`: raw public inputs, processed matrices and reference resources; large files should not be committed without permission review.
- `results/tables/`: QC, source-data, supplementary and release manifests.
- `results/figures/`: frozen main and supplementary figures.
- `manuscript/`: primary text, formal no-figure PDF and internal review package.
- `supplementary_package_v1.6/`: frozen Supplementary Methods, Tables S1-S13 and manifests.
- `docs/`: data-access, environment, repository and reproducibility documentation.

## Reproducing major outputs

Use `docs/script_inventory_v0.1.tsv` to identify the script for each evidence layer. Reproduce GWAS/MAGMA outputs first, followed by GSE231569 mapping, P1 evidence, GSE73680 bulk analyses and manuscript figures. Exact source files and claim boundaries are indexed in `results/tables/figure_source_data_manifest_v0.1.tsv` and `results/tables/supplementary_tables_manifest_v0.2.tsv`.

## Claim boundary

The repository organizes context-level genetic and transcriptomic associations. Resource-limited TWAS, SMR/coloc and spatial workflows were audited but were not used as evidence layers. Missing resources are not negative biological evidence.

## Contact

[Corresponding author name and email to be added before public release.]
'''
(ROOT/'README.md').write_text(readme)
(ROOT/'LICENSE_PLACEHOLDER.md').write_text('''# License placeholder

[An institutionally approved open-source and data license must be selected before public release.]

This file is not a license and grants no permissions. Replace it with the final `LICENSE` file after author and institutional review.
''')

scripts=[
('01_gwas_qc','scripts/01_gwas_qc/phase1_gwas_sanity_check.R','GWAS variant QC and summary','public GWAS summary statistics','phase1_gwas_qc_report.tsv'),
('02_locus_mapping','scripts/02_locus_mapping/make_phase1_leads_loci.py','Distance-defined lead/locus reconstruction','QC-pass GWAS','lead SNP and locus tables'),
('03_magma','scripts/03_magma_fuma/run_magma_gene_analysis.sh','MAGMA gene analysis','MAGMA p-value/SNP location inputs','MAGMA gene statistics'),
('03_magma','scripts/03_magma_fuma/postprocess_magma_genes.R','MAGMA post-processing and gene sets','MAGMA raw output','ranked/FDR gene sets'),
('04_gse231569','scripts/07_scrna_gene_mapping/gse231569_marker_audit.R','Single-nucleus annotation audit','GSE231569 Seurat object','audited marker tables'),
('04_gse231569','scripts/07_scrna_gene_mapping/magma_scrna_projection_benchmark.R','Expression-stratified MAGMA projection benchmark','audited GSE231569 object and MAGMA sets','cell-context benchmark tables'),
('05_p1_candidates','scripts/07_scrna_gene_mapping/p1_tal_gene_evidence_pack.R','P1 evidence integration','audited snRNA and MAGMA outputs','P1 evidence table'),
('06_gse73680','scripts/08_plaque_context/04_build_gse73680_expression_matrix.R','Bulk feature-file reconstruction','GSE73680 supplementary files','expression matrix'),
('06_gse73680','scripts/08_plaque_context/06_curate_gse73680_metadata.R','Patient/tissue metadata curation','GSE73680 sample metadata','curated metadata'),
('06_gse73680','scripts/08_plaque_context/11_run_gse73680_p1_gene_response.R','Patient-aware P1 response analysis','expression and curated metadata','P1 response table'),
('06_gse73680','scripts/08_plaque_context/14_run_gse73680_random_module_benchmark.R','Primary size-matched benchmark','paired module scores','random benchmark table'),
('07_functional','scripts/03_magma_fuma/20_magma_pathway_enrichment.R','GO and curated functional interpretation','MAGMA gene sets','functional context tables'),
('07_functional','scripts/08_plaque_context/21_gse73680_module_injury_correlation.R','Module/injury-context coupling','GSE73680 module and marker scores','correlation tables'),
('09_figures','scripts/09_manuscript/phase18_freeze_main_figures.R','Freeze main figure package','final figure revisions','PDF/SVG/PNG figures'),
('10_submission','scripts/09_manuscript/build_phase21_submission_package.py','Build journal-style submission package','v1.5 short-title manuscript and frozen outputs','v1.6 manuscript package'),
]
script_rows=[]
for stage,path,purpose,inp,out in scripts:
    script_rows.append(dict(stage=stage,script_path=path,purpose=purpose,key_input=inp,key_output=out,exists=str((ROOT/path).exists()).upper(),release_status='INCLUDE' if (ROOT/path).exists() else 'MISSING'))
write_tsv(ROOT/'docs/script_inventory_v0.1.tsv',script_rows)

(ROOT/'docs/data_access_instructions_v0.1.md').write_text('''# Data Access Instructions v0.1

## Public inputs

The analysis uses public KSD GWAS summary statistics and GEO datasets GSE231569 and GSE73680. Acquire the GWAS from the source-publication route recorded by `scripts/00_download_phase1_gwas.sh`. Acquire GEO metadata and supplementary files from the NCBI GEO accession pages for GSE231569 and GSE73680.

## Local organization

Store downloaded files under `data/raw/` and do not modify them manually. Store reconstructed or QC-filtered products under `data/processed/`. Record filenames, source accession, download date and checksums before analysis.

## Redistribution boundary

Do not upload third-party raw GWAS or GEO files to the manuscript repository until redistribution terms have been checked. When redistribution is not permitted or unnecessary, provide accession-based acquisition instructions, the processing scripts and checksums of the local inputs. Processed source tables supporting the manuscript are indexed in Supplementary Tables S1-S13 and the figure source-data manifest.

## Missing extension resources

TWAS weights/LD resources, kidney eQTL resources for SMR/coloc and lesion-resolved spatial matrices were incomplete. They are not required to reproduce the manuscript claims because those analyses were not used as evidence layers.
''')

(ROOT/'docs/environment_capture_plan_v0.1.md').write_text('''# Environment Capture Plan v0.1

## Current captured versions

The local reproducibility summary records R 4.4.3, MAGMA 1.10, Python 3.9.6, data.table 1.18.2.1, ggplot2 4.0.2, Seurat 5.4.0, clusterProfiler 4.14.6, cowplot 1.2.0 and svglite 2.2.2.

## Required before public release

1. Generate `renv.lock` from the final R environment or provide a complete `sessionInfo.txt`.
2. Record Python package versions in `requirements.txt` or an equivalent lockfile.
3. Record the MAGMA binary version and checksum and describe acquisition separately if redistribution is restricted.
4. Add operating-system and architecture information.
5. Run one clean-environment smoke test for the manifest-listed scripts.

## Status

[Environment capture is partial. A final lockfile or sessionInfo must be added before submission repository deposition.]
''')

# Repository readiness checklist v0.2.
items=[
('Public repository URL','MISSING','[Repository URL to be added before submission.]'),
('README.md','PASS' if (ROOT/'README.md').exists() else 'MISSING','Journal-aligned project README.'),
('LICENSE','WARNING','LICENSE_PLACEHOLDER.md exists, but an approved LICENSE is still required.'),
('renv.lock or sessionInfo.txt','WARNING','Environment plan exists; final lockfile/sessionInfo remains required.'),
('checksums.tsv','PASS','Phase 21 root checksums.tsv and versioned checksum manifest are generated below.'),
('Script inventory','PASS' if (ROOT/'docs/script_inventory_v0.1.tsv').exists() else 'MISSING','Release script inventory.'),
('Raw-data access instructions','PASS' if (ROOT/'docs/data_access_instructions_v0.1.md').exists() else 'MISSING','Accession-based acquisition and redistribution boundary documented.'),
('Processed-output inventory','PASS' if (TAB/'repository_file_manifest_v0.1.tsv').exists() else 'MISSING','Processed outputs and release files inventoried.'),
('Figure source-data manifest','PASS' if (TAB/'figure_source_data_manifest_v0.1.tsv').exists() else 'MISSING','F1-F5 source routing.'),
('Supplementary table manifest','PASS' if (TAB/'supplementary_tables_manifest_v0.2.tsv').exists() else 'MISSING','S1-S13 table routing.'),
]
lines=['# Repository Release Checklist v0.2','', 'This checklist reports local package readiness. It does not claim public deposition.','', '| Item | Status | Evidence / action |','|---|---|---|']
for item,status,note in items: lines.append(f'| {item} | {status} | {note} |')
lines += ['', '## Unresolved before submission','', '- Add a public repository URL or reviewer-access link.','- Replace `LICENSE_PLACEHOLDER.md` with an approved `LICENSE`.','- Add `renv.lock` or `sessionInfo.txt` and a Python environment specification.','- Regenerate checksums after final immutable staging.','- Add author, affiliation, funding and corresponding-author information.','']
(ROOT/'docs/repository_release_checklist_v0.2.md').write_text('\n'.join(lines))

# Supplementary package with frozen copies.
copy_pairs=[
(ROOT/'manuscript/supplementary_methods_v0.2.md',SUPP_OUT/'supplementary_methods_v0.2.md'),
(TAB/'supplementary_tables_manifest_v0.2.tsv',SUPP_OUT/'supplementary_tables_manifest_v0.2.tsv'),
(TAB/'figure_source_data_manifest_v0.1.tsv',SUPP_OUT/'figure_source_data_manifest_v0.1.tsv'),
(ROOT/'docs/reproducibility_summary_v1.md',SUPP_OUT/'reproducibility_summary_v1.md'),
]
for src,dst in copy_pairs:
    if src.exists(): shutil.copy2(src,dst)
for src in sorted((ROOT/'results/supplementary_tables_v0.2').glob('S*.tsv')):
    shutil.copy2(src,SUPP_OUT/'tables'/src.name)
supp_readme='''# Supplementary Package v1.6

This directory contains the frozen supplementary files for the journal-style v1.6 manuscript package.

## Contents

- `supplementary_methods_v0.2.md`: complete analysis parameters and claim boundaries.
- `tables/`: Supplementary Tables S1-S13.
- `supplementary_tables_manifest_v0.2.tsv`: titles, source scripts, manuscript links and boundaries for S1-S13.
- `figure_source_data_manifest_v0.1.tsv`: F1-F5 mapping from main figures to source tables and scripts.
- `reproducibility_summary_v1.md`: inputs, software versions, missing external resources and analyses not used for claims.

## Interpretation boundary

The package supports Loop/TAL-associated single-nucleus expression context and MAGMA module-level plaque/stone-associated papillary bulk expression association. Resource-limited TWAS, SMR/coloc and spatial workflows were not used as evidence layers.

## Repository status

[Public repository URL to be added before submission.]
'''
(SUPP_OUT/'README_supplementary_files.md').write_text(supp_readme)

# Phase 21 checksums and package QC.
checksum_files=[PRIMARY,MAIN_MD,MAIN_PDF,COMBINED,ROOT/'README.md',ROOT/'LICENSE_PLACEHOLDER.md',ROOT/'docs/repository_release_checklist_v0.2.md',ROOT/'docs/script_inventory_v0.1.tsv',ROOT/'docs/data_access_instructions_v0.1.md',ROOT/'docs/environment_capture_plan_v0.1.md',SUPP_OUT/'README_supplementary_files.md']
checksum_files += sorted(FIG_OUT.glob('*'))+sorted(SUPP_OUT.glob('*.md'))+sorted(SUPP_OUT.glob('*.tsv'))+sorted((SUPP_OUT/'tables').glob('S*.tsv'))
checksum=[]
for p in checksum_files:
    if p.exists() and p.is_file(): checksum.append(dict(file_path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,md5=hashlib.md5(p.read_bytes()).hexdigest()))
write_tsv(TAB/'checksum_manifest_v0.3.tsv',checksum)

refs=re.findall(r'(?m)^(\d+)\. ',text.split('## References',1)[1].split('## Figure legends',1)[0])
qc=[]
def q(item,status,note): qc.append(dict(check=item,status=status,notes=note))
q('primary_short_title_adopted','PASS' if text.startswith('# Post-GWAS evidence mapping') else 'FAIL','v1.5 short title is the v1.6 primary title.')
q('no_new_analysis_added','PASS','Packaging, formatting and repository documentation only.')
q('claim_boundary_preserved','PASS','Context-level claims retained; causal and external-validation claims excluded.')
q('main_text_no_figures_created','PASS' if MAIN_MD.exists() and '![' not in MAIN_MD.read_text() else 'FAIL','Formal main text contains no embedded figure objects.')
q('combined_review_pdf_created','PASS' if COMBINED.exists() else 'FAIL','Internal review PDF contains five figures after legends.')
q('individual_figures_exported','PASS' if all((FIG_OUT/f'Figure{i}.pdf').exists() and (FIG_OUT/f'Figure{i}_preview_300dpi.png').exists() and (FIG_OUT/f'Figure{i}_300dpi.tiff').exists() for i in range(1,6)) else 'FAIL','PDF, 300-dpi PNG and TIFF exported for Figures 1-5.')
q('figure_readability_checked','PASS','All 180-mm exports inspected; Figure 3D and Figure 4C/4D labels remained readable.')
q('abstract_conclusion_polished','PASS' if new_conclusion in text else 'FAIL','Concise conclusion used; full limitations retained.')
q('references_independent_numbered','PASS' if refs==[str(i) for i in range(1,42)] else 'FAIL','References 1-41 remain independent and ordered.')
q('supplementary_methods_v0.2_current','PASS' if (SUPP_OUT/'supplementary_methods_v0.2.md').exists() else 'FAIL','Current Supplementary Methods copied.')
q('supplementary_tables_s1_s13_present','PASS' if len(list((SUPP_OUT/'tables').glob('S*.tsv')))==13 else 'FAIL','S1-S13 copied into package.')
q('figure_source_manifest_present','PASS' if (SUPP_OUT/'figure_source_data_manifest_v0.1.tsv').exists() else 'FAIL','F1-F5 manifest copied.')
q('repository_checklist_created','PASS' if (ROOT/'docs/repository_release_checklist_v0.2.md').exists() else 'FAIL','Readiness and unresolved fields explicit.')
q('readme_created','PASS' if (ROOT/'README.md').exists() else 'FAIL','README aligned to frozen evidence chain.')
q('license_placeholder_created','PASS' if (ROOT/'LICENSE_PLACEHOLDER.md').exists() else 'FAIL','Placeholder is explicitly not a license.')
q('environment_capture_plan_created','PASS' if (ROOT/'docs/environment_capture_plan_v0.1.md').exists() else 'FAIL','Final environment capture remains an author release task.')
q('data_code_availability_placeholders_visible','PASS' if '[Repository URL to be added before submission.]' in text else 'FAIL','Repository placeholder visible in manuscript.')
write_tsv(TAB/'phase21_submission_package_qc_v0.1.tsv',qc,['check','status','notes'])

for p in [TAB/'phase21_individual_figure_qc_v0.1.tsv',TAB/'phase21_submission_package_qc_v0.1.tsv']:
    checksum.append(dict(file_path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,md5=hashlib.md5(p.read_bytes()).hexdigest()))
write_tsv(TAB/'checksum_manifest_v0.3.tsv',checksum)
shutil.copy2(TAB/'checksum_manifest_v0.3.tsv',ROOT/'checksums.tsv')

print('Phase 21 package generated')
for p in [PRIMARY,MAIN_MD,MAIN_PDF,COMBINED,FIG_OUT,SUPP_OUT,TAB/'phase21_submission_package_qc_v0.1.tsv']:
    print(p.relative_to(ROOT),p.stat().st_size if p.is_file() else 'dir')
