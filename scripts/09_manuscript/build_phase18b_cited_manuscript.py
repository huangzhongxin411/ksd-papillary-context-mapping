from pathlib import Path
import csv, re, xml.etree.ElementTree as ET
import html

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "manuscript/manuscript_draft_v1.0_pre.md"
XML = ROOT / "manuscript/references/pubmed_verified_phase18b.xml"
OUT = ROOT / "manuscript/manuscript_draft_v1.2_phase18b_cited.md"
REF_TSV = ROOT / "manuscript/references/reference_verification_phase18b.tsv"
MAP_TSV = ROOT / "manuscript/references/claim_reference_map_phase18b.tsv"
SUPP = ROOT / "docs/supplementary_materials_scaffold_phase18b.md"
RIS = ROOT / "manuscript/references/references_phase18b.ris"
BROWSER = ROOT / "manuscript/references/citation_review_phase18b.html"

order = [
"27188687","16443041","34907378","40216741","37980427","38896256","19561606","30523390",
"12618515","30066583","27297950","27731368","33514941","27441596","37468493","37468583",
"31249312","28628195","24948743","22891322","31205055","26107257","28274923","25885710",
"29632380","29785013","25605792","34062119","34557778","26854917","27019110","24830394",
"11752295","32533118","20962745","31729369"]

root = ET.parse(XML).getroot()
records = {}
for a in root.findall(".//PubmedArticle"):
    pmid = a.findtext(".//MedlineCitation/PMID")
    art = a.find(".//MedlineCitation/Article")
    journal = art.findtext("./Journal/ISOAbbreviation") or art.findtext("./Journal/Title") or ""
    title = "".join(art.find("./ArticleTitle").itertext()) if art.find("./ArticleTitle") is not None else ""
    year = art.findtext("./Journal/JournalIssue/PubDate/Year") or art.findtext("./ArticleDate/Year") or ""
    volume = art.findtext("./Journal/JournalIssue/Volume") or ""
    issue = art.findtext("./Journal/JournalIssue/Issue") or ""
    pages = art.findtext("./Pagination/MedlinePgn") or ""
    authors=[]
    for au in art.findall("./AuthorList/Author"):
        collective=au.findtext("CollectiveName")
        if collective: authors.append(collective); continue
        last=au.findtext("LastName") or ""; ini=au.findtext("Initials") or ""
        if last: authors.append(f"{last} {ini}".strip())
    doi=""
    for eid in art.findall("./ELocationID"):
        if eid.attrib.get("EIdType")=="doi": doi=eid.text or ""
    if not doi:
        for aid in a.findall(".//PubmedData/ArticleIdList/ArticleId"):
            if aid.attrib.get("IdType")=="doi": doi=aid.text or ""
    records[pmid]=dict(pmid=pmid,title=title,journal=journal,year=year,volume=volume,issue=issue,pages=pages,doi=doi,authors=authors)

missing=[x for x in order if x not in records]
if missing: raise RuntimeError(f"Missing verified PubMed records: {missing}")

def ref_line(n, r):
    aa=r["authors"]
    auth=", ".join(aa[:6]) + (", et al." if len(aa)>6 else ".")
    vi=r["volume"] + (f"({r['issue']})" if r["issue"] else "")
    tail=f" {r['year']};{vi}" if vi else f" {r['year']}"
    if r["pages"]: tail += f":{r['pages']}"
    tail += "."
    doi=f" https://doi.org/{r['doi']}" if r["doi"] else ""
    return f"{n}. {auth} {r['title']} {r['journal']}.{tail}{doi} [PubMed](https://pubmed.ncbi.nlm.nih.gov/{r['pmid']}/)"

refs="\n".join(ref_line(i,records[p]) for i,p in enumerate(order,1))

title = "Post-GWAS mapping links kidney stone risk genes to Loop/TAL-associated renal papillary cellular and disease-context programs"
abstract = """## Abstract

**Background:** Kidney stone disease (KSD) has a substantial polygenic component, but translating genome-wide association signals into renal papillary cellular and disease contexts remains challenging.

**Methods:** We integrated public trans-ancestry KSD GWAS summary statistics, MAGMA gene-based prioritization, audited GSE231569 renal papillary single-nucleus annotations, six-gene P1 candidate evidence scoring and patient-aware GSE73680 plaque/stone papilla expression analysis.

**Results:** MAGMA prioritized 94 Bonferroni-significant genes across 57 independent loci. Prioritized gene sets projected preferentially to Loop/TAL cells in the audited GSE231569 atlas, with benchmark percentiles of 0.998-1.000 for the top-ranked modules and preserved support in locus-balanced and leave-one-locus-out analyses. Six P1 genes formed a TAL, transport, calcium/ion-handling and broader epithelial role spectrum rather than a uniform TAL marker or disease-response panel. In GSE73680, comprising 55 samples and 26 paired patients, MAGMA modules, but not individual P1 genes, showed patient-aware plaque/stone papilla disease-context shifts. MAGMA top-50 scores also coupled with injury/remodeling programs in paired-delta (Spearman rho = 0.82, FDR = 4.5 x 10^-6) and patient/group residual analyses (rho = 0.79, FDR = 7.2 x 10^-12).

**Conclusion:** MAGMA-prioritized KSD genes support a Loop/TAL-associated papillary cellular context and module-level plaque/stone disease-context association. This claim-bounded framework nominates genes and modules for future testing but does not establish causal mediation, P1 single-gene disease validation, TWAS convergence, SMR/coloc support, spatial validation or external molecular validation.

**Keywords:** kidney stone disease; MAGMA; renal papilla; Loop of Henle; thick ascending limb; single-nucleus RNA sequencing"""

intro = """## Introduction

Kidney stone disease is common, recurrent and biologically heterogeneous, reflecting interactions among urinary chemistry, renal epithelial transport, papillary microenvironments and inherited susceptibility.[1-3] Genome-wide association studies have expanded the number of KSD susceptibility loci across ancestries, including loci near genes involved in calcium handling and tubular transport.[4-7] However, most association signals do not directly identify the effector gene, relevant renal cell type or disease-stage tissue context, leaving a substantial post-GWAS interpretation gap.[3,8]

The renal papilla is a plausible tissue context for stone initiation and growth. Randall plaque develops within papillary interstitium and can originate near basement membranes of thin limbs, while plaque-associated tissue shows inflammatory, oxidative-stress and injury-related molecular features.[9-14] Recent single-nucleus and spatially anchored transcriptomic studies have further resolved epithelial, stromal, vascular and immune niches in human kidney and papilla.[15-17] Within this anatomy, the Loop of Henle and thick ascending limb (Loop/TAL) coordinate salt transport, urinary concentration and paracellular calcium and magnesium handling through uromodulin-, claudin- and calcium-sensing pathways.[18-23] These functions make Loop/TAL-associated programs biologically plausible contexts for KSD genetic risk, but biological plausibility alone does not establish causal cell-type mediation.

Post-GWAS gene-based and cell-type enrichment approaches can connect polygenic association signals to tissues and cellular expression programs.[24-26] Single-cell and single-nucleus resources are particularly useful for localization, provided that annotation quality, gene-set background and donor structure are considered. Such analyses remain contextual: they do not substitute for transcriptome-wide association testing, summary-data Mendelian randomization, genetic colocalization or spatial validation.[30-32]

We therefore asked three questions. First, do MAGMA-prioritized KSD genes converge on an audited renal papillary cellular context? Second, do six P1 candidates form a coherent and interpretable TAL, transport and calcium/ion-handling spectrum? Third, are the prioritized modules reflected in an independent plaque/stone papilla disease-context dataset? To address these questions, we integrated the 2025 trans-ancestry KSD GWAS,[4] GSE231569 single-nucleus data,[15] and GSE73680 plaque-associated microarray data.[11] The study explicitly separates supported cellular and disease-context association from causality, TWAS convergence, SMR/coloc support, spatial validation and P1 single-gene disease validation."""

methods = """## Methods

### Study design and public resources

The analysis integrated the 2025 trans-ancestry KSD GWAS,[4] the GSE231569 renal papillary single-nucleus dataset and its associated spatially anchored atlas,[15] and GSE73680 Randall plaque expression profiles.[11,33] The main workflow comprised GWAS quality control, MAGMA prioritization, audited cell-context projection, P1 candidate interpretation, patient-aware plaque/stone disease-context analysis, functional enrichment and injury/remodeling coupling. Detailed file-level reconstruction audits and resource manifests are assigned to Supplementary Methods and Tables.

### GWAS quality control and MAGMA prioritization

Summary statistics were checked for variant identifiers, coordinates, alleles, effect estimates and valid P values before lead-locus construction. Genome-wide significant variants were distance-pruned into independent lead signals and loci. MAGMA v1.10 was used for GRCh37/hg19-compatible gene-based analysis,[24] and fixed top-50, top-100, FDR-significant and suggestive gene sets were carried forward. Detailed row counts, duplicate handling, lead-SNP parameters and QC diagnostics are reported in Supplementary Tables 1-3.

### GSE231569 annotation audit and gene-set projection

GSE231569 nuclei were analyzed using audited broad cell-type labels informed by marker expression and renal epithelial identity.[15-17] Loop of Henle and thick ascending limb labels were harmonized as Loop/TAL. MAGMA gene sets were scored across detected genes, summarized by audited cell type and compared with size-matched random sets. Locus-balanced, conservative-annotation and leave-one-locus-out analyses assessed robustness. Cell-level UMAP module scores were generated for visualization, whereas donor-cell-type summaries and random benchmarks provided the inferential layer. The workflow followed established principles for single-cell integration and cell-type enrichment without treating localization as causal mediation.[25,26,28]

### P1 candidate evidence spectrum

UMOD, CASR, CLDN14, CLDN10, HIBADH and PKD2 were evaluated for TAL expression rank, mean expression, detection frequency, donor-level detection, TAL specificity, TAL-program correlation and biological role. Roles were interpreted in the context of TAL transport, claudin biology, mineral handling and broader epithelial expression.[18-23,35,36]

### GSE73680 patient-aware disease-context analysis

GSE73680 Agilent feature files were reconstructed and linked to patient and tissue metadata.[11] After prespecified inclusion filtering, 55 samples from 29 patients, including 26 paired patients, entered analysis. Repeated sampling was handled with limma duplicateCorrelation,[27] with paired patient-level deltas as sensitivity analyses. P1 single-gene responses were evaluated separately from predefined MAGMA and contextual modules. Size-matched random sets, leave-one-gene-out analysis, MAGMA-without-P1 sensitivity and paired direction consistency assessed robustness. File parsing, feature mapping and reconstruction counts are reported in Supplementary Tables 7-9.

### Functional enrichment and injury/remodeling coupling

GO Biological Process enrichment used clusterProfiler with the MAGMA-tested gene universe as background.[29] Redundant terms were reduced for display, and curated nephron marker sets were assessed by hypergeometric testing. GSE73680 injury/remodeling, inflammation/immune, fibrosis/ECM and epithelial-injury programs were summarized as mean z-scored modules. Coupling with MAGMA modules was evaluated by sample-level, paired-delta and patient/group residual Spearman correlations, with paired and residual estimates prioritized because they better address disease-group composition and repeated samples.

### Resource-limited extensions and statistical boundary

TWAS, SMR/coloc and spatial pipelines were scaffolded and resource-audited but were not used as evidence layers because required external weights, eQTL/LD resources or lesion-resolved matrices remained incomplete.[30-32] Multiple testing was controlled within each analysis family using false-discovery-rate procedures. Software versions, full parameterization and resource status are provided in the reproducibility summary and Supplementary Table 12."""

src=SRC.read_text()
results=src.split("## Results",1)[1].split("## Discussion",1)[0].strip()
results = results.replace("The primary public KSD GWAS reconstruction retained 57 independent loci and was carried forward into MAGMA gene-level prioritization.", "MAGMA prioritization defined a reproducible KSD genetic module layer from 57 independent loci.", 1)
results = results.replace("We next projected MAGMA-prioritized KSD gene sets onto the audited GSE231569 single-nucleus atlas (Fig. 2a-d).", "MAGMA-prioritized KSD genes showed a consistent Loop/TAL-associated expression-context pattern in the audited GSE231569 renal papillary atlas (Fig. 2A-D).", 1)
results = results.replace("The six P1 candidate genes did not behave as a uniform TAL-specific marker set.", "The six P1 candidates formed a heterogeneous TAL-to-broader-epithelial role spectrum rather than a uniform TAL marker set.", 1)
results = results.replace("We next examined whether the MAGMA- and single-nucleus-prioritized signals were reflected in an independent papillary plaque/stone papilla disease-context dataset.", "GSE73680 supported the prioritized signal at the MAGMA module level rather than as a uniform P1 single-gene response.", 1)
results = results.replace("Cell-level GSE231569 score projection supported the Fig. 2 localization pattern on the audited UMAP, with MAGMA top-ranked module scores visibly enriched in the Loop/TAL-associated region.", "Functional enrichment and injury coupling linked the prioritized modules to renal tubular and papillary remodeling contexts without establishing pathway activity or mechanism.", 1)

discussion = """## Discussion

### A claim-bounded post-GWAS framework

This study adds an auditable framework for connecting KSD genetic prioritization to renal papillary cellular and disease contexts. The contribution is not the declaration of a single causal gene or cell type. Instead, MAGMA-defined genetic modules were localized to an audited Loop/TAL-associated single-nucleus context, interpreted through a heterogeneous P1 gene spectrum and evaluated in an independent plaque/stone papilla expression dataset. This layered design addresses a recurring post-GWAS problem: gene-based association and cell-type enrichment can organize evidence, but they remain distinct from causal mediation, TWAS, colocalization and experimental validation.[24-26,30-32]

### Why Loop/TAL is biologically relevant

The Loop/TAL localization is compatible with established renal physiology. TAL epithelia coordinate urinary concentration and transepithelial salt transport, while claudin-dependent paracellular pathways influence calcium and magnesium handling.[18-20,23] The P1 spectrum also aligns with this context: UMOD represents a TAL-associated protein with links to mineral metabolism,[21] CLDN10 and CLDN14 connect epithelial junction biology to ion handling and nephrocalcinosis or stone susceptibility,[7,19,20] and CASR links extracellular calcium sensing to stone-related phenotypes.[22,35,36] HIBADH and PKD2 were retained as broader supporting epithelial-context genes rather than forced into a uniform TAL marker class. Thus, the data support a biologically coherent cellular context while preserving heterogeneity among candidate genes.

### Why module-level disease-context support matters

The absence of FDR-supported P1 single-gene responses in GSE73680 does not conflict with the MAGMA module-level association. KSD susceptibility is polygenic, and individual risk genes need not show uniform differential expression in bulk lesion-associated tissue.[3,4,34] Module aggregation can capture coordinated but modest expression shifts distributed across prioritized genes, whereas single-gene testing is more vulnerable to limited effect size, tissue heterogeneity and multiple-testing burden. The patient-aware paired shifts, random-set benchmarks and MAGMA-without-P1 sensitivity therefore support disease-context association at the module level, not validation of any single P1 gene.

### Papillary injury context and future validation

Randall plaque and stone-forming papilla show inflammatory, immune and remodeling features,[9-14] and the GSE231569 atlas independently described immune injury and matrix remodeling in stone disease.[15] The observed coupling between MAGMA modules and injury/remodeling scores is consistent with this tissue context, but correlation does not establish that genetic risk causes injury or that injury mediates genetic effects. Future work should prioritize papilla-specific eQTL resources, genetic colocalization, lesion-resolved spatial transcriptomics and perturbation of prioritized modules in TAL-relevant epithelial models. These experiments would test whether the mapped cellular context corresponds to causal regulatory mechanisms or spatially localized disease programs."""

limitations = """## Limitations

This study is computational and hypothesis-generating. GSE231569 supports single-nucleus expression-context mapping but not spatial validation or causal cell-type mediation. GSE73680 is a bulk microarray dataset and therefore cannot establish cell-type-specific disease responses. P1 single-gene differential expression was not FDR-supported, and the PKD2 signal remains nominal and exploratory. Functional enrichment and risk-injury coupling are gene-set or module-level interpretations that may depend on gene-set size, background definition, expression detectability and disease-group composition. TWAS, SMR/coloc and spatial validation were not completed because the required external resources were unavailable in the current analysis environment."""

figure_legends = """## Figure legends

### Figure 1. Post-GWAS evidence map for kidney stone papillary cellular context

The framework links four claim-bounded evidence layers: KSD GWAS/MAGMA prioritization, GSE231569 single-nucleus localization, P1 candidate interpretation and GSE73680 plaque/stone papilla analysis. The Manhattan inset summarizes the cleaned KSD GWAS, and the remaining cards summarize Loop/TAL localization, the six-gene role spectrum and paired module-level disease-context shifts. The papillary inset provides anatomical orientation only. Arrows indicate evidence flow rather than causal progression. Supported inferences are separated from unresolved causality, TWAS convergence, SMR/coloc support, spatial validation and P1 single-gene disease validation.

### Figure 2. MAGMA-prioritized KSD genes localize to a Loop/TAL-associated single-nucleus context

(A) KSD GWAS Manhattan summary with representative prioritized loci. (B) MAGMA top-50 module projection on the audited GSE231569 UMAP. (C) Size-matched random-set benchmark percentiles across audited cell contexts; the dashed line marks the 95th percentile. (D) Leading Loop/TAL module contributors. Cell-level maps provide visualization, whereas donor-cell-type summaries and benchmark analyses provide the inferential layer. The figure supports a Loop/TAL-associated expression context but not causal cell-type mediation, TWAS convergence, SMR/coloc support or spatial validation.

### Figure 3. Gene-centric evidence spectrum of P1 candidates

(A) Six P1 genes arranged from TAL identity to broader epithelial context. (B) Expression and detection across six audited GSE231569 contexts. (C) TAL specificity ratios and donor support; the dashed line marks two-fold TAL preference. (D) Evidence matrix summarizing MAGMA prioritization, TAL rank, donor support, specificity, GSE73680 bulk-response status and interpretation role. PKD2 was nominal only, and no P1 gene showed a uniform FDR-supported paired bulk response. The panel is interpretive and does not establish P1 single-gene disease validation.

### Figure 4. GSE73680 supports MAGMA module-level plaque/stone disease-context association

(A) Patient-aware GSE73680 design. (B) Paired deltas for six P1 genes; no gene reached FDR q < 0.05. (C) Paired module-score shifts for MAGMA gene sets and P1 core. (D) Size-matched random-set benchmarks and robustness summaries. MAGMA modules were FDR-supported, whereas P1 core was not. These analyses support module-level disease-context association rather than P1 single-gene validation, causality or cell-type-specific disease expression.

### Figure 5. Functional context, injury coupling and extension audit

(A) Evidence network linking prioritized modules to TAL transport, calcium/ion handling, epithelial junction and injury/remodeling themes. Link width is relative within each evidence layer and is not a causal effect scale. (B) Ranked redundancy-reduced GO terms. (C) Paired-delta correlations between risk modules and injury contexts. (D) Robustness summaries and resource-audit cards. Functional and plaque/stone contexts support module-level interpretation but do not establish pathway activity, causal mediation, TWAS convergence, SMR/coloc support, spatial validation or external molecular validation."""

tail = """## Data availability

This study used public GWAS summary statistics and GEO datasets GSE231569 and GSE73680.[4,11,15,33] Processed tables, QC summaries and figure-source outputs are organized with the analysis code and will require a persistent public repository and accession before submission.

## Code availability

Analysis scripts are organized by workflow stage and can be released with processed outputs at submission. A permanent repository URL and software environment lockfile remain to be added.

## Acknowledgements

[Placeholder for funding, institutional and technical acknowledgements.]

## Author contributions

[Placeholder for a CRediT author-contribution statement.]

## Competing interests

The authors declare no competing interests. [Author confirmation required.]

""" + figure_legends + """

## References

""" + refs + "\n\n## Supplementary materials\n\nThe supplementary package is scaffolded as Supplementary Tables 1-13 and associated methods notes; see `docs/supplementary_materials_scaffold_phase18b.md`."

final = f"# {title}\n\n{abstract}\n\n{intro}\n\n{methods}\n\n## Results\n\n{results}\n\n{discussion}\n\n{limitations}\n\n{tail}\n"
OUT.write_text(final)

with REF_TSV.open("w",newline="") as f:
    fields=["reference_number","pmid","doi","title","journal","year","verification_source","status"]
    w=csv.DictWriter(f,fields,delimiter="\t"); w.writeheader()
    for n,p in enumerate(order,1):
        r=records[p]; w.writerow(dict(reference_number=n,pmid=p,doi=r["doi"],title=r["title"],journal=r["journal"],year=r["year"],verification_source="NCBI PubMed XML",status="verified"))

ris=[]
for p in order:
    r=records[p]; ris += ["TY  - JOUR"]
    ris += [f"AU  - {a}" for a in r["authors"]]
    ris += [f"TI  - {r['title']}", f"JO  - {r['journal']}", f"PY  - {r['year']}"]
    if r["volume"]: ris.append(f"VL  - {r['volume']}")
    if r["issue"]: ris.append(f"IS  - {r['issue']}")
    if r["pages"]: ris.append(f"SP  - {r['pages']}")
    if r["doi"]: ris.append(f"DO  - {r['doi']}")
    ris += [f"UR  - https://pubmed.ncbi.nlm.nih.gov/{p}/", f"AN  - PMID:{p}", "ER  - ", ""]
RIS.write_text("\n".join(ris))

maps=[
("C01","KSD epidemiology, recurrence and inherited susceptibility","1-8","background/strong"),
("C02","Randall plaque and renal papillary disease niche","9-14","strong/background"),
("C03","Human kidney and papilla single-nucleus context","15-17","strong"),
("C04","Loop/TAL transport and calcium/ion handling","18-23,35-36","strong/background"),
("C05","MAGMA and post-GWAS cell-type enrichment methodology","24-26","methods"),
("C06","limma, Seurat and clusterProfiler implementation","27-29","methods"),
("C07","TWAS, SMR and colocalization are distinct evidence layers","30-32","methods/boundary"),
("C08","GEO resource provenance","33","database"),
("C09","Polygenic/module-level interpretation of KSD risk","3-6,34","background")]
with MAP_TSV.open("w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["claim_id","claim","reference_numbers","support_level"]); w.writerows(maps)

rows=[]
for claim_id,claim,refnums,support in maps:
    rows.append(f"<tr><td>{claim_id}</td><td>{html.escape(claim)}</td><td>{refnums}</td><td>{support}</td></tr>")
refrows=[]
for n,p in enumerate(order,1):
    r=records[p]; doi=(f"<a href='https://doi.org/{html.escape(r['doi'])}'>{html.escape(r['doi'])}</a>" if r['doi'] else "No DOI in PubMed")
    refrows.append(f"<tr><td>{n}</td><td>{html.escape(r['title'])}</td><td>{html.escape(r['journal'])} ({r['year']})</td><td>{doi}</td><td><a href='https://pubmed.ncbi.nlm.nih.gov/{p}/'>{p}</a></td></tr>")
BROWSER.write_text("""<!doctype html><meta charset='utf-8'><title>Phase 18B citation review</title><style>body{font-family:Arial,sans-serif;color:#22313B;max-width:1300px;margin:32px auto;line-height:1.45}h1,h2{color:#005A64}table{border-collapse:collapse;width:100%;margin-bottom:30px}th,td{border:1px solid #D3DADC;padding:8px;vertical-align:top}th{background:#EAF1F2;text-align:left}a{color:#005A64}</style><h1>Phase 18B citation review</h1><p>All records were verified from NCBI PubMed XML. Support levels describe how each group is used in the manuscript.</p><h2>Claim-reference map</h2><table><tr><th>ID</th><th>Claim</th><th>References</th><th>Support</th></tr>"""+"".join(rows)+"</table><h2>Verified references</h2><table><tr><th>No.</th><th>Title</th><th>Journal/year</th><th>DOI</th><th>PubMed</th></tr>"+"".join(refrows)+"</table>")

SUPP.write_text("""# Phase 18B Supplementary Materials Scaffold

1. Supplementary Table 1 - GWAS QC and lead-locus summary
2. Supplementary Table 2 - MAGMA gene-level results
3. Supplementary Table 3 - MAGMA gene-set definitions
4. Supplementary Table 4 - GSE231569 cell-type annotation audit
5. Supplementary Table 5 - Loop/TAL benchmark and robustness results
6. Supplementary Table 6 - P1 candidate evidence matrix
7. Supplementary Table 7 - GSE73680 reconstruction and metadata audit
8. Supplementary Table 8 - GSE73680 single-gene and module results
9. Supplementary Table 9 - Random-set benchmark results
10. Supplementary Table 10 - Functional enrichment results
11. Supplementary Table 11 - Risk-injury coupling results
12. Supplementary Table 12 - TWAS/SMR/coloc/spatial resource audit
13. Supplementary Table 13 - Claim-boundary audit

Detailed GWAS row-level QC, Agilent file reconstruction counts, expression-matched benchmarking implementation, file manifests and external-resource status should be moved from the main Methods into the corresponding supplementary methods notes and tables.
""")
print(OUT)
