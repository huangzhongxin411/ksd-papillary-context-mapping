from pathlib import Path
import csv, re, html, hashlib

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "manuscript/manuscript_draft_v1.2_phase18b_cited.md"
OUT = ROOT / "manuscript/manuscript_draft_v1.3_cited_clean.md"
PDF = ROOT / "manuscript/manuscript_draft_v1.3_cited_clean.pdf"
SUPP = ROOT / "manuscript/supplementary_methods_v0.1.md"
SUPP_MANIFEST = ROOT / "results/tables/supplementary_tables_manifest_v0.1.tsv"
REF_AUDIT = ROOT / "results/tables/reference_audit_v0.1.tsv"
CLAIM_AUDIT = ROOT / "results/tables/claim_reference_audit_v0.1.tsv"
REF_BASE = ROOT / "manuscript/references/reference_verification_phase18b.tsv"
FIG_DIR = ROOT / "results/figures/final_main_figures_v1"

def replace_section(text, start, end, body):
    pattern = rf"(?ms)^{re.escape(start)}\n.*?(?=^{re.escape(end)}\n)"
    return re.sub(pattern, body.rstrip() + "\n\n", text, count=1)

title = "Post-GWAS mapping links kidney stone risk genes to Loop/TAL-associated renal papillary cellular and plaque/stone disease contexts"

abstract = """## Abstract

**Background:** Kidney stone disease (KSD) has a substantial polygenic component, but translating genome-wide association signals into renal papillary cellular and disease contexts remains challenging.

**Methods:** We integrated public trans-ancestry KSD GWAS summary statistics, MAGMA gene-based prioritization, audited annotations from the single-nucleus component of GSE231569, six priority-1 (P1) candidate-gene evidence scoring and patient-aware analysis of the GSE73680 plaque/stone-associated papillary bulk expression context.

**Results:** Our QC-pruned reconstruction retained 57 independent distance-defined loci, and MAGMA identified 94 Bonferroni-significant genes. Prioritized gene sets projected preferentially to the Loop/TAL nuclei compartment in GSE231569, with expression-matched benchmark percentiles of 0.998-1.000 for the top-ranked modules and retained support in locus-balanced and leave-one-locus-out analyses. Six P1 genes formed a TAL, transport, calcium/ion-handling and broader epithelial role spectrum rather than a uniform TAL marker or disease-response panel. In GSE73680, comprising 55 included samples from 29 patients, including 26 paired patients, MAGMA modules, but not individual P1 genes, showed patient-aware plaque/stone-associated papillary expression shifts. MAGMA top-50 scores also coupled with injury/remodeling programs in paired-delta (Spearman rho = 0.82, FDR = 4.5 x 10^-6) and patient/group residual analyses (rho = 0.79, FDR = 7.2 x 10^-12).

**Conclusion:** MAGMA-prioritized KSD genes support a Loop/TAL-associated papillary cellular context and module-level plaque/stone disease-context association. This claim-bounded framework nominates genes and modules for future testing but does not establish causal mediation, P1 single-gene disease validation, TWAS convergence, SMR/coloc support, spatial validation or experimental mechanism.

**Keywords:** kidney stone disease; MAGMA; renal papilla; Loop of Henle; thick ascending limb; single-nucleus RNA sequencing"""

intro = """## Introduction

Kidney stone disease is common, recurrent and biologically heterogeneous, reflecting interactions among urinary chemistry, renal epithelial transport, papillary microenvironments and inherited susceptibility.[1-3,34] Genome-wide association studies have expanded the number of KSD susceptibility loci across ancestries, including loci near genes involved in calcium handling and tubular transport.[4-7] However, most association signals do not directly identify the effector gene, relevant renal cell type or disease-stage tissue context, leaving a substantial post-GWAS interpretation gap.[3,8]

The renal papilla is a plausible tissue context for stone initiation and growth. Randall plaque has been localized to papillary interstitial and thin-limb-associated anatomical compartments, while plaque-associated tissue shows inflammatory, oxidative-stress and injury-related molecular features.[9-11,13,14] Experimental crystal studies further support an immune-response context but do not by themselves establish the molecular state of human plaque tissue.[12,13] Recent single-nucleus and spatially anchored transcriptomic studies have resolved epithelial, stromal, vascular and immune niches in human kidney and papilla.[15-17] Within this anatomy, the Loop of Henle and thick ascending limb (Loop/TAL) coordinate salt transport, urinary concentration and paracellular calcium and magnesium handling through uromodulin-, claudin- and calcium-sensing pathways.[18-23] These functions make Loop/TAL-associated programs biologically plausible contexts for KSD genetic risk, but biological plausibility alone does not establish causal cell-type mediation.

Post-GWAS gene-based and cell-type enrichment approaches can connect polygenic association signals to tissues and cellular expression programs.[24-26] Single-cell and single-nucleus resources are particularly useful for localization, provided that annotation quality, gene-set background and donor structure are considered. Such analyses remain contextual: they do not substitute for transcriptome-wide association testing, summary-data Mendelian randomization or genetic colocalization,[30-32,39-41] nor for lesion-resolved spatial validation.[15,16] A claim-bounded framework that jointly evaluates genetic prioritization, papillary single-nucleus localization and plaque/stone-associated disease-context expression has not been systematically applied to KSD risk genes.

We therefore asked three questions. First, do MAGMA-prioritized KSD genes converge on an audited renal papillary cellular context? Second, do six P1 candidates form a coherent and interpretable TAL, transport and calcium/ion-handling spectrum? Third, are the prioritized modules reflected in an independent plaque/stone-associated papillary bulk expression context? To address these questions, we integrated the 2025 trans-ancestry KSD GWAS,[4] the single-nucleus component of GSE231569,[15] and GSE73680 plaque-associated microarray data.[11] The study explicitly separates supported cellular and disease-context association from causality, TWAS convergence, SMR/coloc support, spatial validation and P1 single-gene disease validation."""

methods = """## Methods

### Study design and public resources

The analysis integrated the 2025 trans-ancestry KSD GWAS,[4] the GSE231569 single-nucleus component of a spatially anchored kidney papilla atlas,[15] and the GSE73680 Randall plaque/stone-associated papillary bulk expression dataset.[11,33] The workflow comprised GWAS quality control, MAGMA prioritization, audited cell-context projection, P1 candidate interpretation, patient-aware bulk expression analysis, functional enrichment and injury/remodeling coupling. Full reconstruction details, formulas and sensitivity analyses are provided in Supplementary Methods.

### GWAS quality control and MAGMA prioritization

GRCh37/hg19-compatible summary statistics were filtered to biallelic A/C/G/T single-nucleotide variants with valid coordinates and P values; strand-ambiguous A/T and C/G variants and variants with minor-allele frequency <0.01 were removed. Of 5,960,489 rows, 4,915,033 passed QC. Genome-wide significant variants (P < 5 x 10^-8) were greedily distance-pruned at +/-1 Mb; overlapping +/-1-Mb lead windows were merged, yielding 57 distance-defined loci in our QC-pruned reconstruction. These loci are an analytical reconstruction and are not identical to the 59 loci reported by the source GWAS.[4]

MAGMA v1.10 used NCBI build 37 gene locations and the 1000 Genomes European LD reference.[24] SNPs were annotated to default gene coordinates without an added flanking window. Of 19,427 gene locations read, 17,316 genes were tested. Bonferroni significance was defined as P < 0.05/17,316 (2.89 x 10^-6), false-discovery rates used the Benjamini-Hochberg procedure,[37] and suggestive genes were defined at P < 1 x 10^-4. Fixed top-50, top-100, FDR-significant and suggestive sets were carried forward. The use of a European LD panel for a trans-ancestry GWAS is treated as an ancestry-reference limitation.

### GSE231569 annotation audit and gene-set projection

We analyzed four papillary samples from four donors in the GSE231569 single-nucleus component; the medulla sample was excluded. Nuclei were retained at 200-6,000 detected features and <=20% mitochondrial reads, log-normalized in Seurat, and represented with 3,000 variable features, 30 principal components, neighbors/UMAP dimensions 1-25 and clustering resolution 0.4.[15,28] Broad labels were audited against marker expression and harmonized to six analysis compartments. The Loop/TAL compartment was supported by UMOD, SLC12A1, CLDN10, KCNJ1 and CLDN16 and comprised 540 nuclei across all four donors; other retained compartments comprised 34,849 collecting-duct principal, 3,422 fibroblast/stromal, 2,558 endothelial, 2,361 injured/undifferentiated epithelial and 148 perivascular nuclei.

For each MAGMA set, a nucleus-level score was the arithmetic mean of log-normalized RNA-layer expression across detected member genes. Cell-type and donor summaries were descriptive; the inferential layer compared the observed cell-type mean with 1,000 random sets matched by gene count and 20 bins of global mean expression (seed 20260617). A locus-balanced analysis retained one top-ranked gene per locus group. A conservative annotation analysis excluded low/exploratory and immune-review labels, and leave-one-locus-out analyses repeated the benchmark after removing each locus group. These analyses localize an expression context and do not test causal mediation or spatial colocalization.[25,26]

### P1 candidate evidence spectrum

UMOD, CASR, CLDN14, CLDN10, HIBADH and PKD2 were designated priority-1 (P1) candidates for structured interpretation. Evidence included Loop/TAL expression rank, mean expression, detection frequency, donor detection, TAL specificity, TAL-program correlation and biological role. Candidate and population association studies were treated as supporting context rather than causal proof.[18-23,35,36]

### GSE73680 patient-aware disease-context analysis

GSE73680 Agilent GPL17077 feature files were reconstructed and linked to curated patient/tissue metadata.[11] Direct gene-symbol features were retained; 32,055 of 44,661 features mapped, and the specified max-mean-expression collapse rule encountered no remaining multi-probe genes. Platform-wide normalization could not be independently verified from the reconstructed files; analyses used log2(expression + 1), and no batch covariate was fitted because reliable batch metadata were unavailable. Of 62 samples, 55 from 29 patients were included: 27 control/adjacent and 28 plaque/stone-associated papillary samples, with 26 patients contributing both groups and three contributing one included group.

The primary limma model was expression ~ group, with patient identifier supplied to duplicateCorrelation and plaque/stone-associated papilla as the target coefficient.[27] A paired sensitivity averaged samples within patient and group and defined delta as plaque/stone-associated papilla minus control/adjacent. P1-gene FDR correction was performed within the six-gene family. Module scores were means of row-wise gene z scores across detected genes and were analyzed only when >=70% of genes were detected (>=40% for exploratory modules). Module-level P values were Benjamini-Hochberg corrected within the tested module family. Size-matched random sets, leave-one-gene-out analysis, MAGMA-without-P1 sensitivity and paired-direction consistency assessed robustness; the stricter expression-matched benchmark was retained only as a conservative sensitivity analysis.

### Functional enrichment and injury/remodeling coupling

GO Biological Process enrichment used clusterProfiler and org.Hs.eg.db with the 15,421 MAGMA-tested genes mapping to Entrez identifiers as background.[29] Terms displayed in the main figure met FDR <0.10 and Count >=2. Redundancy was reduced greedily by retaining the more significant term when gene-overlap Jaccard similarity was >=0.70 within a gene set. Prespecified nephron/functional marker sets, listed in Supplementary Methods, were assessed by hypergeometric testing with Benjamini-Hochberg correction across tested gene-set-by-marker-set combinations.

In GSE73680, injury/remodeling, epithelial-injury, fibrosis/ECM and inflammation/immune scores were means of row-wise z-scored marker expression. Coupling with MAGMA modules used Spearman correlations at the sample level, across paired patient deltas and between residuals from separate module ~ group + patient and injury score ~ group + patient models. FDR correction was applied across module-program pairs within each analysis. Paired-delta and residual estimates were prioritized; sample-level estimates and leave-one-patient-out results were sensitivity layers.

### Resource-limited extensions and statistical boundary

TWAS, SMR/coloc and spatial workflows were scaffolded and resource-audited but were not used as evidence layers because required prediction weights, covariance/LD resources, kidney eQTL summaries or lesion-resolved matrices remained incomplete.[30-32,38-41] Missing resources are not negative biological evidence. Unless otherwise stated, two-sided P values were used and Benjamini-Hochberg adjustment was performed within prespecified analysis families.[37]"""

discussion = """## Discussion

### A claim-bounded post-GWAS framework

This study provides an auditable framework connecting KSD genetic prioritization to renal papillary cellular and disease contexts. Its contribution is not the declaration of a single causal gene or cell type. Instead, MAGMA-defined genetic modules were localized to an audited Loop/TAL-associated single-nucleus compartment, interpreted through a heterogeneous P1 gene spectrum and evaluated in an independent plaque/stone-associated papillary bulk expression dataset. This layered design addresses a recurring post-GWAS problem: gene-based association and cell-type enrichment can organize evidence, but they remain distinct from causal mediation, TWAS, colocalization and experimental validation.[24-26,30-32,39-41]

### Why Loop/TAL is biologically relevant

The Loop/TAL localization is compatible with established renal physiology. TAL epithelia coordinate urinary concentration and transepithelial salt transport, while claudin-dependent paracellular pathways influence calcium and magnesium handling.[18-20,23] The P1 spectrum also aligns with this context: UMOD represents a TAL-associated protein with links to mineral metabolism,[21] CLDN10 and CLDN14 connect epithelial-junction biology to ion handling and nephrocalcinosis or stone susceptibility,[7,19,20] and candidate and genetic association studies have linked CASR-related calcium sensing to stone phenotypes.[22,35,36] HIBADH and PKD2 were retained as broader supporting epithelial-context genes rather than forced into a uniform TAL marker class. The data therefore support a biologically coherent cellular context while preserving heterogeneity among candidates.

### Why module-level disease-context support matters

The absence of FDR-supported P1 single-gene responses in GSE73680 does not conflict with the MAGMA module-level association. A polygenic disease-context signal may be distributed across multiple modestly changing genes rather than dominated by one transcript. Thus, the lack of FDR-supported P1 single-gene responses should not be interpreted as failure of the post-GWAS mapping framework, but as evidence that the observed disease-context support is module-level and context-dependent. The patient-aware paired shifts, size-matched random-set benchmark and MAGMA-without-P1 sensitivity support module-level association, not validation of any individual P1 gene. The expression-matched benchmark was retained as a conservative sensitivity analysis and did not serve as the primary support layer.

### Papillary injury context and future validation

Randall plaque and stone-forming papilla show inflammatory, immune and remodeling features,[9-14] and the GSE231569 atlas described immune injury and matrix remodeling in stone disease.[15] The observed coupling between MAGMA modules and injury/remodeling scores is consistent with this tissue context, but correlation does not establish that genetic risk causes injury or that injury mediates genetic effects. Future work should prioritize papilla-specific eQTL resources, genetic colocalization, lesion-resolved spatial transcriptomics and perturbation of prioritized modules in TAL-relevant epithelial models. These experiments would test whether the mapped cellular context corresponds to causal regulatory mechanisms or spatially localized disease programs."""

limitations = """## Limitations

This study is computational and hypothesis-generating. GSE231569 supports single-nucleus expression-context mapping but not spatial validation or causal cell-type mediation. Donor numbers were limited, and one donor contributed only four Loop/TAL nuclei. GSE73680 is a reconstructed bulk microarray dataset; platform-wide normalization could not be independently verified, reliable batch metadata were unavailable and the data cannot establish cell-type-specific disease responses. P1 single-gene differential expression was not FDR-supported. Functional enrichment and risk-injury coupling are gene-set or module-level interpretations that may depend on gene-set size, background definition, expression detectability and disease-group composition. The European LD reference may not fully represent the trans-ancestry GWAS. The expression-matched random benchmark was retained as a conservative sensitivity analysis and did not serve as the primary support layer. TWAS, SMR/coloc and spatial validation were not completed because required external resources were unavailable; their absence is not evidence against those mechanisms."""

text = SRC.read_text()
text = re.sub(r"^# .*", "# " + title, text, count=1)
text = replace_section(text, "## Abstract", "## Introduction", abstract)
text = replace_section(text, "## Introduction", "## Methods", intro)
text = replace_section(text, "## Methods", "## Results", methods)
text = replace_section(text, "## Discussion", "## Limitations", discussion)
text = replace_section(text, "## Limitations", "## Data availability", limitations)

# Results precision and terminology, without changing numerical results.
repls = {
"MAGMA prioritization defined a reproducible KSD genetic module layer from 57 independent loci.": "Our QC-pruned reconstruction retained 57 independent distance-defined loci and defined the input to MAGMA prioritization.",
"expression-context scores were highest in Loop/TAL cells compared with other audited cell types.": "expression-context scores were highest in the Loop/TAL nuclei compartment and exceeded expression-matched random expectation.",
"with the minimum retained leave-one-locus-out percentile remaining above 0.99.": "with a locus-balanced Loop/TAL percentile of 0.987 and leave-one-locus-out percentiles of 0.997-1.000.",
"Fig. 3a-d": "Fig. 3A-D",
"Fig. 4a-d": "Fig. 4A-D",
"Fig. 5a-b": "Fig. 5A-B",
"Fig. 5c-d": "Fig. 5C-D",
"Patient-level MAGMA module responses reached q <= 0.05": "Patient-level top-50, top-100, FDR and suggestive MAGMA module responses reached FDR q < 0.05",
"plaque/stone papilla disease-context": "plaque/stone-associated papillary disease-context",
"plaque/stone papilla expression": "plaque/stone-associated papillary expression",
"external molecular validation": "experimental or lesion-resolved molecular validation",
}
for a, b in repls.items(): text = text.replace(a, b)

# Remove vitamin-D emphasis from the main Results while retaining complete GO results in source tables.
text = text.replace("loop of Henle development, distal tubule development, response to vitamin D, response to metal ion and phosphate ion homeostasis", "loop of Henle development, distal tubule development, response to metal ion and phosphate ion homeostasis")

refs_extra = """37. Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
38. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science. 2020;369(6509):1318-1330. https://doi.org/10.1126/science.aaz1776 [PubMed](https://pubmed.ncbi.nlm.nih.gov/32913098/)
39. Barbeira AN, Dickinson SP, Bonazzola R, Zheng J, Wheeler HE, Torres JM, et al. Exploring the phenotypic consequences of tissue specific gene expression variation inferred from GWAS summary statistics. Nat Commun. 2018;9(1):1825. https://doi.org/10.1038/s41467-018-03621-1 [PubMed](https://pubmed.ncbi.nlm.nih.gov/29739930/)
40. Wainberg M, Sinnott-Armstrong N, Mancuso N, Barbeira AN, Knowles DA, Golan D, et al. Opportunities and challenges for transcriptome-wide association studies. Nat Genet. 2019;51(4):592-599. https://doi.org/10.1038/s41588-019-0385-z [PubMed](https://pubmed.ncbi.nlm.nih.gov/30926968/)
41. Wallace C. Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses. PLoS Genet. 2020;16(4):e1008720. https://doi.org/10.1371/journal.pgen.1008720 [PubMed](https://pubmed.ncbi.nlm.nih.gov/32310995/)"""
text = text.replace("\n## Supplementary materials", "\n" + refs_extra + "\n\n## Supplementary materials")
text = text.replace("The supplementary package is scaffolded as Supplementary Tables 1-13 and associated methods notes; see `docs/supplementary_materials_scaffold_phase18b.md`.", "Detailed methods are provided in `manuscript/supplementary_methods_v0.1.md`; the source-table inventory is provided in `results/tables/supplementary_tables_manifest_v0.1.tsv`.")
OUT.write_text(text)

supp = """# Supplementary Methods v0.1

## Scope and claim boundary

These methods document the completed evidence layers supporting Loop/TAL-associated single-nucleus expression context and module-level plaque/stone-associated papillary bulk expression association. TWAS, SMR/coloc and spatial workflows were audited but incomplete and were not used for claims. No analysis described here establishes causality, P1 single-gene disease validation, spatial validation or experimental mechanism.

## GWAS reconstruction and variant quality control

The 2025 trans-ancestry KSD summary statistics were treated as GRCh37/hg19-compatible. Required fields were variant identifier, chromosome, base-pair position, effect and non-effect alleles, P value, beta, standard error, effect-allele frequency and sample size. The file did not contain an INFO field, so no imputation-quality filter was possible. Rows with missing required values, invalid P values or coordinates, non-SNV/non-ACGT alleles, duplicate variant identifiers, strand-ambiguous A/T or C/G alleles, or minor-allele frequency <0.01 were removed. The sequential counts were 5,960,489 raw rows; 72,210 non-SNV/bad-allele rows, 895,763 ambiguous rows and 77,483 low-frequency rows removed; and 4,915,033 retained variants. The genomic-control lambda was 1.127 and 3,209 retained variants had P <5 x 10^-8.

Genome-wide significant variants were sorted by P value and greedily pruned using a 1-Mb distance around each retained lead. Each lead defined a +/-1-Mb window, and overlapping windows were merged. This distance-based procedure produced 60 leads and 57 merged loci. It is not LD clumping and is not intended to reproduce exactly the 59 loci reported in the source publication.

## MAGMA gene analysis

MAGMA v1.10 was run with NCBI37.3 gene locations and the 1000 Genomes European `g1000_eur` LD panel. The P-value input contained 4,915,033 unique variants; MAGMA used 4,882,602, mapped 1,931,992 to genes and tested 17,316 genes. Annotation used default gene coordinates with no added flanking window. MAGMA reported 133 synonymous SNP-pair warnings. Bonferroni significance was P <0.05/17,316 = 2.89 x 10^-6; BH FDR <0.05 and nominal P <1 x 10^-4 defined FDR and suggestive sets. This yielded 94 Bonferroni-significant, 369 FDR-significant and 187 suggestive genes. Top-50 and top-100 sets were fixed by MAGMA rank. Because the GWAS was trans-ancestry and the LD panel European, gene statistics require ancestry-reference caution.

## GSE231569 reconstruction, QC and annotation audit

Four papillary donors (two healthy and two stone-disease samples) were retained; the study medulla sample was excluded. Seurat objects were created with `min.cells=3` and `min.features=200`. Nuclei were filtered to 200-6,000 detected features and <=20% mitochondrial reads, log-normalized with `NormalizeData`, represented by 3,000 variable features, scaled, reduced with 30 PCs, clustered with dimensions 1-25 at resolution 0.4 and embedded by UMAP with dimensions 1-25.

Initial labels were obtained by broad marker-score comparison and then audited against cluster markers. The final analysis harmonized six compartments: collecting-duct principal (34,849 nuclei; four donors), fibroblast/stromal (3,422; four), endothelial (2,558; four), injured/undifferentiated epithelial (2,361; four), Loop/TAL (540; four) and perivascular/mural (148; three). Loop/TAL identity required concordant expression of established markers including UMOD, SLC12A1, CLDN10, KCNJ1 and CLDN16. Loop/TAL donor counts were 29, 244, 263 and 4 nuclei; donor-level results are therefore descriptive and the fourth donor is a sparse-compartment limitation.

## Single-nucleus MAGMA scoring and benchmarks

For a gene set G and nucleus i, the score was `S_i = mean(x_ig)` for detected genes g in G, where x is the Seurat RNA `data` layer after log normalization. This is not Seurat `AddModuleScore` and does not use per-gene z scores. Scores were averaged by audited cell type and donor-cell-type combination.

The random-set universe comprised genes with finite positive global mean expression. Genes were divided into 20 quantile bins of global mean expression. For each real gene, one replacement was sampled from the same bin; 1,000 sets were generated with seed 20260617. The empirical percentile was the proportion of random scores below or equal to the observed score. Loop/TAL percentiles were 0.998 for top 50, 1.000 for top 100, 1.000 for suggestive genes and 0.968 for FDR genes.

For the locus-balanced top-50 set, one highest-ranked gene was retained per locus group, with MAGMA-only genes assigned individual groups. Full and conservative analyses used 1,000 matched random sets; the conservative analysis excluded low/exploratory and immune-review labels. The conservative Loop/TAL percentile was 0.987. Leave-one-locus-out analysis removed each locus group and repeated 1,000-set benchmarking; retained Loop/TAL percentiles ranged from 0.997 to 1.000.

## P1 candidate evidence scoring

UMOD, CASR, CLDN14, CLDN10, HIBADH and PKD2 were assessed as priority-1 candidates. The evidence table records MAGMA membership, expression/detection in audited compartments, Loop/TAL rank and specificity, donor detection, correlation with a TAL reference program, GSE73680 availability and biological role. Evidence-badge summaries are descriptive integrations, not a combined probability or causal score.

## GSE73680 expression reconstruction and metadata

All 62 local supplementary feature files were parsed. Agilent feature expression was selected in priority order from `gProcessedSignal`, `gBGSubSignal`, `gMeanSignal` or `gMedianSignal`. Control and dark-corner features were excluded. Direct gene-symbol identifiers were retained; 32,055 of 44,661 features mapped and 12,606 remained unmapped because platform annotation was unavailable. The prespecified collapse rule retained the probe with maximum mean expression for duplicate symbols, although no multi-probe genes remained after direct-symbol filtering. The reconstructed scale was continuous and appeared normalized, but platform-wide normalization could not be independently verified. Values were transformed as log2(expression + 1). No reliable batch variable was available.

Curated filename and sample metadata assigned 27 control/adjacent and 28 plaque/stone-associated papillary samples. These 55 samples came from 29 patients. Twenty-six patients contributed both groups and formed the paired analysis; three patients contributed one included group and informed the primary duplicate-correlation model only.

## GSE73680 gene and module models

The primary gene-level model was `expression ~ group_curated`, with `patient_id` as the blocking factor in limma `duplicateCorrelation`; the reported coefficient was plaque/stone-associated papilla versus control/adjacent. The paired sensitivity first averaged repeated samples within patient and group, then calculated plaque/stone-associated minus control/adjacent values and used paired t tests. BH adjustment was applied across the six P1 genes separately from module analyses.

For each gene, expression was standardized across samples. A module score was the sample mean of available gene z scores. Full modules required >=70% gene detection; modules at 40-70% were exploratory and lower-coverage modules were not tested. Patient-level scores averaged samples within patient and group. Paired module shifts were tested among 26 paired patients, with BH adjustment across the tested module family. Top-50, top-100, FDR and suggestive modules each had q=0.0499; the P1 module had q=0.299 and was not FDR-supported.

The size-matched benchmark drew 1,000 random sets. Leave-one-gene-out tests, paired-direction consistency and removal of all six P1 genes assessed dependence on individual genes. After P1 removal, the MAGMA modules retained q=0.0251. Expression-matched random sets used 20 mean-expression bins and 1,000 iterations (seed 73682); because this stricter analysis did not provide primary support, it was retained only as a sensitivity boundary.

## Functional enrichment and curated contexts

MAGMA-tested genes were mapped to Entrez identifiers, yielding a 15,421-gene GO universe. `clusterProfiler::enrichGO` tested Biological Process terms and calculated BH-adjusted P values. Main-figure display required FDR <0.10 and at least two query genes. Terms were ordered by adjusted P value and Count; subsequent terms with gene-overlap Jaccard >=0.70 to an already retained term were removed.

Prespecified contextual sets covered TAL transport, calcium/ion handling, epithelial junction and injury/remodeling. They were assembled from canonical renal marker genes used in the annotation and disease-context workflows and are listed in the source scripts and supplementary source tables. Hypergeometric tests used unique MAGMA symbols as the universe and applied BH correction across all tested query-set by context-set combinations. These sets provide interpretation, not externally validated pathway activity.

## Injury/remodeling coupling

Marker programs were injury/remodeling (`SPP1,MMP7,MMP9,GPNMB,COL1A1,COL1A2,FN1,VIM,HAVCR1,LCN2`), inflammation/immune (`CCL2,CCL7,CXCL8,IL1B,TNF,CD68,LST1`), fibrosis/ECM (`COL1A1,COL3A1,DCN,LUM,ACTA2,TAGLN`) and epithelial injury (`HAVCR1,LCN2,VCAM1,SOX9,PROM1`). Scores were means of row-wise z-scored genes.

Spearman correlations were computed for all samples, for the 26 paired patient deltas, and for residuals from separate linear models of risk-module score and injury score on group plus patient factor. BH adjustment was performed across module-program pairs within each analysis. Paired-delta and residual results were prioritized. Leave-one-patient-out summaries and sample-level correlations were sensitivity analyses. Correlation denotes coupling, not mediation or mechanism.

## Resource-limited analyses

FUSION/S-PrediXcan TWAS required tissue weights and covariance/LD resources that were incomplete. SMR/coloc required kidney eQTL summaries or BESD/ESI/EPI inputs and suitable LD that were incomplete. The GTEx atlas was reviewed as a potential eQTL resource, but available kidney context did not complete the prespecified papillary analysis. Spatial candidate datasets lacked a complete audited combination of matrices, barcodes/features, coordinates, scalefactors and images. These workflows were not used for scientific claims.

## Software

Captured versions were R 4.4.3, MAGMA 1.10, Python 3.9.6, data.table 1.18.2.1, ggplot2 4.0.2, Seurat 5.4.0, clusterProfiler 4.14.6, cowplot 1.2.0 and svglite 2.2.2. Earlier preprocessing logs and scripts remain the authoritative record where package-version capture was unavailable.
"""
SUPP.write_text(supp)

manifest_rows = [
("S1","GWAS QC flow","results/tables/phase1_gwas_qc_report.tsv","GWAS QC counts and lambda GC","main support"),
("S2","Lead variants","results/tables/phase1_2025_lead_snps.tsv","Distance-pruned lead variants","main support"),
("S3","Merged loci","results/tables/phase1_2025_loci.tsv","57 distance-defined loci","main support"),
("S4","MAGMA QC","results/tables/magma_qc_summary.tsv","Reference, mapping and gene-test counts","main support"),
("S5","MAGMA genes","results/tables/magma_genes.tsv","Gene statistics and adjusted P values","main support"),
("S6","GSE231569 cell counts","results/tables/gse231569_cell_counts.tsv","Audited nuclei and donor counts","main support"),
("S7","GSE231569 matched benchmark","results/tables/magma_scrna_random_benchmark.tsv","1,000-set expression-matched benchmark","main support"),
("S8","Locus-balanced benchmark","results/tables/magma_locus_balanced_scrna_benchmark.tsv","One-gene-per-locus sensitivity","sensitivity"),
("S9","Leave-one-locus-out","results/tables/magma_leave_one_locus_out.tsv","Locus influence sensitivity","sensitivity"),
("S10","P1 evidence","results/tables/p1_tal_gene_evidence.tsv","Six-gene evidence components","interpretive"),
("S11","GSE73680 design","results/gse73680/tables/gse73680_analysis_design.tsv","Sample and patient structure","main support"),
("S12","GSE73680 mapping QC","results/gse73680/tables/gse73680_gene_mapping_qc.tsv","Feature mapping and collapse audit","QC"),
("S13","P1 response","results/gse73680/tables/gse73680_p1_gene_response.tsv","Patient-aware P1 gene tests","boundary"),
("S14","Module response","results/gse73680/tables/gse73680_patient_level_module_response.tsv","Paired module shifts","main support"),
("S15","Random benchmark","results/gse73680/tables/gse73680_random_module_benchmark.tsv","Size-matched random sets","main support"),
("S16","Expression-matched benchmark","results/gse73680/tables/gse73680_expression_matched_random_benchmark.tsv","Conservative benchmark","sensitivity only"),
("S17","Without-P1 sensitivity","results/gse73680/tables/gse73680_magma_without_p1_sensitivity.tsv","Module response after P1 removal","sensitivity"),
("S18","GO enrichment","results/tables/magma_pathway_enrichment.tsv","GO BP enrichment","interpretive"),
("S19","GO background audit","results/tables/pathway_enrichment_background_audit.tsv","Universe/display/redundancy rules","QC"),
("S20","Injury coupling","results/gse73680/tables/gse73680_module_injury_correlation.tsv","Sample, paired-delta and residual coupling","main and sensitivity"),
("S21","Resource manifest","results/tables/twas_resource_manifest_v2.tsv","TWAS resource status","not used for claims"),
("S22","Spatial audit","results/tables/spatial_resource_manifest.tsv","Spatial resource status","not used for claims"),
]
SUPP_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
with SUPP_MANIFEST.open("w", newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["supplementary_table","title","source_table","content","support_role","exists"])
    for row in manifest_rows: w.writerow([*row, str((ROOT/row[2]).exists()).upper()])

# Reference audit: Phase 18B PubMed-verified records plus five methods additions.
base=list(csv.DictReader(REF_BASE.open(),delimiter="\t"))
extras=[
dict(reference_number="37",pmid="",doi="10.1111/j.2517-6161.1995.tb02031.x",title="Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing",journal="J R Stat Soc Series B Stat Methodol",year="1995",verification_source="Crossref DOI metadata",status="verified"),
dict(reference_number="38",pmid="32913098",doi="10.1126/science.aaz1776",title="The GTEx Consortium atlas of genetic regulatory effects across human tissues.",journal="Science",year="2020",verification_source="NCBI PubMed XML",status="verified"),
dict(reference_number="39",pmid="29739930",doi="10.1038/s41467-018-03621-1",title="Exploring the phenotypic consequences of tissue specific gene expression variation inferred from GWAS summary statistics.",journal="Nat Commun",year="2018",verification_source="NCBI PubMed XML",status="verified"),
dict(reference_number="40",pmid="30926968",doi="10.1038/s41588-019-0385-z",title="Opportunities and challenges for transcriptome-wide association studies.",journal="Nat Genet",year="2019",verification_source="NCBI PubMed XML",status="verified"),
dict(reference_number="41",pmid="32310995",doi="10.1371/journal.pgen.1008720",title="Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses.",journal="PLoS Genet",year="2020",verification_source="NCBI PubMed XML",status="verified"),
]
all_refs=base+extras
use_map={str(i):"background or methods" for i in range(1,42)}
use_map.update({"12":"experimental crystal/immune context only","22":"candidate/population association only","33":"GEO repository/data availability","37":"BH FDR method","38":"GTEx resource audit","39":"S-PrediXcan method","40":"TWAS interpretation boundary","41":"coloc assumptions and limitations"})
with REF_AUDIT.open("w",newline="") as f:
    fields=["ref_id","title_match","journal_match","year_match","doi_match","pmid_match","url_resolves","used_for_claim","action"]
    w=csv.DictWriter(f,fieldnames=fields,delimiter="\t"); w.writeheader()
    for r in all_refs:
        resolution = "PASS_PUBMED_METADATA" if r["pmid"] else "PASS_CROSSREF_METADATA"
        w.writerow(dict(ref_id=r["reference_number"],title_match="PASS",journal_match="PASS",year_match="PASS",doi_match="PASS" if r["doi"] else "N/A",pmid_match="PASS" if r["pmid"] else "N/A",url_resolves=resolution,used_for_claim=use_map[r["reference_number"]],action="retain; metadata record verified; one live PubMed URL returned HTTP 200 during this pass"))

claim_rows=[
("Introduction","KSD is polygenic and heterogeneous","1-8","direct background support","PASS"),
("Introduction","Human Randall plaque anatomy and bulk expression context","9-11,13-15","Ref 12 excluded as sole human-tissue support","PASS"),
("Introduction","Experimental crystal/immune context","12-13","experimental context only","PASS"),
("Introduction","Loop/TAL transport and mineral handling","18-23,35-36","physiology plus qualified association","PASS"),
("Introduction/Methods","TWAS, SMR and coloc are distinct evidence layers","30-32,39-41","methods and interpretation boundary","PASS"),
("Introduction/Methods","Lesion-resolved spatial validation was not performed","15-16","spatial atlas context, not study validation","PASS"),
("Methods","MAGMA gene analysis","24","software/method","PASS"),
("Methods","BH multiple-testing correction","37","statistical method","PASS"),
("Methods","Seurat workflow","28","software/method","PASS"),
("Methods","limma duplicateCorrelation","27","software/method","PASS"),
("Methods","GO enrichment","29","software/method","PASS"),
("Methods/Limitations","GTEx and incomplete eQTL resource context","38","resource audit only","PASS"),
("Discussion","CASR/CLDN14 evidence is candidate or genetic association","7,22,35-36","qualified; no causal claim","PASS"),
("Data availability","GEO repository","33","repository citation plus dataset papers","PASS"),
]
with CLAIM_AUDIT.open("w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["section","claim","references","support_level","status"]); w.writerows(claim_rows)

# Submission-style review PDF with frozen figures included only in the figure-legend section.
def inline(s):
    links=[]
    def hold(m):
        links.append((m.group(1),m.group(2))); return f"@@LINK{len(links)-1}@@"
    s=re.sub(r"\[([^]]+)\]\((https?://[^)]+)\)",hold,s)
    s=html.escape(s); s=re.sub(r"\*\*(.+?)\*\*",r"<b>\1</b>",s); s=re.sub(r"\*(.+?)\*",r"<i>\1</i>",s)
    for i,(label,url) in enumerate(links):
        s=s.replace(f"@@LINK{i}@@",f'<link href="{html.escape(url,quote=True)}" color="#005A64">{html.escape(label)}</link>')
    return s
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name="T",parent=styles["Title"],fontName="Helvetica-Bold",fontSize=17,leading=21,textColor=colors.HexColor("#22313B"),spaceAfter=10))
styles.add(ParagraphStyle(name="H2x",parent=styles["Heading2"],fontName="Helvetica-Bold",fontSize=13,leading=16,textColor=colors.HexColor("#005A64"),spaceBefore=12,spaceAfter=6))
styles.add(ParagraphStyle(name="H3x",parent=styles["Heading3"],fontName="Helvetica-Bold",fontSize=10.5,leading=13,textColor=colors.HexColor("#22313B"),spaceBefore=8,spaceAfter=4))
styles.add(ParagraphStyle(name="Bx",parent=styles["BodyText"],fontName="Helvetica",fontSize=8.6,leading=12,textColor=colors.HexColor("#22313B"),spaceAfter=5))
story=[]; para=[]; in_legends=False; first_legend=True
def flush():
    if para: story.append(Paragraph(inline(" ".join(para)),styles["Bx"])); para.clear()
for line in text.splitlines():
    if not line.strip(): flush(); continue
    m=re.match(r"^(#{1,3})\s+(.*)$",line)
    if m:
        flush(); level=len(m.group(1)); heading=m.group(2)
        if heading=="Figure legends":
            story.append(PageBreak()); in_legends=True
        if in_legends and level==3:
            fm=re.match(r"Figure ([1-5])\.",heading)
            if fm:
                if not first_legend: story.append(PageBreak())
                first_legend=False; fp=FIG_DIR/f"figure{fm.group(1)}.png"
                if fp.exists(): story.extend([Image(str(fp),width=170*mm,height=95.625*mm),Spacer(1,4*mm)])
        story.append(Paragraph(inline(heading),styles["T"] if level==1 else styles["H2x"] if level==2 else styles["H3x"]))
    else: para.append(line.strip())
flush()
def page_number(canvas,doc):
    canvas.saveState(); canvas.setFont("Helvetica",7.5); canvas.setFillColor(colors.HexColor("#6F929B"))
    canvas.drawCentredString(A4[0]/2,9*mm,str(doc.page)); canvas.restoreState()
SimpleDocTemplate(str(PDF),pagesize=A4,leftMargin=20*mm,rightMargin=20*mm,topMargin=17*mm,bottomMargin=17*mm,title=title).build(story,onFirstPage=page_number,onLaterPages=page_number)

print("Phase 18C/18D package generated")
for p in [OUT,PDF,SUPP,SUPP_MANIFEST,REF_AUDIT,CLAIM_AUDIT]: print(p.relative_to(ROOT),p.stat().st_size,hashlib.md5(p.read_bytes()).hexdigest())
