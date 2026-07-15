# Post-GWAS mapping of kidney stone disease genetic risk to renal papillary contexts

**Article type:** Research article

**Authors:** Zhongxin Huang^1^, Xiaolu Duan^2^, Guohua Zeng^3^, Minglian Pan^2^, Cong Wang^2^, Tianyu Huang^2^, Mingyong Li^1^

**Affiliations:**

1. Department of Urology, The First Affiliated Hospital, Hengyang Medical School, University of South China, Hengyang, China.
2. Guangdong Provincial Key Laboratory of Urological Diseases, Department of Urology, The First Affiliated Hospital of Guangzhou Medical University, Guangzhou, Guangdong, People's Republic of China.
3. Department of Urology, Minimally Invasive Surgery Center, Guangdong Key Laboratory of Urology, Guangzhou Urology Research Institute, The First Affiliated Hospital of Guangzhou Medical University, Guangzhou, China.

**Corresponding authors:** Xiaolu Duan (94302304@qq.com) and Guohua Zeng (gzgyzgh@vip.sina.com)

## Abstract

**Background:** Kidney stone disease (KSD) has a substantial polygenic component, but the renal cell and tissue contexts of associated genes remain incompletely resolved. We integrated trans-ancestry KSD genome-wide association summary statistics with renal papillary single-nucleus, spatial and paired bulk transcriptomic resources. MAGMA defined an EUR-linkage-disequilibrium-reference-based genetic-priority layer, and GTEx Kidney_Cortex S-PrediXcan provided proxy annotation.

**Results:** In GSE231569, 43,878 nuclei from four papillary donors, including 540 Loop of Henle/thick ascending limb (Loop/TAL) nuclei, placed genetically prioritized modules in a donor-level Loop/TAL-associated context in expression- and detection-matched random benchmarks; the broader pattern attenuated but remained after known-driver removal. This cellular context was anchored by a genetic-priority layer in which quality control retained 4,915,033 of 5,960,489 GWAS rows, including 3,209 genome-wide significant variants, and identified 60 lead variants in 57 reconstructed loci; MAGMA tested 17,316 genes and identified 94 Bonferroni-significant, 369 false-discovery-rate (FDR)-significant and 187 suggestive genes. Across ten complete spatial sections, comprising four GSE206306 sections and six GSE231630 sections, label transfer provided broad-compartment orientation but not Loop/TAL or lesion enrichment. Kidney_Cortex S-PrediXcan provided proxy annotation across 5,989 tested genes, with 51 FDR-supported models comprising 42 one-SNP and nine multi-SNP models. The evidence-stratified 235-gene reporting set comprised R1-R6 groups of 68, 1, 2, 141, 16 and 7 genes, respectively. In GSE73680, 55 samples from 29 patients, including 26 paired patients, showed concordant positive module coefficients that did not remain FDR-significant and attenuated in sensitivity analyses using a KIM1/LCN2-like injury proxy and related tissue-state scores.

**Conclusions:** Genetically prioritized KSD genes map to a donor-level Loop/TAL-associated renal papillary context and are embedded in a broader KIM1/LCN2-like injury-associated bulk disease context. This post-GWAS framework provides cellular and tissue-context interpretation while leaving causal genes, causal cell types, plaque localization and papilla-specific regulation unresolved.

## Keywords

kidney stone disease; genome-wide association study; MAGMA; renal papilla; Loop of Henle; thick ascending limb; single-nucleus RNA sequencing; spatial transcriptomics; TWAS

## Background

Kidney stone disease is common, recurrent and biologically heterogeneous. Its pathogenesis reflects interactions among urinary chemistry, epithelial transport, papillary microenvironments and inherited susceptibility.[1-4] Genome-wide association studies (GWAS) have identified susceptibility loci across populations, including signals near genes involved in calcium handling and tubular transport.[5-8] However, association loci do not by themselves resolve the relevant genes, renal cell populations or disease-associated tissue states.[3,9] Connecting genetic priority to these biological contexts remains a central challenge in KSD genomics.

The renal papilla provides a biologically relevant setting for this interpretation. Randall plaque develops within papillary structures, and papillary tissue from stone formers displays inflammatory, oxidative-stress and injury-related features.[10-15] Single-nucleus and spatial transcriptomic studies have further resolved the epithelial, stromal, vascular and immune organization of the human kidney and papilla.[16-18] Within this anatomy, the Loop of Henle and thick ascending limb (Loop/TAL) coordinate salt transport, urinary concentration and paracellular calcium and magnesium handling through uromodulin-, claudin- and calcium-sensing pathways.[19-24] These functions make Loop/TAL biology relevant to KSD, while remaining distinct from lesion-specific localization or mechanism.

Gene-based association and cell-expression mapping offer complementary routes from polygenic signals to biological programs.[25-27] Yet each evidence type answers a different question. Gene-based association depends on the ancestry and linkage-disequilibrium (LD) reference; transcriptome-wide association studies (TWAS) depend on prediction-model architecture and tissue match; and single-cell or spatial projections depend on donor structure, expression background and anatomical resolution.[28-31] Combining these layers without distinguishing their inferential roles can make contextual agreement appear stronger than the underlying data support.

Here, we developed an evidence-stratified post-GWAS framework to ask whether MAGMA-prioritized KSD genes map to coherent renal papillary contexts. We distinguished EUR-LD-reference-based genetic priority from Kidney_Cortex TWAS proxy annotation, evaluated donor-level single-nucleus patterns, examined paired bulk module responses with tissue-state sensitivity analyses, and used spatial data for supplementary broad-compartment projection. This design connects genetic priority to cellular and tissue contexts while preserving the distinct inferential role of each evidence layer.

## Results

### GWAS and MAGMA analyses define an EUR-LD-reference-based genetic-priority layer

The GWAS input contained 5,960,489 rows. Quality control retained 4,915,033 rows, including 3,209 variants at genome-wide significance. Figure 1 provides a quality-control Manhattan visualization of the cleaned summary statistics, with genome-wide significant variants retained; locus-level and MAGMA-ranked gene information is provided in Supplementary Figure S1 and Supplementary Table 1. Distance pruning identified 60 lead SNPs, and merging overlapping windows produced 57 reconstructed loci. The source publication reported 59 loci,[5] but the publication-level lead-SNP or locus table required to reconcile the two counts was unavailable. We therefore report the 57-locus result as a reproducible reconstruction rather than a replacement for the published locus definition.

MAGMA v1.10 tested 17,316 genes using GRCh37 gene locations and the 1000 Genomes European reference panel. Ninety-four genes met Bonferroni significance, 369 met FDR significance and 187 met the prespecified suggestive threshold of P < 1 x 10^-4 (Supplementary Figure S1; Supplementary Table 2). The top-50, top-100, Bonferroni, FDR and suggestive modules were defined before downstream scoring and retained without post hoc membership changes. Because the source GWAS was trans-ancestry whereas the linkage-disequilibrium reference was European, these results constitute an EUR-LD-reference-based genetic-priority layer rather than ancestry-generalizable fine mapping or causal-gene evidence.

### Donor-level single-nucleus analyses map prioritized modules to a Loop/TAL-associated context

The GSE231569 papillary single-nucleus dataset contained 43,878 nuclei from four donors. Harmonized broad-compartment annotation identified 540 nuclei in the Loop/TAL compartment. Donor-by-compartment summaries placed the prespecified MAGMA modules consistently toward the Loop/TAL end of the expression distribution, with donors serving as the biological support units rather than individual nuclei (Figure 2; Supplementary Table 3).

The direction of the Loop/TAL-versus-other-compartment contrast remained stable in leave-one-donor-out analyses. Expression- and detection-matched random benchmarks supported the primary donor-level Loop/TAL-minus-other contrast across canonical modules, whereas rank statistics showed ceiling effects. Panel-level removal of prespecified Loop/TAL and renal-transport drivers attenuated but did not abolish the broader module-level pattern; retention after removal does not prove independence from Loop/TAL biology. Taken together, these analyses support a moderate donor-level Loop/TAL-associated context. They do not establish Loop/TAL as a unique or causal cell type for KSD genetic effects.

### Spatial data provide descriptive broad-compartment tissue-context projection

Ten complete spatial sections were retained for descriptive projection, comprising four GSE206306 sections and six GSE231630 sections. Single-nucleus-to-spatial label transfer was usable for broad-compartment orientation, although two sections showed lower confidence. Loop/TAL predictions were zero or sparse across sections, so the spatial data could not support Loop/TAL enrichment, Loop/TAL co-distribution or cross-section replication claims. No plaque, mineral, calcification, fibrosis or lesion region-of-interest annotation was available.

MAGMA module scores were overlaid descriptively on the tissue sections and interpreted only against broad transferred compartments (Supplementary Figure S2). These overlays provide supplementary tissue-context projection within the assayed papillary sections; the absence of lesion regions and usable Loop/TAL transfer precludes plaque localization, disease-control comparison or papilla-specific regulatory inference.

### An evidence-stratified candidate model separates MAGMA priority from Kidney_Cortex proxy annotation

S-PrediXcan with GTEx v8 Kidney_Cortex models tested 5,989 genes and identified 51 FDR-supported models. Forty-two were one-SNP models and nine were multi-SNP models, making support predominantly variant-limited. Because Kidney_Cortex does not match the renal papilla, TWAS was retained as proxy annotation and did not confer papilla-specific regulatory or causal status (Supplementary Figure S3; Supplementary Table 4).

The candidate universe contained 235 genes assigned deterministically to six mutually exclusive reporting groups: R1, 68; R2, 1; R3, 2; R4, 141; R5, 16; and R6, 7 (Figure 3; Supplementary Figure S3; Supplementary Table 5). R1-R6 are reporting groups, not causal tiers. They organize combinations of MAGMA priority, proxy TWAS support and contextual review; curated exemplars remain interpretive annotations and do not alter assignment. Bulk status records paired disease-context and tissue-state review without upgrading a reporting group. Complete resources for formal SMR/HEIDI or colocalization were unavailable, and no such result was used.

### Paired bulk analyses place prioritized modules in a KIM1/LCN2-like injury-associated disease context

The reconstructed GSE73680 dataset comprised 55 papillary samples from 29 patients, including 26 patients with paired control/adjacent and plaque/stone-associated tissue. Group labels were derived from filenames, a metadata limitation considered in interpretation. Module scores were calculated as the arithmetic mean of gene-wise z scores across detected members, without an additional log transformation. In patient fixed-effect models, all five prespecified MAGMA modules had positive group coefficients (0.169-0.174), but none remained significant after FDR correction (FDR 0.071-0.082) (Figure 4; Supplementary Table 6).

The paired module changes tracked broader tissue-state signatures. The KIM1/LCN2-like injury proxy showed the strongest nominal positive shift (beta = 0.375, P = 0.0057, FDR = 0.057); epithelial and mineralization/remodeling proxies were weaker (beta = 0.348, P = 0.0318, FDR = 0.106; and beta = 0.292, P = 0.0210, FDR = 0.105, respectively). Adjustment and paired-delta sensitivity analyses attenuated the module effects. Because these signatures are bulk expression scores and not estimated cell fractions, the results indicate embedding within a KIM1/LCN2-like injury-associated tissue-state context, not cell-composition-independent regulation, deconvolution or causal validation.

### Integrated evidence supports context mapping while preserving inferential boundaries

Across layers, MAGMA provided the genetic-priority anchor; donor-level snRNA data supplied the clearest renal cell-context signal; spatial data contributed supplementary broad-compartment anatomy; Kidney_Cortex TWAS supplied limited proxy annotation; and paired bulk data linked the modules to a broader disease-associated tissue state. The layers are complementary but not interchangeable, and agreement between them does not constitute causal triangulation.

The integrated result is therefore a post-GWAS renal papillary context map. It prioritizes a Loop/TAL-associated single-nucleus context and an attenuated KIM1/LCN2-like injury-associated bulk background while explicitly leaving causal genes, causal cell populations, plaque-specific localization and papilla-specific regulatory mechanisms unresolved.

## Discussion

This study establishes an evidence-stratified post-GWAS context map for KSD that connects gene-based genetic priority to renal papillary cellular and tissue states without conflating distinct evidence classes. Donor-level snRNA data provided the strongest contextual signal by placing prioritized modules in a Loop/TAL-associated setting; MAGMA anchored genetic priority, spatial data supplied supplementary anatomical orientation, Kidney_Cortex TWAS contributed proxy annotation and paired bulk data situated the modules within a disease-associated tissue background. This explicit separation of roles makes cross-layer agreement interpretable while preventing contextual concordance from being mistaken for causal evidence.

The Loop/TAL-associated pattern is biologically compatible with renal transport and ion-handling pathways implicated in stone susceptibility. Established genes including *UMOD*, *CASR*, *CLDN14* and *CLDN10* provide biological anchors,[19-24] although their familiarity cannot identify the gene mediating a GWAS locus. Expression- and detection-matched random benchmarks supported the donor-level Loop/TAL-minus-other contrast across canonical modules, and known-driver removal attenuated but did not abolish the broader pattern. Together, these findings provide moderate contextual support; the four-donor design, including one donor with only four Loop/TAL nuclei, limits stronger generalization and does not establish independence from Loop/TAL biology.

The bulk results provide a complementary tissue-state view. Concordant positive coefficients across five modules, although not FDR-significant, placed genetically prioritized programs within a broader papillary response. Their attenuation in KIM1/LCN2-like injury-proxy and epithelial-state sensitivity analyses supports tissue-state embedding rather than an isolated disease-specific module. Because these expression scores are not cell fractions, the data cannot distinguish within-cell transcriptional regulation from changes in tissue composition.

Spatial and TWAS evidence extended the context map while retaining narrower inferential roles. Ten spatial sections from GSE206306 and GSE231630 supported broad-compartment orientation, but sparse or absent Loop/TAL transfer and missing lesion annotations precluded plaque- or Loop/TAL-specific inference. Kidney_Cortex TWAS likewise remained proxy annotation because most supported models used one SNP and the tissue did not match the papilla. Complete SMR/HEIDI or colocalization inputs were unavailable, identifying the locus-resolved and papilla-relevant resources needed for stronger regulatory interpretation.

Future work should combine ancestry-matched genetic analyses, locus-level fine mapping and papilla-relevant eQTL resources with lesion-resolved spatial profiling. Donor-aware perturbation in relevant tubular systems could then test whether selected genes and modules exert direct effects or mark downstream tissue responses. These steps are required to move from contextual prioritization to mechanism.

## Limitations

The genetic layer was constrained by using a 1000 Genomes European LD reference for a trans-ancestry GWAS and by the inability to reconcile 57 reconstructed loci with the 59 published loci without the source lead-SNP or locus table. The donor-level cell-context result was supported by four donors and 540 Loop/TAL nuclei, with one donor contributing only four such nuclei, and matched-random support attenuated after known-driver removal. Spatial projection remained descriptive: two sections had lower transfer confidence, Loop/TAL prediction was unusable and lesion annotations were absent. These limitations narrow ancestry, cell-type and spatial generalization but do not alter the observed donor-level Loop/TAL-associated context.

Regulatory and disease-context interpretation was also bounded by the available resources. Kidney_Cortex TWAS was tissue-mismatched and dominated by one-SNP models, and complete inputs for formal SMR/HEIDI and colocalization were unavailable. The bulk dataset was reconstructed from microarray data using filename-derived group labels; tissue-state scores were expression proxies rather than cell fractions, and module associations did not survive FDR correction. Because all analyses were observational, perturbational experiments remain necessary to test whether prioritized genes, modules or cellular contexts have direct mechanistic roles.

## Methods

### Study design and evidence hierarchy

The study integrated the 2025 trans-ancestry KSD GWAS,[5] GSE231569 papillary snRNA data,[16] ten spatial sections comprising four GSE206306 and six GSE231630 sections,[16] GSE73680 paired papillary bulk expression,[12] and GTEx v8 Kidney_Cortex prediction models.[29] Analyses were assigned prospectively by inferential role: GWAS/MAGMA for genetic priority, snRNA for donor-level cell context, spatial data for supplementary broad-compartment projection, Kidney_Cortex TWAS for proxy annotation and bulk data for paired disease-context sensitivity.

### GWAS quality control and locus reconstruction

GRCh37-compatible summary statistics were restricted to valid biallelic A/C/G/T SNP records with valid coordinates and P values; strand-ambiguous variants and duplicate or malformed records were removed, with low-frequency filtering applied where allele frequency was available. Genome-wide significant variants were greedily pruned within +/-1 Mb and overlapping lead windows were merged. The procedure retained 4,915,033 of 5,960,489 rows, identified 3,209 genome-wide significant variants and 60 lead SNPs, and reconstructed 57 loci. The reconstruction was not used to explain the difference from the 59 loci reported by the source study.

### MAGMA analysis and canonical module freeze

MAGMA v1.10 used NCBI build 37 gene locations and the 1000 Genomes European LD reference.[25] Of 17,316 tested genes, significance was defined using Bonferroni correction, Benjamini-Hochberg FDR,[32] and a suggestive threshold of P < 1 x 10^-4. Top-50, top-100, Bonferroni, FDR and suggestive lists were defined before downstream scoring. Executable, version, gene-location, LD-prefix, SNP/P input, gene output and log provenance were recorded for reproducibility. MAGMA results were interpreted as EUR-LD-reference-based genetic priority.

### Single-nucleus processing, scoring and sensitivity analyses

Four papillary donors from GSE231569 were retained, yielding 43,878 nuclei after quality control and broad-compartment harmonization. Seurat was used for single-nucleus object handling, normalization and broad-compartment workflows.[33] The Loop/TAL compartment contained 540 nuclei. For each prespecified module, nucleus-level scores were calculated as arithmetic means of log-normalized expression among detected genes and summarized by donor and compartment. Donors were biological support units. Leave-one-donor-out summaries, 1,000 expression- and detection-matched random gene sets per module, low-cell-count sensitivity, and panel-level removal of prespecified Loop/TAL/transport drivers evaluated donor dependence, matched expectation and driver dependence.

### Spatial broad-compartment projection and descriptive module overlays

Ten complete sections with expression, coordinate and image assets were processed as Seurat spatial objects, comprising four GSE206306 and six GSE231630 sections.[33] The GSE231569 reference supported label transfer only at the broad-compartment level. Transfer confidence and prediction sparsity were assessed by section. Because transferred Loop/TAL prediction scores were zero or nonzero-sparse and no lesion-region annotations were available, module scores were displayed as descriptive overlays without formal Loop/TAL, plaque, disease-control or cross-section enrichment tests. Sections, not spots, defined the highest available biological unit.

### Kidney_Cortex proxy TWAS and candidate reporting groups

S-PrediXcan used GTEx v8 Kidney_Cortex MASHR models.[29-30] Benjamini-Hochberg correction was applied across 5,989 tested genes, and 51 FDR-supported models were classified by the number of contributing variants as one-SNP (42) or multi-SNP (9). The tissue was treated as a cortex proxy, and TWAS did not confer causal or papilla-specific regulatory status. A deterministic union and deduplication procedure produced 235 candidate genes and mutually exclusive R1-R6 reporting groups (68, 1, 2, 141, 16 and 7 genes). The groups are reporting categories, not causal tiers. Feasibility of SMR/HEIDI[34] and colocalization[35-36] was evaluated, but incomplete inputs precluded formal analyses.

### Paired bulk reconstruction, module scoring and tissue-state sensitivity

GSE73680 feature-level expression and filename-derived group labels yielded 55 samples from 29 patients, including 26 complete pairs.[12,37] Filename-based metadata provenance was retained as an interpretive limitation. Module scores were arithmetic means of gene-wise z scores across the 55 samples for detected module genes. Paired deltas were defined as plaque/stone-associated minus control/adjacent, and patient fixed-effect models used `module_score ~ group_binary + factor(patient_id)`. Benjamini-Hochberg correction was applied across the five canonical modules.

KIM1/LCN2-like injury, epithelial, extracellular-matrix/fibrosis, immune/inflammatory and mineralization/remodeling marker scores were used as tissue-state expression proxies. Paired-delta and patient-aware covariate models tested whether module effects persisted after adjustment. These scores were not interpreted as cell fractions, and the analyses were not described as deconvolution or mediation.

### Statistical reporting and reproducibility

Tests were two-sided where applicable, with correction families defined within each analysis layer. Donors, patients and sections, rather than nuclei, samples or spots, were treated as biological support units for cross-unit interpretation. Machine-readable source data, manifests, logs and scripts were retained for reported tables and figures. No SMR, colocalization, plaque-localization or spatial Loop/TAL result was included among the reported evidence.

## Conclusions

Genetically prioritized KSD genes map to a donor-level Loop/TAL-associated renal papillary context and an attenuated KIM1/LCN2-like injury-associated bulk disease context. This framework provides post-GWAS functional interpretation while preserving the distinction between contextual support and causal gene, cell-type, lesion-localization or papilla-specific regulatory evidence.

## List of abbreviations

BH: Benjamini-Hochberg; ECM: extracellular matrix; eQTL: expression quantitative trait locus; FDR: false discovery rate; GEO: Gene Expression Omnibus; GWAS: genome-wide association study; KSD: kidney stone disease; LD: linkage disequilibrium; MAGMA: Multi-marker Analysis of GenoMic Annotation; QC: quality control; ROI: region of interest; SMR: summary-data Mendelian randomization; SNP: single-nucleotide polymorphism; snRNA-seq: single-nucleus RNA sequencing; TAL: thick ascending limb; TWAS: transcriptome-wide association study; UMAP: uniform manifold approximation and projection.

## Declarations

### Ethics approval and consent to participate

This study re-analysed publicly available, de-identified GWAS summary statistics and transcriptomic datasets. No new human participants, human tissue samples or animal experiments were involved. The original studies obtained ethics approval and informed consent as described in their respective publications. Institutional ethics review was not required for this secondary analysis of publicly available de-identified data. [TO VERIFY WITH INSTITUTION]

### Consent for publication

Not applicable.

### Availability of data and materials

The Cao et al. KSD GWAS summary statistics are available from Zenodo (https://doi.org/10.5281/zenodo.14790324), and the associated GWAS Catalog study record is GCST90652506. Public transcriptomic data are available from the Gene Expression Omnibus under accessions GSE231569, GSE73680, GSE206306 and GSE231630. GTEx v8/PredictDB Kidney_Cortex MASHR models are available from their original distribution resource. Analysis code, derived result tables, figure source data, Supplementary Tables, manifests and selected logs are available at the GitHub repository (https://github.com/huangzhongxin411/ksd-papillary-context-mapping). The current manuscript-synchronized release will be v1.0.1 after final approval; GitHub Release v1.0.0 is retained as a pre-Zenodo checkpoint. A manuscript-synchronized repository archive and DOI will be added after final release and Zenodo minting.

### Competing interests

The authors declare that they have no competing interests.

### Funding

This work was supported by the National Natural Science Foundation of China (No. 82270798). The funder had no role in study design, data analysis, interpretation of data, manuscript preparation or the decision to submit the work for publication.

### Authors' contributions

[TO FILL: author initials and CRediT contribution statement]

### Acknowledgements

The authors thank the investigators who generated and shared the public GWAS, single-nucleus, bulk and spatial transcriptomic datasets re-analysed in this study.

## References

1. Khan SR, Pearle MS, Robertson WG, Gambaro G, Canales BK, Doizi S, et al. Kidney stones. Nat Rev Dis Primers. 2016;2:16008. https://doi.org/10.1038/nrdp.2016.8

2. Moe OW. Kidney stones: pathophysiology and medical management. Lancet. 2006;367(9507):333-44. https://doi.org/10.1016/S0140-6736(06)68071-9

3. Singh P, Harris PC, Sas DJ, Lieske JC. The genetics of kidney stone disease and nephrocalcinosis. Nat Rev Nephrol. 2022;18(4):224-240. https://doi.org/10.1038/s41581-021-00513-4

4. Howles SA, Thakker RV. Genetics of kidney stone disease. Nat Rev Urol. 2020;17(7):407-421. https://doi.org/10.1038/s41585-020-0332-x

5. Cao X, Jiang M, Guan Y, Li S, Duan C, Gong Y, et al. Trans-ancestry GWAS identifies 59 loci and improves risk prediction and fine-mapping for kidney stone disease. Nat Commun. 2025;16(1):3473. https://doi.org/10.1038/s41467-025-58782-7

6. Hao X, Shao Z, Zhang N, Jiang M, Cao X, Li S, et al. Integrative genome-wide analyses identify novel loci associated with kidney stones and provide insights into its genetic architecture. Nat Commun. 2023;14(1):7498. https://doi.org/10.1038/s41467-023-43400-1

7. Chen WC, Chen YC, Chen YH, Liu TY, Tsai CH, Tsai FJ. Identification of novel genetic susceptibility loci for calcium-containing kidney stone disease by genome-wide association study and polygenic risk score in a Taiwanese population. Urolithiasis. 2024;52(1):94. https://doi.org/10.1007/s00240-024-01577-0

8. Thorleifsson G, Holm H, Edvardsson V, Walters GB, Styrkarsdottir U, Gudbjartsson DF, et al. Sequence variants in the CLDN14 gene associate with kidney stones and bone mineral density. Nat Genet. 2009;41(8):926-30. https://doi.org/10.1038/ng.404

9. Palsson R, Indridason OS, Edvardsson VO, Oddsson A. Genetics of common complex kidney stone disease: insights from genome-wide association studies. Urolithiasis. 2019;47(1):11-21. https://doi.org/10.1007/s00240-018-1094-2

10. Evan AP, Lingeman JE, Coe FL, Parks JH, Bledsoe SB, Shao Y, et al. Randall's plaque of patients with nephrolithiasis begins in basement membranes of thin loops of Henle. J Clin Invest. 2003;111(5):607-16. https://doi.org/10.1172/JCI17038

11. Evan AP, Coe FL, Lingeman J, Bledsoe S, Worcester EM. Randall's plaque in stone formers originates in ascending thin limbs. Am J Physiol Renal Physiol. 2018;315(5):F1236-F1242. https://doi.org/10.1152/ajprenal.00035.2018

12. Taguchi K, Hamamoto S, Okada A, Unno R, Kamisawa H, Naiki T, et al. Genome-Wide Gene Expression Profiling of Randall's Plaques in Calcium Oxalate Stone Formers. J Am Soc Nephrol. 2017;28(1):333-347. https://doi.org/10.1681/ASN.2015111271

13. Taguchi K, Okada A, Hamamoto S, Unno R, Moritoki Y, Ando R, et al. M1/M2-macrophage phenotypes regulate renal calcium oxalate crystal development. Sci Rep. 2016;6:35167. https://doi.org/10.1038/srep35167

14. Khan SR, Canales BK, Dominguez-Gutierrez PR. Randall's plaque and calcium oxalate stone formation: role for immunity and inflammation. Nat Rev Nephrol. 2021;17(6):417-433. https://doi.org/10.1038/s41581-020-00392-1

15. Marien TP, Miller NL. Characteristics of renal papillae in kidney stone formers. Minerva Urol Nefrol. 2016;68(6):496-515. https://pubmed.ncbi.nlm.nih.gov/27441596/

16. Canela VH, Bowen WS, Ferreira RM, Syed F, Lingeman JE, Sabo AR, et al. A spatially anchored transcriptomic atlas of the human kidney papilla identifies significant immune injury in patients with stone disease. Nat Commun. 2023;14(1):4140. https://doi.org/10.1038/s41467-023-38975-8

17. Lake BB, Menon R, Winfree S, Hu Q, Melo Ferreira R, Kalhor K, et al. An atlas of healthy and injured cell states and niches in the human kidney. Nature. 2023;619(7970):585-594. https://doi.org/10.1038/s41586-023-05769-3

18. Lake BB, Chen S, Hoshi M, Plongthongkum N, Salamon D, Knoten A, et al. A single-nucleus RNA-sequencing pipeline to decipher the molecular anatomy and pathophysiology of human kidneys. Nat Commun. 2019;10(1):2832. https://doi.org/10.1038/s41467-019-10861-2

19. Bleich M, Wulfmeyer VC, Himmerkus N, Milatz S. Heterogeneity of tight junctions in the thick ascending limb. Ann N Y Acad Sci. 2017;1405(1):5-15. https://doi.org/10.1111/nyas.13400

20. Yu AS. Claudins and the kidney. J Am Soc Nephrol. 2015;26(1):11-9. https://doi.org/10.1681/ASN.2014030284

21. Breiderhoff T, Himmerkus N, Stuiver M, Mutig K, Will C, Meij IC, et al. Deletion of claudin-10 (Cldn10) in the thick ascending limb impairs paracellular sodium permeability and leads to hypermagnesemia and nephrocalcinosis. Proc Natl Acad Sci U S A. 2012;109(35):14241-6. https://doi.org/10.1073/pnas.1203834109

22. Wolf MTF, Zhang J, Nie M. Uromodulin in mineral metabolism. Curr Opin Nephrol Hypertens. 2019;28(5):481-489. https://doi.org/10.1097/MNH.0000000000000522

23. Guha M, Bankura B, Ghosh S, Pattanayak AK, Ghosh S, Pal DK, et al. Polymorphisms in CaSR and CLDN14 Genes Associated with Increased Risk of Kidney Stone Disease in Patients from the Eastern Part of India. PLoS One. 2015;10(6):e0130790. https://doi.org/10.1371/journal.pone.0130790

24. Alexander RT, Dimke H. Effect of diuretics on renal tubular transport of calcium and magnesium. Am J Physiol Renal Physiol. 2017;312(6):F998-F1015. https://doi.org/10.1152/ajprenal.00032.2017

25. de Leeuw CA, Mooij JM, Heskes T, Posthuma D. MAGMA: generalized gene-set analysis of GWAS data. PLoS Comput Biol. 2015;11(4):e1004219. https://doi.org/10.1371/journal.pcbi.1004219

26. Finucane HK, Reshef YA, Anttila V, Slowikowski K, Gusev A, Byrnes A, et al. Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. Nat Genet. 2018;50(4):621-629. https://doi.org/10.1038/s41588-018-0081-4

27. Skene NG, Bryois J, Bakken TE, Breen G, Crowley JJ, Gaspar HA, et al. Genetic identification of brain cell types underlying schizophrenia. Nat Genet. 2018;50(6):825-833. https://doi.org/10.1038/s41588-018-0129-5

28. Gusev A, Ko A, Shi H, Bhatia G, Chung W, Penninx BW, et al. Integrative approaches for large-scale transcriptome-wide association studies. Nat Genet. 2016;48(3):245-52. https://doi.org/10.1038/ng.3506

29. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science. 2020;369(6509):1318-1330. https://doi.org/10.1126/science.aaz1776

30. Barbeira AN, Dickinson SP, Bonazzola R, Zheng J, Wheeler HE, Torres JM, et al. Exploring the phenotypic consequences of tissue specific gene expression variation inferred from GWAS summary statistics. Nat Commun. 2018;9(1):1825. https://doi.org/10.1038/s41467-018-03621-1

31. Wainberg M, Sinnott-Armstrong N, Mancuso N, Barbeira AN, Knowles DA, Golan D, et al. Opportunities and challenges for transcriptome-wide association studies. Nat Genet. 2019;51(4):592-599. https://doi.org/10.1038/s41588-019-0385-z

32. Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

33. Hao Y, Hao S, Andersen-Nissen E, Mauck WM, Zheng S, Butler A, et al. Integrated analysis of multimodal single-cell data. Cell. 2021;184(13):3573-3587.e29. https://doi.org/10.1016/j.cell.2021.04.048

34. Zhu Z, Zhang F, Hu H, Bakshi A, Robinson MR, Powell JE, et al. Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. Nat Genet. 2016;48(5):481-7. https://doi.org/10.1038/ng.3538

35. Giambartolomei C, Vukcevic D, Schadt EE, Franke L, Hingorani AD, Wallace C, et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet. 2014;10(5):e1004383. https://doi.org/10.1371/journal.pgen.1004383

36. Wallace C. Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses. PLoS Genet. 2020;16(4):e1008720. https://doi.org/10.1371/journal.pgen.1008720

37. Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. Nucleic Acids Res. 2002;30(1):207-10. https://doi.org/10.1093/nar/30.1.207

## Figures and legends

### Figure 1. GWAS quality-control Manhattan plot

GWAS quality-control Manhattan plot generated from downsampled cleaned trans-ancestry KSD summary statistics, with genome-wide significant variants retained. The horizontal reference line marks genome-wide significance. The panel visualizes GWAS input behavior and is not a fine-mapping or causal-gene inference display. Locus-level and MAGMA-ranked gene information is provided in Supplementary Figure S1 and Supplementary Table 1.

### Figure 2. Donor-level snRNA context of genetically prioritized modules

Donor-aware analysis of 43,878 nuclei from four GSE231569 papillary donors, including 540 Loop/TAL nuclei. (A) snRNA resource summary. (B) Broad-compartment UMAP. (C) Loop/TAL marker-expression context. (D) Donor-by-compartment scores for prespecified MAGMA modules. (E) Expression- and detection-matched random benchmark for the Loop/TAL-minus-other contrast. (F) Sensitivity after removal of prespecified Loop/TAL and renal-transport drivers. Donor x compartment is the biological unit. The benchmark and driver-removal analyses support a moderate Loop/TAL-associated context and do not establish a causal cell type.

### Figure 3. Evidence-stratified candidate reporting model

Evidence-stratified organization of the candidate reporting model. (A) Parallel roles of genetic priority, cell context, proxy annotation, tissue context and disease-state review. (B) Representative rows from the candidate evidence matrix. (C) Kidney_Cortex proxy-model burden, including one-SNP and multi-SNP classes. (D) Counts across R1-R6. (E) Interpretation boundary. R1-R6 are reporting groups, not causal tiers. Kidney_Cortex TWAS is proxy annotation, and contextual evidence does not upgrade causal status. The complete 235-gene matrix is shown in Supplementary Figure S3 and Supplementary Table 5.

### Figure 4. Paired bulk disease-context and tissue-state analysis

Paired GSE73680 analysis of 55 samples from 29 patients, including 26 complete pairs. (A) Compact paired-study design. (B) Selected paired module scores. (C) Base disease-context coefficients for prespecified MAGMA modules. (D) Paired tissue-state proxy shifts. (E) Coupling between module and tissue-state paired differences. (F) Attenuation of disease-context coefficients after tissue-state adjustment. Positive module coefficients did not remain FDR-significant and attenuated in sensitivity analyses, placing prioritized modules in a broader KIM1/LCN2-like injury-associated tissue state. Expression scores are not cell fractions, and the analysis does not provide causal validation.

### Supplementary Figure S1. GWAS and MAGMA diagnostics

GWAS and MAGMA diagnostic panels complement the main-text quality-control Manhattan plot. The panels show GWAS quantile-quantile and P-value behavior together with MAGMA significance and rank summaries; the top 20 MAGMA-ranked genes are displayed in panel S1C. These genes represent gene-based priority and should not be interpreted as causal genes.

### Supplementary Figure S2. Spatial broad-compartment tissue-context projection

Five-page supplementary spatial projection across ten complete sections, comprising four GSE206306 and six GSE231630 sections. S2A presents spatial-resource and label-transfer quality control. S2B1-S2B3 show Bonferroni, Top 100 and suggestive module overlays, and S2B4 summarizes section-level module scores by predicted broad compartment. The transferred labels provide broad tissue context only: no lesion ROI was available, Loop/TAL prediction was unusable, and no Loop/TAL enrichment, disease-control test or plaque localization was performed.

### Supplementary Figure S3. Kidney_Cortex TWAS proxy and candidate evidence matrix

Thirteen-page supplement separating TWAS proxy quality from the complete candidate matrix. S3A summarizes 5,989 tested Kidney_Cortex models and 51 FDR-supported results, including 42 one-SNP and nine multi-SNP models. S3B paginates the complete 235-gene matrix by R group, with every cell rendered directly from Supplementary Table 5. R1-R6 are reporting groups, not causal tiers, and Kidney_Cortex TWAS remains proxy annotation rather than papilla-specific regulatory evidence.
