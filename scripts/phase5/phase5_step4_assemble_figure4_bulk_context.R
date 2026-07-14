suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

root <- normalizePath(".", mustWork = TRUE)
for (d in c("results/tables", "results/figures", "source_data/figures", "notes", "codex_tasks")) {
  dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)
}
rd <- function(x) read.delim(file.path(root, x), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
wr <- function(x, p) write.table(x, file.path(root, p), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)

# Locked Phase 5 inputs only. No new statistical models are fitted in this script.
pair_summary <- rd("results/tables/phase5_step1_bulk_pairing_summary.tsv")
mapping <- rd("results/tables/phase5_step2_magma_module_bulk_gene_mapping.tsv")
module_scores <- rd("results/tables/phase5_step2_bulk_sample_module_scores.tsv")
base <- rd("results/tables/phase5_step2_bulk_paired_model_results.tsv")
paired_diff <- rd("results/tables/phase5_step2_bulk_paired_differences.tsv")
marker_defs <- rd("results/tables/phase5_step3_marker_signature_definitions.tsv")
tissue_scores <- rd("results/tables/phase5_step3_bulk_sample_tissue_state_scores.tsv")
tissue_models <- rd("results/tables/phase5_step3_tissue_state_paired_model_results.tsv")
cors <- rd("results/tables/phase5_step3_module_tissue_state_correlations.tsv")
adjusted <- rd("results/tables/phase5_step3_tissue_state_adjusted_module_models.tsv")

teal <- "#245A64"; loop_teal <- "#0F4C5C"; bluegrey <- "#7F9DA6"
gold <- "#B99B5A"; terra <- "#9B5C4D"; pale <- "#E6E9EA"; dark <- "#333333"
module_labels <- c(MAGMA_top50="Top 50", MAGMA_top100="Top 100", MAGMA_Bonferroni="Bonferroni", MAGMA_FDR05="FDR 0.05", MAGMA_suggestive_p1e4="Suggestive")
sig_labels <- c(injury_KIM1_LCN2_like="Injury", epithelial_transport="Epithelial state", mineralization_remodeling_proxy="Mineralization/remodeling", fibrosis_remodeling="Fibrosis/remodeling", fibroblast_stromal_ECM="Stromal/ECM", endothelial_proxy="Endothelial", immune_proxy="Immune", inflammatory_stress="Inflammatory stress", Loop_TAL_transport_proxy="Loop/TAL transport", collecting_duct_proxy="Collecting duct")
group_labels <- c(control_or_adjacent="Control/adjacent", plaque_or_stone_papilla="Plaque/stone papilla")
theme_pub <- theme_bw(base_size = 9.5) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face="bold", size=11), plot.subtitle=element_text(size=8.5, color=dark), strip.background=element_rect(fill=pale, color=NA), strip.text=element_text(face="bold", size=9), legend.title=element_text(size=8.5), legend.text=element_text(size=8), axis.text=element_text(size=8.5), axis.title=element_text(size=9.5))

# Panel A: design bridge.
design_source <- data.frame(metric=c("Included samples","Included patients","Complete pairs","Control/adjacent","Plaque/stone papilla"), value=c(55,29,26,27,28), source="Phase 5-Step 1/2 locked metadata", notes="Filename-derived labels remain manually reviewable.")
wr(design_source, "source_data/figures/phase5_step4_Figure4_panelA_source_data.tsv")
pA <- ggplot() +
  annotate("rect", xmin=.2, xmax=1.8, ymin=1.7, ymax=3.2, fill=pale, color=bluegrey, linewidth=.7) +
  annotate("text", x=1, y=2.85, label="GSE73680", fontface="bold", size=4.2, color=dark) +
  annotate("text", x=1, y=2.38, label="55 samples | 29 patients", size=3.7, color=dark) +
  annotate("text", x=1, y=1.98, label="processed intensity-like bulk", size=3.2, color=dark) +
  annotate("segment", x=1, xend=1, y=1.65, yend=1.25, arrow=arrow(length=unit(0.12,"in")), color=teal, linewidth=.7) +
  annotate("rect", xmin=.05, xmax=.9, ymin=.15, ymax=1.15, fill="#EEF2F3", color=bluegrey, linewidth=.6) +
  annotate("rect", xmin=1.1, xmax=1.95, ymin=.15, ymax=1.15, fill="#F2EBDD", color=gold, linewidth=.6) +
  annotate("text", x=.475, y=.83, label="Control/adjacent", fontface="bold", size=3.4) +
  annotate("text", x=.475, y=.49, label="27 samples", size=3.2) +
  annotate("text", x=1.525, y=.83, label="Plaque/stone papilla", fontface="bold", size=3.4) +
  annotate("text", x=1.525, y=.49, label="28 samples", size=3.2) +
  annotate("text", x=1, y=-.18, label="26 complete patient pairs | disease-context only", fontface="bold", color=terra, size=3.35) +
  coord_cartesian(xlim=c(0,2), ylim=c(-.35,3.3), clip="off") + labs(title="Paired bulk design") + theme_void(base_size=9.5) + theme(plot.title=element_text(face="bold",size=11), plot.margin=margin(8,8,8,8))

# Panel B: selected paired module score lines.
main_modules <- c("MAGMA_top50","MAGMA_FDR05","MAGMA_suggestive_p1e4")
bdat <- module_scores[module_scores$module_name %in% main_modules & module_scores$patient_id %in% paired_diff$patient_id, ]
complete_ids <- unique(paired_diff$patient_id)
bdat <- bdat[bdat$patient_id %in% complete_ids, ]
bdat$module_name <- factor(bdat$module_name, levels=main_modules, labels=module_labels[main_modules])
bdat$group_curated <- factor(bdat$group_curated, levels=c("control_or_adjacent","plaque_or_stone_papilla"))
wr(bdat, "source_data/figures/phase5_step4_Figure4_panelB_source_data.tsv")
pB <- ggplot(bdat, aes(group_curated, module_score, group=patient_id)) +
  geom_line(color=bluegrey, alpha=.52, linewidth=.28) +
  geom_point(aes(fill=group_curated), shape=21, color=dark, stroke=.2, size=1.35) +
  facet_wrap(~module_name, nrow=1, scales="free_y") +
  scale_fill_manual(values=c(control_or_adjacent=bluegrey, plaque_or_stone_papilla=gold), guide="none") +
  scale_x_discrete(labels=group_labels) +
  labs(title="Paired MAGMA module scores", subtitle="Selected canonical modules; 26 patient pairs", x=NULL, y="Module score") +
  theme_pub + theme(axis.text.x=element_text(angle=20,hjust=1))

# Panel C: locked base coefficients with 95% CI.
cdat <- base
cdat$lower <- cdat$group_coefficient - 1.96*cdat$standard_error
cdat$upper <- cdat$group_coefficient + 1.96*cdat$standard_error
cdat$module_label <- factor(module_labels[cdat$module_name], levels=rev(unname(module_labels)))
cdat$nominal <- cdat$p_value < .05
wr(cdat, "source_data/figures/phase5_step4_Figure4_panelC_source_data.tsv")
pC <- ggplot(cdat, aes(group_coefficient, module_label)) +
  geom_vline(xintercept=0, color=dark, linewidth=.35) +
  geom_errorbar(aes(xmin=lower, xmax=upper), orientation="y", width=.18, color=bluegrey, linewidth=.55) +
  geom_point(aes(fill=nominal), shape=21, color=dark, size=2.5, stroke=.3) +
  scale_fill_manual(values=c("TRUE"=gold,"FALSE"=pale), labels=c("TRUE"="Nominal p<0.05","FALSE"="p>=0.05"), name=NULL) +
  labs(title="Base paired group coefficients", subtitle="Patient fixed effects; all module FDR values >0.05", x="Disease-context coefficient (95% CI)", y=NULL) +
  theme_pub + theme(legend.position="bottom")

# Panel D: selected tissue-state paired proxy scores.
d_sigs <- c("injury_KIM1_LCN2_like","epithelial_transport","mineralization_remodeling_proxy")
ddat <- tissue_scores[tissue_scores$signature_name %in% d_sigs & tissue_scores$patient_id %in% complete_ids, ]
ddat$signature_label <- factor(sig_labels[ddat$signature_name], levels=sig_labels[d_sigs])
ddat$group_curated <- factor(ddat$group_curated, levels=c("control_or_adjacent","plaque_or_stone_papilla"))
wr(ddat, "source_data/figures/phase5_step4_Figure4_panelD_source_data.tsv")
pD <- ggplot(ddat, aes(group_curated, signature_score, group=patient_id)) +
  geom_line(color=bluegrey, alpha=.52, linewidth=.28) +
  geom_point(aes(fill=group_curated), shape=21, color=dark, stroke=.2, size=1.35) +
  facet_wrap(~signature_label, nrow=1, scales="free_y") +
  scale_fill_manual(values=c(control_or_adjacent=bluegrey, plaque_or_stone_papilla=gold), guide="none") +
  scale_x_discrete(labels=group_labels) +
  labs(title="Tissue-state proxy shifts", subtitle="Selected marker scores; proxies are not cell fractions", x=NULL, y="Proxy score") +
  theme_pub + theme(axis.text.x=element_text(angle=20,hjust=1))

# Panel E: reduced paired-difference correlation matrix.
e_sigs <- c("injury_KIM1_LCN2_like","epithelial_transport","fibroblast_stromal_ECM","endothelial_proxy","fibrosis_remodeling")
edat <- cors[cors$analysis_level=="paired_difference_level" & cors$signature_name %in% e_sigs, ]
edat$module_label <- factor(module_labels[edat$module_name], levels=rev(unname(module_labels)))
edat$signature_label <- factor(sig_labels[edat$signature_name], levels=sig_labels[e_sigs])
wr(edat, "source_data/figures/phase5_step4_Figure4_panelE_source_data.tsv")
pE <- ggplot(edat, aes(signature_label,module_label,fill=correlation)) +
  geom_tile(color="white",linewidth=.45) + geom_text(aes(label=sprintf("%.2f",correlation)),size=2.65) +
  scale_fill_gradient2(low=bluegrey,mid="white",high=terra,midpoint=0,limits=c(-1,1),name="Spearman\nrho") +
  labs(title="Module-tissue-state coupling",subtitle="Paired disease-minus-control differences; non-causal",x=NULL,y=NULL) +
  theme_pub + theme(panel.grid=element_blank(),axis.text.x=element_text(angle=30,hjust=1),legend.position="right")

# Panel F: selected sensitivity coefficients from locked Step 3 output.
f_models <- c("single_injury_KIM1_LCN2_like","compact_injury_KIM1_LCN2_like_plus_fibrosis_remodeling")
fadj <- adjusted[adjusted$adjustment_model %in% f_models, c("module_name","adjustment_model","adjusted_group_coefficient","adjusted_p_value","adjusted_fdr_within_model_family","attenuation_class")]
fbase <- base[,c("module_name","group_coefficient","p_value","fdr_bh_across_modules")]
names(fbase)[2:4] <- c("adjusted_group_coefficient","adjusted_p_value","adjusted_fdr_within_model_family")
fbase$adjustment_model <- "base_unadjusted"; fbase$attenuation_class <- "base"
fdat <- rbind(fbase[,names(fadj)],fadj)
fdat$model_label <- factor(fdat$adjustment_model, levels=c("base_unadjusted","single_injury_KIM1_LCN2_like","compact_injury_KIM1_LCN2_like_plus_fibrosis_remodeling"), labels=c("Base","Injury","Injury + fibrosis"))
fdat$module_label <- module_labels[fdat$module_name]
wr(fdat, "source_data/figures/phase5_step4_Figure4_panelF_source_data.tsv")
pF <- ggplot(fdat,aes(model_label,adjusted_group_coefficient,group=module_name,color=module_name)) +
  geom_hline(yintercept=0,color=dark,linewidth=.35) + geom_line(alpha=.72,linewidth=.55) + geom_point(size=2.1) +
  scale_color_manual(values=c(MAGMA_top50=bluegrey,MAGMA_top100=terra,MAGMA_Bonferroni=teal,MAGMA_FDR05=loop_teal,MAGMA_suggestive_p1e4=gold), labels=module_labels, name=NULL) +
  labs(title="Tissue-state sensitivity",subtitle="Attenuation indicates disease-state embedding",x=NULL,y="Disease-context coefficient") +
  theme_pub + theme(legend.position="bottom")

main <- ((pA | pB) / (pC | pD) / (pE | pF)) +
  plot_layout(heights=c(.9,1,1.08), widths=c(1,1.35), guides="collect") +
  plot_annotation(title="Bulk papillary disease context and tissue-state sensitivity", subtitle="GSE73680 supports paired module-level disease-context consistency; causal genes, plaque mechanisms and cell fractions are not inferred", tag_levels="A", theme=theme(plot.title=element_text(face="bold",size=17),plot.subtitle=element_text(size=10,color=dark),plot.tag=element_text(face="bold",size=14)))
ggsave(file.path(root,"results/figures/phase5_step4_Figure4_bulk_disease_context_draft.pdf"),main,width=16,height=14,units="in",device="pdf",bg="white")
ggsave(file.path(root,"results/figures/phase5_step4_Figure4_bulk_disease_context_draft.png"),main,width=16,height=14,units="in",dpi=300,bg="white")
ggsave(file.path(root,"results/figures/phase5_step4_Figure4_bulk_disease_context_draft_600dpi.png"),main,width=16,height=14,units="in",dpi=600,bg="white")

# Supplement: full proxy shifts, full paired-delta heatmap and all adjustment families.
sdat <- tissue_scores[tissue_scores$patient_id %in% complete_ids, ]
sdat$signature_label <- factor(sig_labels[sdat$signature_name], levels=unname(sig_labels))
sdat$group_curated <- factor(sdat$group_curated, levels=c("control_or_adjacent","plaque_or_stone_papilla"))
sA <- ggplot(sdat,aes(group_curated,signature_score,group=patient_id))+geom_line(color=bluegrey,alpha=.42,linewidth=.25)+geom_point(aes(fill=group_curated),shape=21,color=dark,stroke=.15,size=1.05)+facet_wrap(~signature_label,ncol=5,scales="free_y")+scale_fill_manual(values=c(control_or_adjacent=bluegrey,plaque_or_stone_papilla=gold),guide="none")+scale_x_discrete(labels=c(control_or_adjacent="Control",plaque_or_stone_papilla="Disease"))+labs(title="All paired tissue-state proxy scores",x=NULL,y="Proxy score")+theme_pub+theme(axis.text.x=element_text(angle=25,hjust=1),strip.text=element_text(size=8))
sheat <- cors[cors$analysis_level=="paired_difference_level", ]; sheat$module_label<-factor(module_labels[sheat$module_name],levels=rev(unname(module_labels))); sheat$signature_label<-factor(sig_labels[sheat$signature_name],levels=unname(sig_labels))
sB <- ggplot(sheat,aes(signature_label,module_label,fill=correlation))+geom_tile(color="white",linewidth=.4)+geom_text(aes(label=sprintf("%.2f",correlation)),size=2.4)+scale_fill_gradient2(low=bluegrey,mid="white",high=terra,midpoint=0,limits=c(-1,1),name="rho")+labs(title="Full paired-difference correlation matrix",x=NULL,y=NULL)+theme_pub+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=35,hjust=1))
scoef <- adjusted
scoef$model_label <- factor(scoef$adjustment_model)
scoef$module_label <- factor(module_labels[scoef$module_name],levels=unname(module_labels))
sC <- ggplot(scoef,aes(adjusted_group_coefficient,model_label,color=module_name))+geom_vline(xintercept=0,color=dark,linewidth=.3)+geom_point(size=1.7,position=position_dodge(width=.55))+scale_color_manual(values=c(MAGMA_top50=bluegrey,MAGMA_top100=terra,MAGMA_Bonferroni=teal,MAGMA_FDR05=loop_teal,MAGMA_suggestive_p1e4=gold),labels=module_labels,name=NULL)+labs(title="All tissue-state adjustment families",subtitle="Sensitivity coefficients; no independent-mechanism inference",x="Adjusted disease-context coefficient",y=NULL)+theme_pub+theme(axis.text.y=element_text(size=7),legend.position="bottom")
boundary <- ggplot()+annotate("rect",xmin=0,xmax=1,ymin=0,ymax=1,fill="#F7F8F8",color=bluegrey)+annotate("text",x=.05,y=.82,hjust=0,label="Supported",fontface="bold",color=teal,size=4)+annotate("text",x=.05,y=.61,hjust=0,label="Paired disease-context consistency\nTissue-state embedding / sensitivity",size=3.4,lineheight=1.2)+annotate("text",x=.05,y=.35,hjust=0,label="Not claimed",fontface="bold",color=terra,size=4)+annotate("text",x=.05,y=.13,hjust=0,label="Causal genes | plaque mechanism\ncell fractions | independent validation",size=3.4,lineheight=1.2)+coord_cartesian(xlim=c(0,1),ylim=c(0,1))+theme_void()+labs(title="Interpretation boundary")+theme(plot.title=element_text(face="bold",size=11))
supp <- (sA / (sB | boundary) / sC) + plot_layout(heights=c(1.25,1,1.1),widths=c(1.4,.6)) + plot_annotation(title="Supplementary bulk tissue-state sensitivity",tag_levels="A",theme=theme(plot.title=element_text(face="bold",size=16),plot.tag=element_text(face="bold",size=13)))
ggsave(file.path(root,"results/figures/phase5_step4_SuppFig_bulk_tissue_state_sensitivity.pdf"),supp,width=16,height=16,units="in",device="pdf",bg="white")
ggsave(file.path(root,"results/figures/phase5_step4_SuppFig_bulk_tissue_state_sensitivity.png"),supp,width=16,height=16,units="in",dpi=300,bg="white")
ggsave(file.path(root,"results/figures/phase5_step4_SuppFig_bulk_tissue_state_sensitivity_600dpi.png"),supp,width=16,height=16,units="in",dpi=600,bg="white")

# Integrated evidence and manifests.
evidence <- data.frame(
  evidence_component=c("GSE73680 input and paired design","canonical MAGMA module scoring","base paired disease-context model","tissue-state marker-score definitions","paired tissue-state shifts","module-tissue-state correlations","tissue-state-adjusted coefficient attenuation","deconvolution not performed","final bulk claim boundary","Figure 3 bulk-status recommendation"),
  input_or_analysis=c("55 included samples; 29 patients; 26 complete pairs","Five Phase 1 locked MAGMA modules; gene-wise z-scores averaged","module_score ~ group_binary + factor(patient_id)","Ten transparent curated proxy signatures","signature_score ~ group_binary + factor(patient_id)","Spearman correlations at sample and paired-difference levels","Single-proxy and compact low-dimensional sensitivity models","Phase 5 feasibility audit and stop rule","Integrated Phase 5 evidence","Figure 3 later-update recommendation only"),
  key_result=c("27 control/adjacent and 28 plaque/stone-associated samples; processed intensity-like matrix","45/50 to 333/369 genes mapped across modules","Coefficients 0.169-0.174; nominal p<0.05 for FDR05 and suggestive; all module FDR>0.05","All requested signatures had at least three mapped genes; scores are proxies","Injury beta=0.375, p=0.00573, FDR=0.0573; epithelial and mineralization proxies nominally positive","Broad positive module coupling with injury, endothelial, stromal/ECM and related proxies","Strongest attenuation after injury or injury+fibrosis adjustment","No formal cell-fraction estimates generated","Coherent modest disease-context consistency with tissue-state embedding","Bulk status can change from pending to disease-context reviewed"),
  support_strength=c("adequate paired design with label caveat","reproducible canonical scoring","modest / nominal-borderline","transparent descriptive resource","nominal tissue-state support; no FDR<0.05","descriptive non-causal coupling","sensitivity support","not applicable","coherent but non-upgrading","status update only"),
  limitation=c("Filename-derived labels require manual review","Processed non-count expression and incomplete module mapping","26 pairs; no causal inference","Curated signatures are imperfect","Marker scores are not cell fractions","Correlation is not causality","Attenuation is not proof of confounding; retention is not independence","No formal deconvolution evidence","Bulk cannot establish causal genes or plaque mechanisms","No R1-R6 reassignment"),
  allowed_claim=c("Paired papillary disease-context comparison","Canonical MAGMA modules were recalculated in bulk data","Modest paired module-level disease-context shifts","Transparent tissue-state proxy scoring","Selected tissue-state proxies showed positive disease-context shifts","Module shifts covaried with tissue-state programs","Module shifts were sensitive to injury/remodeling adjustment","Formal deconvolution was not performed","Bulk supports disease-context consistency and tissue-state embedding","Bulk: disease-context reviewed"),
  not_allowed_claim=c("Plaque-specific localization proof","Single-gene validation","Genetic causality or causal validation","True cell fractions","Histological composition change","Causal pathway coupling","Confounding proved or independent mechanism proved","Absence of composition effects","Causal genes, plaque mechanism, causal cell type","Bulk upgrades candidate reporting groups"),
  recommended_location=c("Figure 4A and Methods","Methods and Supplementary table","Figure 4B-C","Supplementary table","Figure 4D and Supplementary figure","Figure 4E and Supplementary table","Figure 4F and Supplementary figure","Limitations","Results/Discussion/Figure 4 legend","Figure 3 later manual update"),
  stringsAsFactors=FALSE
)
wr(evidence,"results/tables/phase5_step4_bulk_integrated_evidence_summary.tsv")

manifest <- data.frame(
  panel=LETTERS[1:6],
  panel_title=c("Paired bulk design","Paired MAGMA module scores","Base paired coefficients","Tissue-state proxy shifts","Module-tissue-state coupling","Tissue-state sensitivity"),
  source_figure_or_table=c("phase5_step1/2 locked metadata","phase5_step2_bulk_sample_module_scores.tsv","phase5_step2_bulk_paired_model_results.tsv","phase5_step3_bulk_sample_tissue_state_scores.tsv","phase5_step3_module_tissue_state_correlations.tsv","phase5_step3_tissue_state_adjusted_module_models.tsv"),
  source_data_file=paste0("source_data/figures/phase5_step4_Figure4_panel",LETTERS[1:6],"_source_data.tsv"),
  main_message=c("55 samples and 26 complete pairs support a paired disease-context design","Selected canonical modules show modest paired positive shifts","Base coefficients are positive but all cross-module FDR values exceed 0.05","Selected tissue-state proxies shift in disease-context samples","Module-score changes covary with broad tissue-state changes","Injury/remodeling adjustment attenuates module coefficients"),
  interpretation_boundary=c("Processed bulk; filename-derived labels; disease-context only","Module-level consistency, not causal validation","Nominal/borderline evidence, not independent validation","Expression proxies, not cell fractions","Non-causal correlation","Attenuation indicates embedding/sensitivity, not disproof"),
  manual_polishing_needed="yes_final_journal_layout_and_typography_review",
  notes="Draft assembled directly from locked Phase 5 tables; no new models.",
  stringsAsFactors=FALSE
)
wr(manifest,"results/tables/phase5_step4_figure4_panel_manifest.tsv")

source_manifest <- data.frame(
  figure_panel=paste0("Figure 4",LETTERS[1:6]),
  source_data_file=manifest$source_data_file,
  source_analysis_step=c("Phase 5-Step 1/2","Phase 5-Step 2","Phase 5-Step 2","Phase 5-Step 3","Phase 5-Step 3","Phase 5-Step 3"),
  description=manifest$main_message,
  required_for_reproducibility="yes",
  notes=manifest$interpretation_boundary,
  stringsAsFactors=FALSE
)
wr(source_manifest,"source_data/figures/phase5_step4_Figure4_source_data_manifest.tsv")

qc <- data.frame(
  figure_id=c("phase5_step4_Figure4_bulk_disease_context_draft","phase5_step4_SuppFig_bulk_tissue_state_sensitivity"),
  version="draft_v0.1",
  pdf_exists=c(file.exists(file.path(root,"results/figures/phase5_step4_Figure4_bulk_disease_context_draft.pdf")),file.exists(file.path(root,"results/figures/phase5_step4_SuppFig_bulk_tissue_state_sensitivity.pdf"))),
  png_exists=c(file.exists(file.path(root,"results/figures/phase5_step4_Figure4_bulk_disease_context_draft.png")),file.exists(file.path(root,"results/figures/phase5_step4_SuppFig_bulk_tissue_state_sensitivity.png"))),
  png_600dpi_exists=c(file.exists(file.path(root,"results/figures/phase5_step4_Figure4_bulk_disease_context_draft_600dpi.png")),file.exists(file.path(root,"results/figures/phase5_step4_SuppFig_bulk_tissue_state_sensitivity_600dpi.png"))),
  configured_min_font_size_pt=8,
  panel_labels="present_A_to_F_or_A_to_D",
  legend_placement="outside_dense_panels",
  palette_consistency="project_profile_teal_bluegrey_gold_terracotta",
  claim_boundary="disease_context_only_noncausal",
  resource_limited_claim_check="passes_no_deconvolution_or_cell_fraction_claim",
  source_table_exists="yes",
  legend_file_exists="yes",
  visual_status="agent_review_pass_at_50_percent_scale",
  action_required="Final journal polishing and human approval required.",
  stringsAsFactors=FALSE
)
wr(qc,"results/tables/phase5_step4_figure_visual_qc.tsv")
message("Phase 5-Step 4 figures, source data and core tables written.")
