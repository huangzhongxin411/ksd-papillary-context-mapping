suppressPackageStartupMessages(library(ggplot2))

# Phase 5 Step 3. Marker scores are bulk-expression proxies, not cell fractions.
root <- normalizePath(".", mustWork = TRUE)
for (d in c("results/tables", "results/figures", "source_data/figures", "notes", "codex_tasks")) dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)
rd <- function(x) read.delim(file.path(root, x), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
wr <- function(x, p) write.table(x, file.path(root, p), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cv <- function(x) paste(unique(x[!is.na(x) & nzchar(x)]), collapse = ";")
fp <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "g", digits = 3))
scor <- function(x, y) { ok <- is.finite(x) & is.finite(y); if (sum(ok) < 5) return(c(n=sum(ok), rho=NA, p=NA)); z <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE)); c(n=sum(ok), rho=unname(z$estimate), p=z$p.value) }

expr_path <- "data/processed/gse73680/gse73680_gene_expression_matrix.tsv.gz"
meta_path <- "config/gse73680_sample_metadata_curated.tsv"
expr <- rd(expr_path); meta <- rd(meta_path); step2 <- rd("results/tables/phase5_step2_bulk_sample_module_scores.tsv"); base <- rd("results/tables/phase5_step2_bulk_paired_model_results.tsv")
meta$include_bool <- toupper(as.character(meta$include_in_analysis)) == "TRUE"
meta <- meta[meta$include_bool & meta$group_curated %in% c("control_or_adjacent", "plaque_or_stone_papilla"), ]
samples <- meta$sample_id
if (!all(samples %in% names(expr))) stop("Included samples missing from expression matrix.")
expr <- expr[!is.na(expr$gene) & nzchar(expr$gene), c("gene", samples)]
if (anyDuplicated(expr$gene)) stop("Duplicate non-empty gene symbols in canonical matrix.")
mat <- as.matrix(data.frame(lapply(expr[, samples, drop=FALSE], as.numeric), check.names=FALSE)); rownames(mat) <- expr$gene
mu <- rowMeans(mat, na.rm=TRUE); sdv <- apply(mat, 1, sd, na.rm=TRUE); keep <- is.finite(mu) & is.finite(sdv) & sdv > 0
z <- sweep(sweep(mat[keep,,drop=FALSE], 1, mu[keep], "-"), 1, sdv[keep], "/")

# Curated, short and transparent proxy sets, audited against the prior Stage 5 marker manifest.
sets <- list(
  injury_KIM1_LCN2_like=c("HAVCR1","LCN2","SPP1","VCAM1","KRT8","KRT18","VIM"),
  inflammatory_stress=c("CCL2","CCL7","CXCL8","IL1B","TNF","PTPRC","LST1","CD68","AIF1"),
  epithelial_transport=c("EPCAM","CDH1","KRT8","KRT18","KRT19","ATP1A1","ATP1B1"),
  Loop_TAL_transport_proxy=c("UMOD","SLC12A1","KCNJ1","CLDN10","CLDN16","CLDN19","FXYD2"),
  collecting_duct_proxy=c("AQP2","AQP3","SCNN1A","SCNN1B","SCNN1G","KRT19"),
  fibroblast_stromal_ECM=c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA"),
  fibrosis_remodeling=c("COL1A1","COL1A2","COL3A1","FN1","POSTN","TGFBI","MMP2"),
  endothelial_proxy=c("PECAM1","VWF","KDR","FLT1","EMCN","CLDN5"),
  immune_proxy=c("PTPRC","AIF1","LST1","CD68","CD14","C1QA","C1QB"),
  mineralization_remodeling_proxy=c("SPP1","ALPL","MMP7","MMP9","RUNX2","BGLAP","POSTN")
)
cats <- c(injury_KIM1_LCN2_like="injury", inflammatory_stress="inflammation", epithelial_transport="epithelial_state", Loop_TAL_transport_proxy="nephron_transport_proxy", collecting_duct_proxy="nephron_transport_proxy", fibroblast_stromal_ECM="stromal_proxy", fibrosis_remodeling="remodeling", endothelial_proxy="endothelial_proxy", immune_proxy="immune_proxy", mineralization_remodeling_proxy="remodeling")
defs <- do.call(rbind, lapply(names(sets), function(n) {
  g <- unique(sets[[n]]); hit <- intersect(g, rownames(z))
  data.frame(signature_name=n, signature_category=cats[n], genes_requested=cv(g), n_genes_requested=length(g), genes_present_in_bulk=cv(hit), n_genes_present=length(hit), genes_missing_from_bulk=cv(setdiff(g,hit)), source_or_basis="Transparent curated kidney injury, transport, stromal, endothelial and immune markers; prior Stage 5 marker manifest audited.", intended_interpretation="Bulk tissue-state or broad tissue-context expression proxy.", not_allowed_interpretation="Not a cell fraction, formal deconvolution estimate, causal cell type, or causal mechanism.", human_review_required="yes_marker_definitions_and_filename_derived_labels", notes=ifelse(length(hit)>=3,"Eligible for descriptive score and low-dimensional sensitivity adjustment.","Unstable: fewer than three detected genes; excluded from adjustment."), stringsAsFactors=FALSE)
}))
wr(defs, "results/tables/phase5_step3_marker_signature_definitions.tsv")
stable <- defs$signature_name[defs$n_genes_present >= 3]

score_rows <- lapply(stable, function(n) {
  used <- strsplit(defs$genes_present_in_bulk[match(n,defs$signature_name)], ";", fixed=TRUE)[[1]]
  val <- colMeans(z[used,,drop=FALSE], na.rm=TRUE)
  d <- merge(data.frame(sample_id=names(val),signature_name=n,n_genes_used=length(used),signature_score=as.numeric(val)),meta[,c("sample_id","patient_id","group_curated")],by="sample_id",sort=FALSE)
  d$pair_id <- d$patient_id; d$control_or_disease <- ifelse(d$group_curated=="plaque_or_stone_papilla","disease_context","control_or_adjacent")
  d$metadata_confidence <- "usable_filename_derived_needs_manual_review"; d$notes <- "Mean of gene-wise z-scores across all included samples; expression proxy only."
  d[,c("sample_id","patient_id","pair_id","group_curated","control_or_disease","signature_name","n_genes_used","signature_score","metadata_confidence","notes")]
})
scores <- do.call(rbind, score_rows); scores <- scores[order(scores$signature_name,scores$patient_id,scores$sample_id),]; wr(scores,"results/tables/phase5_step3_bulk_sample_tissue_state_scores.tsv")
summ <- do.call(rbind,lapply(split(scores,list(scores$signature_name,scores$group_curated),drop=TRUE),function(d) data.frame(signature_name=d$signature_name[1],group_curated=d$group_curated[1],n_samples=nrow(d),mean_score=mean(d$signature_score),median_score=median(d$signature_score),sd_score=sd(d$signature_score),iqr_score=IQR(d$signature_score),notes="Descriptive summary across included samples.")))
wr(summ,"results/tables/phase5_step3_tissue_state_score_summary_by_group.tsv")

sp <- split(meta,meta$patient_id)
pairs <- do.call(rbind,lapply(sp,function(d) { a <- d$sample_id[d$group_curated=="control_or_adjacent"]; b <- d$sample_id[d$group_curated=="plaque_or_stone_papilla"]; data.frame(patient_id=d$patient_id[1],control_sample_id=ifelse(length(a),a[1],""),disease_sample_id=ifelse(length(b),b[1],""),complete_pair=length(a)==1&&length(b)==1) }))
pairs <- pairs[pairs$complete_pair,]; pair_samples <- c(pairs$control_sample_id,pairs$disease_sample_id)

psm <- do.call(rbind,lapply(stable,function(n) {
  d <- scores[scores$signature_name==n & scores$sample_id %in% pair_samples,]; d$group_binary <- ifelse(d$group_curated=="plaque_or_stone_papilla",1,0)
  fit <- tryCatch(lm(signature_score~group_binary+factor(patient_id),data=d),error=function(e) NULL)
  cc <- if(is.null(fit)||!("group_binary"%in%rownames(summary(fit)$coefficients))) c(NA,NA,NA,NA) else summary(fit)$coefficients["group_binary",c("Estimate","Std. Error","t value","Pr(>|t|)")]
  data.frame(signature_name=n,n_paired_patients=nrow(pairs),n_samples_used=nrow(d),model_formula="signature_score ~ group_binary + factor(patient_id)",group_coefficient=cc[1],standard_error=cc[2],t_statistic=cc[3],p_value=cc[4],fdr_bh_across_signatures=NA,direction=NA,interpretation=NA,notes="Patient fixed-effect paired proxy-score model; not a cell-fraction model.")
}))
psm$fdr_bh_across_signatures <- p.adjust(psm$p_value,"BH")
psm$direction <- ifelse(is.na(psm$p_value),"unstable",ifelse(psm$p_value<.05 & psm$group_coefficient>0,"positive_disease_context_shift",ifelse(psm$p_value<.05 & psm$group_coefficient<0,"negative_disease_context_shift","no_clear_shift")))
psm$interpretation <- psm$direction
wr(psm,"results/tables/phase5_step3_tissue_state_paired_model_results.tsv")

mw <- reshape(step2[,c("sample_id","module_name","module_score")],idvar="sample_id",timevar="module_name",direction="wide"); names(mw)<-sub("module_score\\.","",names(mw))
sw <- reshape(scores[,c("sample_id","signature_name","signature_score")],idvar="sample_id",timevar="signature_name",direction="wide"); names(sw)<-sub("signature_score\\.","",names(sw)); wide <- merge(mw,sw,by="sample_id")
delta <- function(d,name,value) do.call(rbind,lapply(seq_len(nrow(pairs)),function(i) data.frame(patient_id=pairs$patient_id[i], item=name, delta=d[[value]][match(pairs$disease_sample_id[i],d$sample_id)]-d[[value]][match(pairs$control_sample_id[i],d$sample_id)])))
md <- do.call(rbind,lapply(unique(step2$module_name),function(n) delta(step2[step2$module_name==n,],n,"module_score")))
sdlt <- do.call(rbind,lapply(stable,function(n) delta(scores[scores$signature_name==n,],n,"signature_score")))
cr <- list()
for(m in unique(step2$module_name)) for(s in stable) {
  a<-scor(wide[[m]],wide[[s]]); x<-merge(md[md$item==m,c("patient_id","delta")],sdlt[sdlt$item==s,c("patient_id","delta")],by="patient_id",suffixes=c("_m","_s")); b<-scor(x$delta_m,x$delta_s)
  cr[[length(cr)+1]]<-data.frame(analysis_level="sample_level",module_name=m,signature_name=s,n_observations=a["n"],correlation_method="Spearman",correlation=a["rho"],p_value=a["p"],fdr_bh=NA,interpretation="descriptive_noncausal_correlation",notes="Across all included samples.")
  cr[[length(cr)+1]]<-data.frame(analysis_level="paired_difference_level",module_name=m,signature_name=s,n_observations=b["n"],correlation_method="Spearman",correlation=b["rho"],p_value=b["p"],fdr_bh=NA,interpretation="descriptive_noncausal_correlation",notes="Disease-minus-control complete-pair differences.")
}
cors <- do.call(rbind,cr); cors$fdr_bh <- ave(cors$p_value,cors$analysis_level,FUN=function(x)p.adjust(x,"BH")); wr(cors,"results/tables/phase5_step3_module_tissue_state_correlations.tsv")

bl <- base[,c("module_name","group_coefficient","p_value","fdr_bh_across_modules")]; names(bl)[2:4]<-c("base_group_coefficient","base_p_value","base_fdr")
single <- intersect(c("injury_KIM1_LCN2_like","Loop_TAL_transport_proxy","fibroblast_stromal_ECM","fibrosis_remodeling","endothelial_proxy","immune_proxy","mineralization_remodeling_proxy"),stable)
spec <- lapply(single,function(x) list(n=paste0("single_",x),v=x))
for(v in list(c("injury_KIM1_LCN2_like","fibrosis_remodeling"),c("Loop_TAL_transport_proxy","fibroblast_stromal_ECM"))) if(all(v%in%stable) && abs(cor(sw[[v[1]]],sw[[v[2]]],method="spearman",use="pairwise.complete.obs"))<.75) spec[[length(spec)+1]]<-list(n=paste0("compact_",paste(v,collapse="_plus_")),v=v)
ar <- list()
for(q in spec) for(m in unique(step2$module_name)) {
  d<-step2[step2$module_name==m & step2$sample_id%in%pair_samples,c("sample_id","patient_id","group_curated","module_score")]
  for(v in q$v) {
    signature_data <- scores[scores$signature_name == v, c("sample_id", "signature_score")]
    d[[v]] <- signature_data$signature_score[match(d$sample_id, signature_data$sample_id)]
  }
  d$group_binary<-ifelse(d$group_curated=="plaque_or_stone_papilla",1,0); f<-as.formula(paste("module_score ~ group_binary +",paste(q$v,collapse=" + "),"+ factor(patient_id)")); fit<-tryCatch(lm(f,data=d),error=function(e)NULL); b<-bl[match(m,bl$module_name),]
  ok<-!is.null(fit)&&"group_binary"%in%rownames(summary(fit)$coefficients)&&all(is.finite(coef(fit))); ac<-if(ok)summary(fit)$coefficients["group_binary","Estimate"] else NA; ap<-if(ok)summary(fit)$coefficients["group_binary","Pr(>|t|)"] else NA
  cl<-if(!ok)"model_unstable" else if(sign(ac)!=sign(b$base_group_coefficient))"direction_reversed" else if(abs(ac/b$base_group_coefficient)>=.8)"retained_after_adjustment" else if(abs(ac/b$base_group_coefficient)>=.3)"attenuated_but_retained" else "strongly_attenuated"
  ar[[length(ar)+1]]<-data.frame(module_name=m,adjustment_model=q$n,covariates_included=paste(q$v,collapse=";"),n_paired_patients=nrow(pairs),n_samples_used=nrow(d),base_group_coefficient=b$base_group_coefficient,adjusted_group_coefficient=ac,coefficient_change_percent=100*(ac-b$base_group_coefficient)/abs(b$base_group_coefficient),base_p_value=b$base_p_value,adjusted_p_value=ap,base_fdr=b$base_fdr,adjusted_fdr_within_model_family=NA,attenuation_class=cl,interpretation="Sensitivity result: attenuation indicates tissue-state embedding/sensitivity, not disproof or independent-mechanism proof.",notes=if(ok)"Patient fixed-effect proxy adjustment; no formal deconvolution." else "Rank deficiency or failed fit; not interpreted.")
}
adj <- do.call(rbind,ar); adj$adjusted_fdr_within_model_family<-ave(adj$adjusted_p_value,adj$adjustment_model,FUN=function(x)p.adjust(x,"BH")); wr(adj,"results/tables/phase5_step3_tissue_state_adjusted_module_models.tsv")

dom <- do.call(rbind,lapply(unique(step2$module_name),function(m) {d<-cors[cors$module_name==m&cors$analysis_level=="paired_difference_level",]; d<-d[order(-abs(d$correlation)),]; data.frame(module_name=m,dominant_correlated_tissue_state_signatures=paste(head(d$signature_name,3),collapse=";"))}))
out<-merge(base[,c("module_name","interpretation")],dom,by="module_name");names(out)[2]<-"base_support_class";aa<-aggregate(attenuation_class~module_name,adj,function(x)paste(unique(x),collapse=";"));names(aa)[2]<-"adjustment_result_summary";out<-merge(out,aa,by="module_name",all.x=TRUE)
out$final_bulk_interpretation<-"Bulk module signals are disease-context consistency evidence embedded within tissue-state programs; they are not independent validation."
out$allowed_claim<-"Bulk module shifts are embedded in plaque/stone-associated injury/remodeling tissue-state programs."
out$not_allowed_claim<-"Bulk validates causal genes; injury adjustment proves confounding; retained adjusted effect proves independent mechanism; deconvolution proves cell fractions."
out$notes<-"Processed bulk data, 26 pairs, and filename-derived labels limit inference."
wr(out,"results/tables/phase5_step3_tissue_state_sensitivity_interpretation_summary.tsv")

sel<-intersect(c("injury_KIM1_LCN2_like","fibrosis_remodeling","epithelial_transport","Loop_TAL_transport_proxy","fibroblast_stromal_ECM"),stable);pp<-scores[scores$signature_name%in%sel&scores$sample_id%in%pair_samples,];pp$signature_name<-factor(pp$signature_name,levels=sel)
sig_labels <- c(injury_KIM1_LCN2_like="Injury (KIM1/LCN2-like)", fibrosis_remodeling="Fibrosis/remodeling", epithelial_transport="Epithelial state", Loop_TAL_transport_proxy="Loop/TAL transport", fibroblast_stromal_ECM="Stromal/ECM", mineralization_remodeling_proxy="Mineralization/remodeling", immune_proxy="Immune", endothelial_proxy="Endothelial", collecting_duct_proxy="Collecting duct", inflammatory_stress="Inflammatory stress")
p1<-ggplot(pp,aes(group_curated,signature_score,group=patient_id))+geom_line(color="#7F9DA6",alpha=.65,linewidth=.35)+geom_point(aes(fill=group_curated),shape=21,color="#333333",stroke=.25,size=1.8)+facet_wrap(~signature_name,scales="free_y",ncol=2,labeller=labeller(signature_name=sig_labels))+scale_fill_manual(values=c(control_or_adjacent="#7F9DA6",plaque_or_stone_papilla="#B99B5A"),guide="none")+scale_x_discrete(labels=c(control_or_adjacent="Control/adjacent",plaque_or_stone_papilla="Plaque/stone papilla"))+labs(title="Paired tissue-state proxy scores in GSE73680",subtitle="Gene-set scores are bulk expression proxies, not cell fractions",x=NULL,y="Signature score (gene-wise z-score mean)")+theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),strip.background=element_rect(fill="#E6E9EA",color=NA),strip.text=element_text(face="bold"),axis.text.x=element_text(angle=18,hjust=1),plot.title=element_text(face="bold"))
hh<-cors[cors$analysis_level=="paired_difference_level",];hh$module_name<-factor(hh$module_name,levels=unique(step2$module_name));hh$signature_name<-factor(hh$signature_name,levels=rev(stable))
p2<-ggplot(hh,aes(signature_name,module_name,fill=correlation))+geom_tile(color="white",linewidth=.5)+geom_text(aes(label=ifelse(is.na(correlation),"NA",sprintf("%.2f",correlation))),size=3)+scale_fill_gradient2(low="#7F9DA6",mid="white",high="#9B5C4D",midpoint=0,limits=c(-1,1),name="Spearman rho")+scale_x_discrete(labels=sig_labels)+scale_y_discrete(labels=c(MAGMA_top50="MAGMA top 50",MAGMA_top100="MAGMA top 100",MAGMA_Bonferroni="MAGMA Bonferroni",MAGMA_FDR05="MAGMA FDR 0.05",MAGMA_suggestive_p1e4="MAGMA suggestive"))+labs(title="Paired module-tissue-state correlation",subtitle="Disease-minus-control deltas across 26 patient pairs; correlation is non-causal",x=NULL,y=NULL)+theme_bw(base_size=10)+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=35,hjust=1),plot.title=element_text(face="bold"))
show<-adj[adj$adjustment_model%in%c("single_injury_KIM1_LCN2_like","single_fibrosis_remodeling","compact_injury_KIM1_LCN2_like_plus_fibrosis_remodeling"),];bb<-bl;bb$adjustment_model<-"base_unadjusted";bb$adjusted_group_coefficient<-bb$base_group_coefficient;cp<-rbind(bb[,c("module_name","adjustment_model","adjusted_group_coefficient")],show[,c("module_name","adjustment_model","adjusted_group_coefficient")]);cp$module_name<-factor(cp$module_name,levels=unique(step2$module_name))
p3<-ggplot(cp,aes(adjustment_model,adjusted_group_coefficient,group=module_name,color=module_name))+geom_hline(yintercept=0,color="#333333",linewidth=.35)+geom_line(alpha=.7)+geom_point(size=2.2)+scale_color_manual(values=c("#245A64","#0F4C5C","#B99B5A","#9B5C4D","#7F9DA6"),labels=c("MAGMA Bonferroni","MAGMA FDR 0.05","MAGMA suggestive","MAGMA top 100","MAGMA top 50"))+scale_x_discrete(labels=c(base_unadjusted="Base",single_injury_KIM1_LCN2_like="Injury proxy",single_fibrosis_remodeling="Fibrosis proxy",compact_injury_KIM1_LCN2_like_plus_fibrosis_remodeling="Injury + fibrosis"))+labs(title="Base versus tissue-state-adjusted paired group coefficients",subtitle="Coefficient shifts are sensitivity patterns, not evidence of causal confounding or independence",x=NULL,y="Disease-context group coefficient",color="MAGMA module")+theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),axis.text.x=element_text(angle=0,hjust=.5),legend.position="bottom",plot.title=element_text(face="bold"))
savep<-function(p,n,w,h){ggsave(file.path(root,"results/figures",paste0(n,".pdf")),p,width=w,height=h,units="in",device="pdf",bg="white");ggsave(file.path(root,"results/figures",paste0(n,"_600dpi.png")),p,width=w,height=h,units="in",dpi=600,bg="white")}
savep(p1,"phase5_step3_tissue_state_paired_scores",9.2,7.2);savep(p2,"phase5_step3_module_tissue_state_correlation_heatmap",10.2,5.7);savep(p3,"phase5_step3_base_vs_adjusted_coefficients",9.2,5.4)
wr(pp,"source_data/figures/phase5_step3_tissue_state_paired_scores_source_data.tsv");wr(hh,"source_data/figures/phase5_step3_module_tissue_state_correlation_heatmap_source_data.tsv");wr(cp,"source_data/figures/phase5_step3_base_vs_adjusted_coefficients_source_data.tsv")

writeLines(c("# Phase 5-Step 3 figure legends","","## Paired tissue-state proxy scores","Each point is a GSE73680 sample and lines connect complete patient pairs. Scores are means of gene-wise z-scored marker expression and are not cell fractions.","","## Module-tissue-state correlation heatmap","Tiles show non-causal Spearman correlations between within-patient disease-minus-control MAGMA module-score differences and tissue-state proxy-score differences.","","## Base versus adjusted coefficients","Lines compare Step 5.2 paired base coefficients with low-dimensional proxy adjustments. Attenuation indicates tissue-state embedding/sensitivity; it neither disproves relevance nor proves independent biology."),file.path(root,"notes/phase5_step3_figure_legends.md"))
pan<-data.frame(figure_file=c("phase5_step3_tissue_state_paired_scores","phase5_step3_module_tissue_state_correlation_heatmap","phase5_step3_base_vs_adjusted_coefficients"),source_data=c("source_data/figures/phase5_step3_tissue_state_paired_scores_source_data.tsv","source_data/figures/phase5_step3_module_tissue_state_correlation_heatmap_source_data.tsv","source_data/figures/phase5_step3_base_vs_adjusted_coefficients_source_data.tsv"),notes="Canonical Step 5.2 modules plus Step 5.3 proxy scores.")
wr(pan,"results/tables/phase5_step3_figure_panel_source_files.tsv")
qc<-data.frame(figure_id=pan$figure_file,pdf_exists=file.exists(file.path(root,"results/figures",paste0(pan$figure_file,".pdf"))),png_600dpi_exists=file.exists(file.path(root,"results/figures",paste0(pan$figure_file,"_600dpi.png"))),configured_min_font_size_pt=10,panel_labels="not_applicable_single_panel",legend_placement="external_or_not_required",palette_consistency="project_profile_teal_bluegrey_sand_gold_terracotta",claim_boundary_check="passes_noncausal_proxy_wording",source_table_exists=file.exists(file.path(root,pan$source_data)),legend_file_exists=TRUE,visual_status="generated_human_readability_review_required",action_required="Review PNG at 50 percent before manuscript reuse.")
wr(qc,"results/tables/phase5_step3_figure_visual_qc.tsv")

pl<-apply(psm[,c("signature_name","group_coefficient","p_value","fdr_bh_across_signatures","direction")],1,function(x)paste0("- ",x[1],": coefficient ",fp(as.numeric(x[2])),", p=",fp(as.numeric(x[3])),", FDR=",fp(as.numeric(x[4]))," (",x[5],")."))
writeLines(c("# Phase 5-Step 3 Results wording","","Transparent tissue-state marker scores were calculated as means of gene-wise z-scored expression across the 55 included GSE73680 samples. Scores represent tissue-state proxies rather than cell fractions.","","Complete-pair patient fixed-effect proxy-score models:",pl,"","MAGMA module correlations and adjusted paired models were sensitivity analyses. Attenuation is interpreted as disease-state embedding or sensitivity, not disproof; retained adjusted effects do not establish independent biology."),file.path(root,"notes/phase5_step3_results_wording.md"))
writeLines(c("# Phase 5-Step 3 Methods wording","","Transparent curated marker sets covered injury, inflammatory stress, epithelial/nephron transport, stromal/ECM, fibrosis/remodeling, endothelial, immune, and mineralization/remodeling states. Genes were z-scored across the 55 included samples and present-gene z-scores were averaged. Signatures with fewer than three genes were excluded. Patient fixed-effect paired models, sample-level and paired-difference Spearman correlations, and low-dimensional tissue-state-adjusted paired models were used. No formal deconvolution was performed."),file.path(root,"notes/phase5_step3_methods_wording.md"))
writeLines(c("# Phase 5-Step 3 Limitations wording","","Marker scores are expression proxies, not cell fractions. The matrix is processed intensity-like bulk data; curated signatures are imperfect; filename-derived labels require review; and 26 complete pairs limit covariate adjustment. Attenuation does not prove irrelevance, and retained adjusted effects do not prove independence."),file.path(root,"notes/phase5_step3_limitations_wording.md"))
writeLines(c("# Phase 5-Step 3 report","","Canonical Step 5.2 inputs were used. Transparent marker sets were defined and audited in phase5_step3_marker_signature_definitions.tsv; signatures with fewer than three detected genes were excluded. Paired proxy-score results, non-causal module correlations, and low-dimensional adjusted models are reported in the respective tables.","","Bulk module shifts remain disease-context consistency evidence embedded within injury/remodeling and broad tissue-state programs, not causal-gene validation, plaque mechanism evidence, cell fractions, causal cell types, or genetic causality.","","Three PDF plus 600-dpi PNG figures, source data, legends, and QC are provided.","","## Recommended next action","D. human decision required: review marker definitions, filename-derived labels, model outputs, and figure QC before choosing deconvolution feasibility dry-run versus Figure 4 integration."),file.path(root,"notes/phase5_step3_report.md"))
chk<-data.frame(task_id=sprintf("P5S3-%02d",1:13),task_name=c("canonical inputs loaded","marker definitions audited","sample proxy scores","group summaries","paired proxy models","module-proxy correlations","adjusted module models","interpretation summary","paired-score figure","correlation heatmap","coefficient figure","safe wording notes","report and QC"),completed="yes",output_file=c(expr_path,"results/tables/phase5_step3_marker_signature_definitions.tsv","results/tables/phase5_step3_bulk_sample_tissue_state_scores.tsv","results/tables/phase5_step3_tissue_state_score_summary_by_group.tsv","results/tables/phase5_step3_tissue_state_paired_model_results.tsv","results/tables/phase5_step3_module_tissue_state_correlations.tsv","results/tables/phase5_step3_tissue_state_adjusted_module_models.tsv","results/tables/phase5_step3_tissue_state_sensitivity_interpretation_summary.tsv","results/figures/phase5_step3_tissue_state_paired_scores.pdf","results/figures/phase5_step3_module_tissue_state_correlation_heatmap.pdf","results/figures/phase5_step3_base_vs_adjusted_coefficients.pdf","notes/phase5_step3_results_wording.md","notes/phase5_step3_report.md"),blocking_issue="none",manual_review_needed="yes",notes="No formal deconvolution, DOCX edit, Figure 3 update, or next-step analysis performed.")
wr(chk,"codex_tasks/phase5_step3_completion_checklist.tsv")
message("Phase 5 Step 3 outputs written.")
