#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(svglite); library(ragg)})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
fig_dir <- file.path(root, "results/figures/revision/stage7E_figure3_data_enriched_v0.4")
table_dir <- file.path(root, "results/tables/revision/stage7E_figure3_data_enriched_v0.4")
log_dir <- file.path(root, "logs/revision/stage7E_figure3_data_enriched_v0.4")
for (x in c(fig_dir, table_dir, log_dir)) dir.create(x, recursive = TRUE, showWarnings = FALSE)

paths <- c(model="results/tables/revision/stage3R_gene_tiering/candidate_gene_evidence_model_v0.2.tsv",
  decision="results/tables/revision/stage3R_gene_tiering/candidate_gene_decision_log_v0.2.tsv",
  counts="results/tables/revision/stage3R_gene_tiering/evidence_model_summary_counts_v0.2.tsv",
  exemplars="results/tables/revision/stage3R_gene_tiering/curated_exemplar_panel_v0.2.tsv",
  twas="results/tables/revision/stage2_genetic/twas_output_audit.tsv",
  one_snp="results/tables/revision/stage2_genetic/twas_one_snp_model_audit.tsv",
  smr="results/tables/revision/stage2_genetic/smr_coloc_feasibility.tsv")
full_paths <- setNames(file.path(root, unname(paths)), names(paths))
if (!all(file.exists(full_paths))) stop("Missing frozen Figure 3 source files")
model <- fread(full_paths[["model"]]); counts <- fread(full_paths[["counts"]])
exemplars <- fread(full_paths[["exemplars"]]); twas <- fread(full_paths[["twas"]])
one_snp <- fread(full_paths[["one_snp"]]); smr <- fread(full_paths[["smr"]])

expected_groups <- c(R1_MAGMA_Bonferroni_only=79L, R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy=1L,
  R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy=13L, R4_MAGMA_lower_priority=33L,
  R5_TWAS_proxy_only=4L, R6_contextual_or_unresolved=0L)
obs <- counts$n_unique_genes[match(names(expected_groups), counts$reporting_group)]
if (!identical(as.integer(obs), as.integer(expected_groups))) stop("Reporting-group counts differ from lock")
if (twas$n_genes_tested[1]!=5989 || twas$n_fdr_significant[1]!=51 || twas$n_one_snp_fdr_significant[1]!=42 || twas$n_multi_snp_fdr_significant[1]!=9) stop("TWAS counts differ from lock")

font_family <- "Helvetica"
pal <- c(dark="#2D3436", grey="#657174", mid="#AAB5B7", line="#DDE3E4", light="#F6F8F8",
  teal="#245A64", teal2="#0F6B73", teal_mid="#7F9DA6", teal_pale="#DCEBEC",
  gold="#B99B5A", gold_dark="#725E28", gold_pale="#F3EBD8", coral="#A76655", coral_pale="#F3E3DF", white="#FFFFFF")
theme_pub <- function(base_size=7.2) theme_classic(base_size=base_size, base_family=font_family)+
  theme(text=element_text(family=font_family, colour=pal[["dark"]]), plot.title=element_text(face="bold", size=base_size+1, margin=margin(b=2)),
    plot.subtitle=element_text(size=base_size-.2, colour=pal[["grey"]], margin=margin(b=3)), axis.title=element_text(size=base_size-.1),
    axis.text=element_text(size=base_size-.5), axis.line=element_line(linewidth=.32), axis.ticks=element_line(linewidth=.32),
    panel.grid=element_blank(), plot.margin=margin(4,4,4,4))
tag_panel <- function(p, tag) p+labs(tag=tag)+theme(plot.tag=element_text(family=font_family, face="bold", size=10), plot.tag.position=c(0,1))

# A: source-backed two-axis count matrix.
cells <- data.table(xmin=c(0,1,0,1), xmax=c(1,2,1,2), ymin=c(0,0,1,1), ymax=c(1,1,2,2),
  label=c("33 / 0\nR4 / R6\nlower / unresolved","79\nR1\nMAGMA only","4\nR5\nTWAS proxy only","1 / 13\nR2 / R3\nMAGMA + proxy TWAS"),
  fill=c(pal[["light"]],pal[["teal_pale"]],pal[["gold_pale"]],"#E7E6D9"))
pA <- ggplot(cells)+geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill),colour=pal[["white"]],linewidth=1.1)+
  geom_text(aes((xmin+xmax)/2,(ymin+ymax)/2,label=label),family=font_family,fontface="bold",size=2.25,lineheight=.9)+
  scale_fill_identity()+annotate("text",x=1.98,y=2.13,label="gene count",hjust=1,family=font_family,size=1.55,colour=pal[["grey"]])+
  scale_x_continuous(breaks=NULL,limits=c(-.05,2.05),expand=c(0,0))+
  scale_y_continuous(breaks=NULL,limits=c(-.05,2.22),expand=c(0,0))+
  labs(title="Two-axis count matrix",x="MAGMA genetic priority",y="Kidney cortex TWAS proxy")+theme_pub(7.1)+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text=element_text(size=5.9),axis.title.x=element_text(face="bold",colour=pal[["teal"]]),
    axis.title.y=element_text(face="bold",colour=pal[["gold_dark"]]),plot.margin=margin(4,6,4,7))

# B: reporting-group counts, not a priority ladder.
group_labels <- c(R1_MAGMA_Bonferroni_only="R1  MAGMA only",
  R2_MAGMA_Bonferroni_plus_multi_snp_TWAS_proxy="R2  MAGMA + multi-SNP proxy",
  R3_MAGMA_Bonferroni_plus_one_snp_TWAS_proxy="R3  MAGMA + one-SNP proxy",
  R4_MAGMA_lower_priority="R4  Lower MAGMA priority",R5_TWAS_proxy_only="R5  TWAS proxy only",R6_contextual_or_unresolved="R6  Contextual / unresolved")
dB <- counts[reporting_group %chin% names(group_labels),.(reporting_group,n_unique_genes)]
dB[,label:=factor(unname(group_labels[reporting_group]),levels=rev(unname(group_labels)))]
dB[,class:=fifelse(grepl("R5|R6",as.character(label)),"context",fifelse(grepl("R2|R3",as.character(label)),"proxy","magma"))]
dB[,label_x:=fifelse(n_unique_genes==0,1.6,n_unique_genes+1.6)]
pB <- ggplot(dB,aes(n_unique_genes,label))+geom_col(aes(fill=class),width=.58)+
  geom_point(data=dB[n_unique_genes==0],aes(x=0,y=label),shape=21,size=2.4,fill=pal[["white"]],colour=pal[["mid"]],stroke=.55,inherit.aes=FALSE)+
  geom_text(aes(x=label_x,label=n_unique_genes),hjust=0,family=font_family,fontface="bold",size=2.1)+
  scale_fill_manual(values=c(magma=pal[["teal_mid"]],proxy=pal[["gold"]],context=pal[["mid"]]),guide="none")+
  scale_x_continuous(limits=c(0,88),breaks=c(0,40,80),expand=c(0,0))+
  labs(title="Reporting group distribution",subtitle="Mutually exclusive reporting groups; not causal tiers",x="Genes",y=NULL)+theme_pub(7.1)+
  theme(panel.grid.major.x=element_line(colour=pal[["line"]],linewidth=.25),axis.line.y=element_blank(),axis.ticks.y=element_blank(),
    axis.text.y=element_text(size=5.7),plot.margin=margin(4,6,4,6))

# C: visual TWAS proxy flow and one-SNP burden.
pC <- ggplot()+
  annotate("rect",xmin=.03,xmax=.25,ymin=.57,ymax=.82,fill=pal[["light"]],colour=pal[["line"]],linewidth=.35)+
  annotate("rect",xmin=.39,xmax=.61,ymin=.57,ymax=.82,fill=pal[["teal_pale"]],colour=pal[["teal_mid"]],linewidth=.35)+
  annotate("segment",x=.27,xend=.37,y=.695,yend=.695,colour=pal[["mid"]],linewidth=.55,arrow=grid::arrow(length=grid::unit(3,"pt"),type="closed"))+
  annotate("text",x=.14,y=.695,label="5,989\ntested",family=font_family,fontface="bold",size=2.35,lineheight=.88)+
  annotate("text",x=.50,y=.695,label="51\nFDR-supported",family=font_family,fontface="bold",size=2.2,lineheight=.88)+
  annotate("rect",xmin=.66,xmax=.91,ymin=.65,ymax=.79,fill=pal[["coral_pale"]],colour=pal[["white"]],linewidth=.3)+
  annotate("rect",xmin=.91,xmax=.985,ymin=.65,ymax=.79,fill=pal[["gold_pale"]],colour=pal[["white"]],linewidth=.3)+
  annotate("segment",x=.63,xend=.65,y=.695,yend=.695,colour=pal[["mid"]],linewidth=.55,arrow=grid::arrow(length=grid::unit(3,"pt"),type="closed"))+
  annotate("text",x=.785,y=.72,label="42 one-SNP",family=font_family,fontface="bold",size=1.75)+
  annotate("text",x=.948,y=.72,label="9",family=font_family,fontface="bold",size=1.75)+
  annotate("text",x=.66,y=.53,label="one-SNP caution  ·  supplementary proxy only",hjust=0,family=font_family,size=1.55,colour=pal[["coral"]])+
  annotate("text",x=.03,y=.22,label="GTEx kidney cortex proxy  ·  no papilla-specific eQTL",hjust=0,family=font_family,fontface="bold",size=1.7,colour=pal[["gold_dark"]])+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1),clip="off")+labs(title="TWAS proxy flow")+theme_void(base_family=font_family)+
  theme(plot.title=element_text(face="bold",size=8.1,margin=margin(b=2,l=13)),plot.margin=margin(4,7,4,5))

# D: curated exemplars as interpretive flags only.
exemplar_names <- c("UMOD","CASR","CLDN14","CLDN10","HIBADH","PKD2")
if (!setequal(exemplars$gene,exemplar_names)) stop("Exemplar set differs from lock")
chip_x<-c(.16,.49,.82,.16,.49,.82); chip_y<-c(.64,.64,.64,.34,.34,.34)
pD <- ggplot()+annotate("rect",xmin=.03,xmax=.97,ymin=.18,ymax=.80,fill=pal[["gold_pale"]],colour=pal[["gold"]],linewidth=.35)+
  annotate("text",x=.5,y=.90,label="Biology-informed flag  ·  No evidence upgrade",family=font_family,fontface="bold",size=1.9,colour=pal[["gold_dark"]])
for(i in seq_along(exemplar_names)) pD<-pD+annotate("rect",xmin=chip_x[i]-.125,xmax=chip_x[i]+.125,ymin=chip_y[i]-.10,ymax=chip_y[i]+.10,fill=pal[["white"]],colour="#D2C39E",linewidth=.28)+
  annotate("text",x=chip_x[i],y=chip_y[i],label=exemplar_names[i],family=font_family,fontface="bold.italic",size=2.05)
pD <- pD+annotate("text",x=.5,y=.08,label="Interpretive anchors only; not validated and not SMR/coloc-supported",family=font_family,size=1.55,colour=pal[["grey"]])+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1),clip="off")+labs(title="Curated exemplar boundary")+theme_void(base_family=font_family)+
  theme(plot.title=element_text(face="bold",size=8.1,margin=margin(b=2,l=13)),plot.margin=margin(4,5,4,5))

# E: compact paired claim bands.
pE<-ggplot()+
  annotate("rect",xmin=.02,xmax=.61,ymin=.22,ymax=.68,fill=pal[["teal_pale"]],colour=pal[["teal2"]],linewidth=.38)+
  annotate("rect",xmin=.64,xmax=.985,ymin=.22,ymax=.68,fill=pal[["coral_pale"]],colour=pal[["coral"]],linewidth=.38)+
  annotate("text",x=.035,y=.80,label="SUPPORTED REPORTING",hjust=0,family=font_family,fontface="bold",size=1.9,colour=pal[["teal"]])+
  annotate("text",x=.655,y=.80,label="NOT CLAIMED",hjust=0,family=font_family,fontface="bold",size=1.9,colour=pal[["coral"]])+
  annotate("text",x=.315,y=.45,label="Transparent reporting  ·  MAGMA priority  ·  kidney cortex proxy context  ·  interpretive exemplars",family=font_family,fontface="bold",size=1.52)+
  annotate("text",x=.8125,y=.45,label="Causal genes  ·  validated genes\nSMR/coloc support  ·  therapeutic targets",family=font_family,fontface="bold",size=1.48,lineheight=.9)+
  coord_cartesian(xlim=c(0,1.025),ylim=c(0,1),clip="off")+labs(title="Claim boundary")+theme_void(base_family=font_family)+
  theme(plot.title=element_text(face="bold",size=8.2,margin=margin(b=1,l=13)),plot.margin=margin(2,4,3,4))

design<-"\nAAAAABBBBBBB\nAAAAABBBBBBB\nCCCCCCDDDDDD\nEEEEEEEEEEEE\n"
fig3<-tag_panel(pA,"A")+tag_panel(pB,"B")+tag_panel(pC,"C")+tag_panel(pD,"D")+tag_panel(pE,"E")+
  plot_layout(design=design,heights=c(1.45,1.45,1.02,.72))+
  plot_annotation(title="Two-axis evidence model for MAGMA-prioritized genes",
    theme=theme(plot.title=element_text(family=font_family,face="bold",size=12.2,colour=pal[["dark"]],margin=margin(b=6))))

width_in<-7.205;height_in<-6.15;stem<-file.path(fig_dir,"figure3_candidate_evidence_model_draft_v0.4")
svglite::svglite(paste0(stem,".svg"),width=width_in,height=height_in,bg="white",system_fonts=list(sans=font_family));print(fig3);dev.off()
grDevices::pdf(paste0(stem,".pdf"),width=width_in,height=height_in,family=font_family,bg="white",useDingbats=FALSE);print(fig3);dev.off()
ragg::agg_png(paste0(stem,".png"),width=width_in,height=height_in,units="in",res=600,background="white");print(fig3);dev.off()
preview_path<-file.path(log_dir,"figure3_candidate_evidence_model_v0.4_qc_50pct.png")
if(requireNamespace("magick",quietly=TRUE)){im<-magick::image_read(paste0(stem,".png"));im<-magick::image_resize(im,"50%");magick::image_write(im,preview_path,format="png")
}else if(requireNamespace("png",quietly=TRUE)){im<-png::readPNG(paste0(stem,".png"));h<-floor(dim(im)[1]/2);w<-floor(dim(im)[2]/2);ro<-seq.int(1,2*h,2);re<-ro+1;co<-seq.int(1,2*w,2);ce<-co+1;half<-(im[ro,co,,drop=FALSE]+im[re,co,,drop=FALSE]+im[ro,ce,,drop=FALSE]+im[re,ce,,drop=FALSE])/4;png::writePNG(half,preview_path)}
log_lines<-c(paste0("timestamp=",format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z")),"backend=R ggplot2/patchwork/svglite/pdf/ragg","figure=Figure 3 v0.4",
  "claim=two-axis evidence-stratified reporting model; not causal-gene identification or validation","reporting_group_counts=79,1,13,33,4,0",
  "twas_tested=5989","twas_fdr=51","one_snp=42","multi_snp=9","curated_exemplars=UMOD,CASR,CLDN14,CLDN10,HIBADH,PKD2",
  "source_data_changed_from_v0.3=no","claim_changed_from_v0.3=no",paste0("canvas_inches=",width_in,"x",height_in),"font=Helvetica","png_dpi=600","figure2_modified=no","figure4_started=no")
writeLines(log_lines,file.path(log_dir,"stage7E_generate_figure3_v0.4.log"));cat(paste(log_lines,collapse="\n"),"\n")
