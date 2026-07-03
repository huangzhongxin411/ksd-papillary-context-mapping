#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(svglite); library(ragg)})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
fig_dir <- file.path(root, "results/figures/revision/stage7E_figure3_graphical_flow_v0.5")
table_dir <- file.path(root, "results/tables/revision/stage7E_figure3_graphical_flow_v0.5")
log_dir <- file.path(root, "logs/revision/stage7E_figure3_graphical_flow_v0.5")
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

# A: proportional evidence-flow mosaic. Width equals gene count.
flow <- data.table(
  group=c("R1","R4","R2","R3","R5"), n=c(79,33,1,13,4),
  label=c("R1\nMAGMA only\n79","R4\nlower MAGMA\n33","R2\nmulti-SNP\n1","R3\none-SNP\n13","R5\nproxy only\n4"),
  fill=c(pal[["teal_mid"]],pal[["line"]],pal[["gold"]],"#D9A99C",pal[["gold_pale"]])
)
flow[,xmax:=cumsum(n)];flow[,xmin:=shift(xmax,fill=0)];flow[,xmid:=(xmin+xmax)/2]
pA <- ggplot()+
  annotate("rect",xmin=-27,xmax=-10,ymin=.31,ymax=.69,fill=pal[["light"]],colour=pal[["line"]],linewidth=.4)+
  annotate("text",x=-18.5,y=.54,label="130",family=font_family,fontface="bold",size=4.0,colour=pal[["teal"]])+
  annotate("text",x=-18.5,y=.39,label="candidate genes",family=font_family,size=1.75,colour=pal[["grey"]])+
  annotate("segment",x=-8,xend=-1.5,y=.50,yend=.50,colour=pal[["mid"]],linewidth=.7,arrow=grid::arrow(length=grid::unit(4,"pt"),type="closed"))+
  geom_rect(data=flow,aes(xmin=xmin,xmax=xmax,ymin=.25,ymax=.75,fill=fill),colour=pal[["white"]],linewidth=.65)+
  geom_text(data=flow[n>=10],aes(x=xmid,y=.50,label=label),family=font_family,fontface="bold",size=2.0,lineheight=.88)+
  annotate("text",x=112.5,y=.50,label="R2\n1",family=font_family,fontface="bold",size=1.65,lineheight=.85)+
  annotate("text",x=128,y=.50,label="R5\n4",family=font_family,fontface="bold",size=1.65,lineheight=.85)+
  annotate("segment",x=0,xend=112,y=.86,yend=.86,colour=pal[["teal"]],linewidth=.65,lineend="round")+
  annotate("text",x=56,y=.94,label="MAGMA evidence  ·  112 without proxy",family=font_family,fontface="bold",size=1.8,colour=pal[["teal"]])+
  annotate("segment",x=112,xend=126,y=.86,yend=.86,colour=pal[["gold"]],linewidth=.65,lineend="round")+
  annotate("text",x=119,y=.94,label="MAGMA + proxy  14",family=font_family,fontface="bold",size=1.55,colour=pal[["gold_dark"]])+
  annotate("segment",x=130.5,xend=130.5,y=.23,yend=.77,colour=pal[["mid"]],linewidth=.7)+
  annotate("text",x=132,y=.12,label="R6 = 0",hjust=1,family=font_family,fontface="bold",size=1.55,colour=pal[["grey"]])+
  scale_fill_identity()+coord_cartesian(xlim=c(-28,133),ylim=c(0,1.06),clip="off")+
  labs(title="Evidence-flow assignment",subtitle="Mutually exclusive reporting groups; not causal tiers")+theme_void(base_family=font_family)+
  theme(plot.title=element_text(face="bold",size=8.6,margin=margin(b=1,l=13)),plot.subtitle=element_text(size=6.6,colour=pal[["grey"]],margin=margin(b=2,l=13)),plot.margin=margin(3,6,3,5))

# B: 51-tile TWAS model burden.
waffle <- data.table(id=1:51)
waffle[,class:=ifelse(id<=42,"One-SNP","Multi-SNP")]
waffle[,x:=(id-1)%%17+1];waffle[,y:=3-floor((id-1)/17)]
pB <- ggplot(waffle,aes(x,y,fill=class))+
  geom_tile(width=.78,height=.78,colour=pal[["white"]],linewidth=.25)+
  scale_fill_manual(values=c("One-SNP"="#D9A99C","Multi-SNP"=pal[["gold"]]),guide="none")+
  annotate("rect",xmin=.4,xmax=5.7,ymin=3.65,ymax=4.35,fill=pal[["light"]],colour=pal[["line"]],linewidth=.35)+
  annotate("rect",xmin=8.0,xmax=13.3,ymin=3.65,ymax=4.35,fill=pal[["teal_pale"]],colour=pal[["teal_mid"]],linewidth=.35)+
  annotate("segment",x=6.0,xend=7.7,y=4,yend=4,colour=pal[["mid"]],linewidth=.55,arrow=grid::arrow(length=grid::unit(3,"pt"),type="closed"))+
  annotate("text",x=3.05,y=4,label="5,989 tested",family=font_family,fontface="bold",size=2.1)+
  annotate("text",x=10.65,y=4,label="51 FDR-supported",family=font_family,fontface="bold",size=2.0)+
  annotate("text",x=1,y=.20,label="42 one-SNP",hjust=0,family=font_family,fontface="bold",size=1.8,colour=pal[["coral"]])+
  annotate("text",x=8,y=.20,label="9 multi-SNP",hjust=0,family=font_family,fontface="bold",size=1.8,colour=pal[["gold_dark"]])+
  annotate("text",x=1,y=-.45,label="GTEx kidney cortex proxy  ·  supplementary only  ·  no papilla-specific eQTL",hjust=0,family=font_family,size=1.55,colour=pal[["grey"]])+
  coord_cartesian(xlim=c(.2,17.8),ylim=c(-.65,4.55),clip="off")+
  labs(title="TWAS proxy burden")+theme_void(base_family=font_family)+
  theme(plot.title=element_text(face="bold",size=8.2,margin=margin(b=2,l=13)),plot.margin=margin(4,7,4,5))

# C: curated exemplars as interpretive flags only.
exemplar_names <- c("UMOD","CASR","CLDN14","CLDN10","HIBADH","PKD2")
if (!setequal(exemplars$gene,exemplar_names)) stop("Exemplar set differs from lock")
chip_x<-c(.16,.49,.82,.16,.49,.82); chip_y<-c(.64,.64,.64,.34,.34,.34)
pC <- ggplot()+annotate("rect",xmin=.03,xmax=.97,ymin=.18,ymax=.80,fill=pal[["gold_pale"]],colour=pal[["gold"]],linewidth=.35)+
  annotate("text",x=.5,y=.90,label="Biology-informed flag  ·  No evidence upgrade",family=font_family,fontface="bold",size=1.9,colour=pal[["gold_dark"]])
for(i in seq_along(exemplar_names)) pC<-pC+annotate("rect",xmin=chip_x[i]-.125,xmax=chip_x[i]+.125,ymin=chip_y[i]-.10,ymax=chip_y[i]+.10,fill=pal[["white"]],colour="#D2C39E",linewidth=.28)+
  annotate("text",x=chip_x[i],y=chip_y[i],label=exemplar_names[i],family=font_family,fontface="bold.italic",size=2.05)
pC <- pC+annotate("text",x=.5,y=.08,label="Interpretive anchors only",family=font_family,size=1.65,colour=pal[["grey"]])+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1),clip="off")+labs(title="Curated exemplar boundary")+theme_void(base_family=font_family)+
  theme(plot.title=element_text(face="bold",size=8.1,margin=margin(b=2,l=13)),plot.margin=margin(4,5,4,5))

# D: compact paired claim bands.
pD<-ggplot()+
  annotate("rect",xmin=.02,xmax=.61,ymin=.22,ymax=.68,fill=pal[["teal_pale"]],colour=pal[["teal2"]],linewidth=.38)+
  annotate("rect",xmin=.64,xmax=.985,ymin=.22,ymax=.68,fill=pal[["coral_pale"]],colour=pal[["coral"]],linewidth=.38)+
  annotate("text",x=.035,y=.80,label="SUPPORTED REPORTING",hjust=0,family=font_family,fontface="bold",size=1.9,colour=pal[["teal"]])+
  annotate("text",x=.655,y=.80,label="NOT CLAIMED",hjust=0,family=font_family,fontface="bold",size=1.9,colour=pal[["coral"]])+
  annotate("text",x=.315,y=.45,label="Transparent reporting  ·  MAGMA priority  ·  kidney cortex proxy context  ·  interpretive exemplars",family=font_family,fontface="bold",size=1.52)+
  annotate("text",x=.8125,y=.45,label="Causal genes  ·  validated genes\nSMR/coloc support  ·  therapeutic targets",family=font_family,fontface="bold",size=1.48,lineheight=.9)+
  coord_cartesian(xlim=c(0,1.025),ylim=c(0,1),clip="off")+labs(title="Claim boundary")+theme_void(base_family=font_family)+
  theme(plot.title=element_text(face="bold",size=8.2,margin=margin(b=1,l=13)),plot.margin=margin(2,4,3,4))

design<-"\nAAAAAAAAAAAA\nAAAAAAAAAAAA\nBBBBBBBCCCCC\nDDDDDDDDDDDD\n"
fig3<-tag_panel(pA,"A")+tag_panel(pB,"B")+tag_panel(pC,"C")+tag_panel(pD,"D")+
  plot_layout(design=design,heights=c(.92,.92,1.30,.58))+
  plot_annotation(title="Evidence-stratified candidate gene reporting model",subtitle="MAGMA genetic priority is separated from kidney cortex proxy TWAS support",
    theme=theme(plot.title=element_text(family=font_family,face="bold",size=12.2,colour=pal[["dark"]],margin=margin(b=2)),plot.subtitle=element_text(family=font_family,size=8,colour=pal[["grey"]],margin=margin(b=6))))

width_in<-7.205;height_in<-5.45;stem<-file.path(fig_dir,"figure3_candidate_evidence_model_draft_v0.5")
svglite::svglite(paste0(stem,".svg"),width=width_in,height=height_in,bg="white",system_fonts=list(sans=font_family));print(fig3);dev.off()
grDevices::pdf(paste0(stem,".pdf"),width=width_in,height=height_in,family=font_family,bg="white",useDingbats=FALSE);print(fig3);dev.off()
ragg::agg_png(paste0(stem,".png"),width=width_in,height=height_in,units="in",res=600,background="white");print(fig3);dev.off()
preview_path<-file.path(log_dir,"figure3_candidate_evidence_model_v0.5_qc_50pct.png")
if(requireNamespace("magick",quietly=TRUE)){im<-magick::image_read(paste0(stem,".png"));im<-magick::image_resize(im,"50%");magick::image_write(im,preview_path,format="png")
}else if(requireNamespace("png",quietly=TRUE)){im<-png::readPNG(paste0(stem,".png"));h<-floor(dim(im)[1]/2);w<-floor(dim(im)[2]/2);ro<-seq.int(1,2*h,2);re<-ro+1;co<-seq.int(1,2*w,2);ce<-co+1;half<-(im[ro,co,,drop=FALSE]+im[re,co,,drop=FALSE]+im[ro,ce,,drop=FALSE]+im[re,ce,,drop=FALSE])/4;png::writePNG(half,preview_path)}
log_lines<-c(paste0("timestamp=",format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z")),"backend=R ggplot2/patchwork/svglite/pdf/ragg","figure=Figure 3 v0.5",
  "claim=two-axis evidence-stratified reporting model; not causal-gene identification or validation","reporting_group_counts=79,1,13,33,4,0",
  "twas_tested=5989","twas_fdr=51","one_snp=42","multi_snp=9","curated_exemplars=UMOD,CASR,CLDN14,CLDN10,HIBADH,PKD2",
  "source_data_changed_from_v0.4=no","claim_changed_from_v0.4=no",paste0("canvas_inches=",width_in,"x",height_in),"font=Helvetica","png_dpi=600","figure2_modified=no","figure4_started=no")
writeLines(log_lines,file.path(log_dir,"stage7E_generate_figure3_v0.5.log"));cat(paste(log_lines,collapse="\n"),"\n")
