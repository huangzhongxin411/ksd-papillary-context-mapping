#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(svglite); library(ragg)})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
fig_dir <- file.path(root, "results/figures/revision/stage7E_figure3_final_micro_polish_v0.6R2")
table_dir <- file.path(root, "results/tables/revision/stage7E_figure3_final_micro_polish_v0.6R2")
log_dir <- file.path(root, "logs/revision/stage7E_figure3_final_micro_polish_v0.6R2")
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

# A: manually drawn three-layer ribbon flow.
ribbon <- function(id,x0,x1,l0,l1,r0,r1,fill){
  x<-seq(x0,x1,length.out=60);t<-(x-x0)/(x1-x0);s<-t*t*(3-2*t)
  top<-l1+(r1-l1)*s;bottom<-l0+(r0-l0)*s
  data.table(id=id,x=c(x,rev(x)),y=c(top,rev(bottom)),fill=fill)
}
ribbons <- rbindlist(list(
  ribbon("L-M1",.18,.43,.080,.681,.085,.643,"#A9C1C6"),
  ribbon("L-M2",.18,.43,.681,.894,.668,.866,"#DDE3E4"),
  ribbon("L-M3",.18,.43,.894,.920,.891,.915,"#E6D8AF"),
  ribbon("M1-R1",.53,.80,.085,.559,.104,.542,"#7F9DA6"),
  ribbon("M1-R2",.53,.80,.559,.565,.560,.566,"#B99B5A"),
  ribbon("M1-R3",.53,.80,.565,.643,.584,.656,"#D9A99C"),
  ribbon("M2-R4",.53,.80,.668,.866,.674,.857,"#DDE3E4"),
  ribbon("M3-R5",.53,.80,.891,.915,.875,.897,"#E6D8AF")
))
left_nodes<-data.table(xmin=.08,xmax=.18,ymin=.08,ymax=.92,fill=pal[["light"]])
mid_nodes<-data.table(xmin=.43,xmax=.53,ymin=c(.085,.668,.891),ymax=c(.643,.866,.915),fill=c(pal[["teal_pale"]],pal[["light"]],pal[["gold_pale"]]))
right_nodes<-data.table(xmin=.80,xmax=.90,ymin=c(.104,.560,.584,.674,.875),ymax=c(.542,.566,.656,.857,.897),fill=c(pal[["teal_mid"]],pal[["gold"]],"#D9A99C",pal[["line"]],pal[["gold_pale"]]))
pA<-ggplot()+
  geom_polygon(data=ribbons,aes(x,y,group=id,fill=fill),alpha=.78,colour=NA)+scale_fill_identity()+
  geom_rect(data=left_nodes,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill),colour=pal[["line"]],linewidth=.35)+
  geom_rect(data=mid_nodes,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill),colour=pal[["white"]],linewidth=.4)+
  geom_rect(data=right_nodes,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill),colour=pal[["white"]],linewidth=.35)+
  annotate("text",x=.13,y=.54,label="130",family=font_family,fontface="bold",size=3.6,colour=pal[["teal"]])+
  annotate("text",x=.13,y=.44,label="candidate genes",family=font_family,size=1.55,colour=pal[["grey"]])+
  annotate("text",x=.48,y=.36,label="MAGMA Bonferroni\n93",family=font_family,fontface="bold",size=1.92,lineheight=.88)+
  annotate("text",x=.48,y=.767,label="Lower MAGMA\n33",family=font_family,fontface="bold",size=1.76,lineheight=.88)+
  annotate("text",x=.48,y=.955,label="TWAS proxy only  4",family=font_family,fontface="bold",size=1.58,colour=pal[["gold_dark"]])+
  annotate("text",x=.85,y=.323,label="R1  MAGMA only\n79",family=font_family,fontface="bold",size=1.9,lineheight=.88)+
  annotate("text",x=.85,y=.620,label="R3  one-SNP\n13",family=font_family,fontface="bold",size=1.55,lineheight=.88)+
  annotate("text",x=.85,y=.765,label="R4  lower MAGMA\n33",family=font_family,fontface="bold",size=1.6,lineheight=.88)+
  annotate("segment",x=.90,xend=.938,y=.563,yend=.563,colour=pal[["gold"]],linewidth=.50)+
  annotate("text",x=.943,y=.563,label="R2  multi-SNP  1",hjust=0,family=font_family,fontface="bold",size=1.52,colour=pal[["gold_dark"]])+
  annotate("segment",x=.90,xend=.938,y=.886,yend=.886,colour=pal[["gold"]],linewidth=.50)+
  annotate("text",x=.943,y=.886,label="R5  proxy only  4",hjust=0,family=font_family,fontface="bold",size=1.52,colour=pal[["gold_dark"]])+
  annotate("segment",x=.90,xend=.90,y=.935,yend=.978,colour=pal[["mid"]],linewidth=.68)+
  annotate("text",x=.92,y=.978,label="R6  contextual / unresolved  0",hjust=0,family=font_family,fontface="bold",size=1.43,colour=pal[["grey"]])+
  annotate("text",x=.13,y=.995,label="CANDIDATES",family=font_family,fontface="bold",size=1.68,colour=pal[["grey"]])+
  annotate("text",x=.48,y=.995,label="EVIDENCE AXIS",family=font_family,fontface="bold",size=1.68,colour=pal[["grey"]])+
  annotate("text",x=.85,y=.995,label="REPORTING GROUP",family=font_family,fontface="bold",size=1.68,colour=pal[["grey"]])+
  coord_cartesian(xlim=c(.04,1.08),ylim=c(.04,1.03),clip="off")+
  labs(title="Candidate gene evidence-flow assignment",subtitle="Mutually exclusive reporting groups; not causal tiers")+theme_void(base_family=font_family)+
  theme(plot.title=element_text(face="bold",size=8.6,margin=margin(b=1,l=13)),plot.subtitle=element_text(size=6.6,colour=pal[["grey"]],margin=margin(b=2,l=13)),plot.margin=margin(3,7,3,5))

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
  annotate("text",x=17.2,y=4.42,label="1 tile = 1 FDR-supported TWAS model",hjust=1,family=font_family,size=1.45,colour=pal[["grey"]])+
  annotate("text",x=1,y=.20,label="42  one-SNP model",hjust=0,family=font_family,fontface="bold",size=1.75,colour=pal[["coral"]])+
  annotate("text",x=8,y=.20,label="9  multi-SNP model",hjust=0,family=font_family,fontface="bold",size=1.75,colour=pal[["gold_dark"]])+
  annotate("text",x=1,y=-.45,label="GTEx kidney cortex proxy  ·  supplementary only  ·  no papilla-specific eQTL",hjust=0,family=font_family,size=1.55,colour=pal[["grey"]])+
  coord_cartesian(xlim=c(.2,17.8),ylim=c(-.65,4.55),clip="off")+
  labs(title="TWAS proxy burden")+theme_void(base_family=font_family)+
  theme(plot.title=element_text(face="bold",size=8.2,margin=margin(b=2,l=13)),plot.margin=margin(4,7,4,5))

# C: curated exemplars as interpretive flags only.
exemplar_names <- c("UMOD","CASR","CLDN14","CLDN10","HIBADH","PKD2")
if (!setequal(exemplars$gene,exemplar_names)) stop("Exemplar set differs from lock")
chip_x<-c(.17,.50,.83,.17,.50,.83); chip_y<-c(.54,.54,.54,.30,.30,.30)
pC <- ggplot()+
  annotate("rect",xmin=.13,xmax=.87,ymin=.76,ymax=.91,fill=pal[["gold_pale"]],colour=pal[["gold"]],linewidth=.32)+
  annotate("text",x=.5,y=.835,label="Interpretive exemplars only",family=font_family,fontface="bold",size=1.9,colour=pal[["gold_dark"]])+
  annotate("text",x=.5,y=.68,label="Biology-informed  ·  no evidence upgrade",family=font_family,size=1.6,colour=pal[["grey"]])
for(i in seq_along(exemplar_names)) pC<-pC+annotate("rect",xmin=chip_x[i]-.098,xmax=chip_x[i]+.098,ymin=chip_y[i]-.065,ymax=chip_y[i]+.065,fill="#FBFAF6",colour="#D2C39E",linewidth=.23)+
  annotate("text",x=chip_x[i],y=chip_y[i],label=exemplar_names[i],family=font_family,fontface="bold.italic",size=1.68)
pC <- pC+annotate("text",x=.5,y=.09,label="not validated  ·  not SMR/coloc-supported",family=font_family,size=1.5,colour=pal[["grey"]])+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1),clip="off")+labs(title="Interpretive exemplar inset")+theme_void(base_family=font_family)+
  theme(plot.title=element_text(face="bold",size=8.1,margin=margin(b=2,l=13)),plot.margin=margin(4,5,4,5))

# D: compact paired claim bands.
pD<-ggplot()+
  annotate("rect",xmin=.02,xmax=.61,ymin=.22,ymax=.68,fill=pal[["teal_pale"]],colour=pal[["teal2"]],linewidth=.38)+
  annotate("rect",xmin=.64,xmax=.985,ymin=.22,ymax=.68,fill=pal[["coral_pale"]],colour=pal[["coral"]],linewidth=.38)+
  annotate("text",x=.035,y=.80,label="SUPPORTED REPORTING",hjust=0,family=font_family,fontface="bold",size=1.9,colour=pal[["teal"]])+
  annotate("text",x=.655,y=.80,label="NOT CLAIMED",hjust=0,family=font_family,fontface="bold",size=1.9,colour=pal[["coral"]])+
  annotate("text",x=.315,y=.45,label="MAGMA priority  ·  kidney cortex proxy  ·  interpretive exemplars",family=font_family,fontface="bold",size=1.72)+
  annotate("text",x=.8125,y=.45,label="causal / validated genes\nSMR/coloc support  ·  therapeutic targets",family=font_family,fontface="bold",size=1.62,lineheight=.9)+
  coord_cartesian(xlim=c(0,1.025),ylim=c(0,1),clip="off")+labs(title="Reporting boundary")+theme_void(base_family=font_family)+
  theme(plot.title=element_text(face="bold",size=8.2,margin=margin(b=1,l=13)),plot.margin=margin(2,4,3,4))

design<-"\nAAAAAAAAAAAA\nAAAAAAAAAAAA\nBBBBBBBCCCCC\nDDDDDDDDDDDD\n"
fig3<-tag_panel(pA,"A")+tag_panel(pB,"B")+tag_panel(pC,"C")+tag_panel(pD,"D")+
  plot_layout(design=design,heights=c(.92,.92,1.30,.58))+
  plot_annotation(title="Evidence-stratified candidate gene reporting model",subtitle="MAGMA genetic priority is separated from kidney cortex proxy TWAS support",
    theme=theme(plot.title=element_text(family=font_family,face="bold",size=12.2,colour=pal[["dark"]],margin=margin(b=2)),plot.subtitle=element_text(family=font_family,size=8,colour=pal[["grey"]],margin=margin(b=6))))

width_in<-7.205;height_in<-5.45;stem<-file.path(fig_dir,"figure3_candidate_evidence_model_draft_v0.6R2")
svglite::svglite(paste0(stem,".svg"),width=width_in,height=height_in,bg="white",system_fonts=list(sans=font_family));print(fig3);dev.off()
grDevices::pdf(paste0(stem,".pdf"),width=width_in,height=height_in,family=font_family,bg="white",useDingbats=FALSE);print(fig3);dev.off()
ragg::agg_png(paste0(stem,".png"),width=width_in,height=height_in,units="in",res=600,background="white");print(fig3);dev.off()
preview_path<-file.path(log_dir,"figure3_candidate_evidence_model_v0.6R2_qc_50pct.png")
if(requireNamespace("magick",quietly=TRUE)){im<-magick::image_read(paste0(stem,".png"));im<-magick::image_resize(im,"50%");magick::image_write(im,preview_path,format="png")
}else if(requireNamespace("png",quietly=TRUE)){im<-png::readPNG(paste0(stem,".png"));h<-floor(dim(im)[1]/2);w<-floor(dim(im)[2]/2);ro<-seq.int(1,2*h,2);re<-ro+1;co<-seq.int(1,2*w,2);ce<-co+1;half<-(im[ro,co,,drop=FALSE]+im[re,co,,drop=FALSE]+im[ro,ce,,drop=FALSE]+im[re,ce,,drop=FALSE])/4;png::writePNG(half,preview_path)}
log_lines<-c(paste0("timestamp=",format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z")),"backend=R ggplot2/patchwork/svglite/pdf/ragg","figure=Figure 3 v0.6R2",
  "claim=two-axis evidence-stratified reporting model; not causal-gene identification or validation","reporting_group_counts=79,1,13,33,4,0",
  "twas_tested=5989","twas_fdr=51","one_snp=42","multi_snp=9","curated_exemplars=UMOD,CASR,CLDN14,CLDN10,HIBADH,PKD2",
  "baseline=Figure 3 v0.6R production proof","structure_changed_from_v0.6R=no","source_data_changed_from_v0.6R=no","claim_changed_from_v0.6R=no",paste0("canvas_inches=",width_in,"x",height_in),"font=Helvetica","png_dpi=600","figure2_modified=no","figure4_started=no")
writeLines(log_lines,file.path(log_dir,"stage7E_generate_figure3_v0.6R2.log"));cat(paste(log_lines,collapse="\n"),"\n")
