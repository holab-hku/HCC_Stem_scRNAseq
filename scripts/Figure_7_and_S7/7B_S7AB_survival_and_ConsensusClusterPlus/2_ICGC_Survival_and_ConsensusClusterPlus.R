rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)
library(glmnet)
library(ROCR)
library(caret)
library(ggpubr)
library(patchwork)
library(survival)
library(survminer)
library(ggthemes)
library(Seurat)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ConsensusClusterPlus)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_7_and_S7/7B_S7AB_survival_and_ConsensusClusterPlus"
outs_subpath <- paste0(dir,"/outs/",figure)
data_subpath <- paste0(dir,"/data/",figure)
data_subpath_0 <- paste0(dir,"/data/Figure_7_and_S7")
dir.create(outs_subpath, recursive=TRUE)

load(paste0(data_subpath,"/ICGC-JP-sur_model.Rdata"))
#load(file = "signature.RData")

meta$time=meta$time/30.5
row.names(meta)=meta$ID
exprSet=log2(edgeR::cpm(exprSet)+1)

mouse_to_human_gene_names <- function(mouse_genes_list) {
  human_to_mouse_gene <- read.table(paste0(data_subpath_0,"/human_to_mouse_gene.txt"), header=T, sep="\t")
  rownames(human_to_mouse_gene) <- human_to_mouse_gene$mouse_gene
  
  intersected_genes <- intersect(mouse_genes_list, rownames(human_to_mouse_gene))
  
  human_to_mouse_gene_intersect <- human_to_mouse_gene[intersected_genes,]
  final_list <- human_to_mouse_gene_intersect$human_gene
  
  no_intersect <- setdiff(mouse_genes_list, intersected_genes)
  if (length(no_intersect) > 0) {
    cat("The following genes can't be found in the MGI table, will convert them anyway from the lowercase format:\n")
    print(no_intersect)
    
    no_intersect_lc <- str_to_upper(no_intersect)
    print(no_intersect_lc)
    
    final_list <- c(final_list, no_intersect_lc)
  }
  return(final_list)
  
}

#Pos_DEGs (ca 800)
Posgenes <- read.table(paste0(data_subpath,"/Table S7_Conserved markers of Prom1+ lineage.txt"), header = FALSE)$V1
genes <- mouse_to_human_gene_names(Posgenes)
genes <- genes[duplicated(genes)==FALSE]

#Select patient
metatime <- filter(meta, treatment=="no treatment")

#ConsensusClusterPlus with selected genes
exprSet_genes=exprSet[(row.names(exprSet) %in% genes) == TRUE, str_split(colnames(exprSet), "_", simplify = T)[,1] %in% row.names(metatime)]
exprSet_genes = sweep(exprSet_genes, 1, apply(exprSet_genes,1,median,na.rm=T))

ICGC_outs_path <- paste0(outs_subpath,"/S7A_consensus_matrix/ICGC")
dir.create(ICGC_outs_path, recursive=TRUE)
title=ICGC_outs_path

title=ICGC_outs_path

results = ConsensusClusterPlus(as.matrix(exprSet_genes),maxK=4,reps=50,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="km",distance="pearson",seed=111,plot="png")

icl = calcICL(results,title=title,plot="png")

results[[2]][["consensusClass"]][1:10]

thisPal <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", 
             "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", 
             "#6A3D9A", "#FFFF99", "#B15928", "#bd18ea", "#2ef4ca", 
             "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106", 
             "#85848f", "#000000", "#076f25", "#93cd7f", "#4d0776", 
             "#ffffff")


survival_outs_dir <- paste0(outs_subpath,"/7B_S7B_survival/ICGC")
dir.create(survival_outs_dir, recursive=TRUE)
#for (i in 2:4){
for (i in 2:2){
  y=metatime
  y$group <- results[[i]][["consensusClass"]]
  
  Groupcol = thisPal[1:i]
  names(Groupcol)=as.character(1:i)
  
  #KMplot
  sfit1=survfit(Surv(time, event)~group, data=y)
  ggsurvplot(sfit1,
             pval =TRUE,
             pval.size = 5,
             data = y, 
             risk.table = F,
             palette = Groupcol,
             xlim = c(0,80),
             break.time.by = 20,
             xlab ="Time (months)",
             title = "ICGC-JP",
             legend = c(0.78,0.8),
             legend.title = "Signature score",
             legend.labs = as.character(c(1:i)),
             font.legend = 12,
             ggtheme = theme_base())
  ggsave(paste0(survival_outs_dir,"/ICGC_signature_KMplot_ConsensusClusterPlus", i, ".pdf"), width = 4, height = 4)
  print(pairwise_survdiff(Surv(time, event)~group, data=y, p.adjust.method = "BH", rho = 0))
  
  #gene_patient heatmap
  x=exprSet[(row.names(exprSet) %in% genes) == TRUE, str_split(colnames(exprSet), "_", simplify = T)[,1] %in% row.names(y)]
  y=y[match(str_split(colnames(x), "_", simplify = T)[,1],y$ID),]
  all(y$ID==str_split(colnames(x), "_", simplify = T)[,1])
  
  test=t(scale(t(x)))
  Groupcol = thisPal[1:i]
  names(Groupcol)=as.character(1:i)
  ha = HeatmapAnnotation(Group=y$group,
                         col = list(Group=Groupcol))
  
  pdf(paste0(survival_outs_dir,"/ICGC_signature_KMplot_ConsensusClusterPlus", i, "_hm.pdf"), width = 12, height = 11)
  print(Heatmap(test, 
                col = colorRamp2(c(-2, 0,2), c("#030389", "white", "#8b0000")), 
                top_annotation = ha,
                cluster_columns = F,
                clustering_distance_rows = "euclidean",
                clustering_method_rows = "ward.D",
                column_split = factor(y$group),
                show_column_names = FALSE,
                show_row_names = FALSE,
                show_row_dend = FALSE))
  dev.off()
  
  #boxplot
  pl=mutate(y, score = colMeans(test))
  pl$group <- as.factor(pl$group)
  pdf(paste0(survival_outs_dir,"/ICGC_signature_KMplot_ConsensusClusterPlus", i, "_box.pdf"), width = 3.5, height = 5)
  print(ggplot(pl, aes(x=group, y=score)) + 
          geom_boxplot(lwd=1) + 
          geom_jitter(aes(color=group),shape=16, position=position_jitter(0.35)) +
          scale_color_manual(values=thisPal[1:i]) +
          theme(legend.position = 'none', 
                axis.title.x = element_blank(), 
                axis.title.y = element_text(size=20), 
                axis.text.x = element_text(angle = 0, hjust=0.5, size = 18, colour = "black"),
                axis.text.y = element_text(size = 18, colour = "black"),
                plot.background = element_rect(color = "white"),
                panel.background = element_rect(fill = NA),
                plot.margin = unit(c(10, 10, 10, 10), "pt"),
                axis.line.x = element_line(colour = "black"),
                axis.line.y = element_line(colour = "black")) + 
          ylab("Expression of signature (z-scored)") +
          stat_compare_means(method = "t.test", label.y= -1, size = 5))
  dev.off()
}
