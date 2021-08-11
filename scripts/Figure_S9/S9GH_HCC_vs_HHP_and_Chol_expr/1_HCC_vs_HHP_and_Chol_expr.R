#rm(list = ls())
#options(stringsAsFactors = F)

library(dplyr)
library(Seurat)
#library(reticulate)
library(ggplot2)
#library(cowplot)
library(gprofiler2)
#library(gProfileR)
#library(RCurl)
library(RColorBrewer)
#library(tidyverse)
library(tidyr)
library(reshape2)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_S9/S9GH_HCC_vs_HHP_and_Chol_expr"
outs_subpath <- paste0(dir,"/outs/",figure)
data_subpath <- paste0(dir,"/data/",figure)
dir.create(outs_subpath, recursive=TRUE)

other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
other_outs_subpath <- paste0(dir,"/outs/",other_figure)

day_type.integrated <- readRDS(paste0(other_outs_subpath,"/6_samples_integrated_f500_npcs50_res2.rds"))

Idents(day_type.integrated) <- "seurat_clusters"

#RenameIdents
day_type.integrated$Cluster_type<-paste(day_type.integrated$type, day_type.integrated$seurat_clusters,  sep = "_")
Idents(day_type.integrated) <- "Cluster_type"

newcluster.ids <- c("5", "2", "1", "3", "1", "1", "7", "11", "9", "6", "1", "11", "1", "4", "4", "1", "2", "1", "3", "1", "2", "2", "1", "10", "12", "1", "2", "2", "13", "2", "13", "6", "5", "1", "1", "2", "1", "2", "1", "1", "1", "2", "2", "2", "2", "1", "2", "1", "8", "2", "1", "11", "1", "1", "1", "9", "1", "2", "1", "1", "10", "6", "6", "6", "6", "6", "6", "12", "6")
names(newcluster.ids) <- levels(day_type.integrated)
day_type.integrated <- RenameIdents(day_type.integrated, newcluster.ids)
day_type.integrated$merged_clusters <- Idents(day_type.integrated)

lineage.ids <- c("MP", "HCC", "HCC", "Endo", "Endo", "Chol", "HHP", "Endo", "MP", "HSC", "B cell", "NP", "HCC")
names(lineage.ids) <- levels(day_type.integrated)
day_type.integrated <- RenameIdents(day_type.integrated, lineage.ids)
day_type.integrated$lineage_clusters <- Idents(day_type.integrated)

DR_Chol <- c("Spp1", "Sox9", "Sox4", "Cd24a", "Alb", "Afp", "Scd1", "Gpc3", "Krt7", "Krt19", "Epcam", "Mki67")

Idents(day_type.integrated) <- "lineage_clusters"
day_type.integrated_epithelial <- subset(day_type.integrated, idents= c("HCC", "Chol", "HHP"))
day_type.integrated_epithelial$lineage_clusters <- factor(day_type.integrated_epithelial$lineage_clusters, levels = c("HCC", "HHP", "Chol"))
day_type.integrated_epithelial$day_type <- factor(day_type.integrated_epithelial$day_type, levels = c("Day3_Pos", "Day3_Neg", "Day10_Pos", "Day10_Neg", "Day30_Pos", "Day30_Neg"))

DefaultAssay(day_type.integrated_epithelial) <- "RNA"

pdf(paste0(outs_subpath,"/S9H_DR_Chol_no_dots.pdf"), width=15, height=8)
print(VlnPlot(day_type.integrated_epithelial, features = DR_Chol, group.by = "lineage_clusters", cols = c("#8DA0CB", "#FFD92F", "#A6D854"), pt.size = 0))
dev.off()

#Find markers
epithelial.markers <- FindAllMarkers(day_type.integrated_epithelial, only.pos = TRUE, min.pct = 0.25)
write.csv(epithelial.markers, file=paste0(data_subpath,"/HCC_Chol_HHP_markers_ALL.csv"), sep = "\t")

##################################################################################################################################
#GO annotation
#Submit the above three DEG list to gProfiler
#Get a multiquery "gProfiler_mmusculus_4-8-2021_11-25-40 PM__multiquery_ALL.txt"
#Get a multiquery "gProfiler_mmusculus_4-8-2021_11-25-40 PM__multiquery_POS.txt"
#Highlight signature annotations
#export to "ALL_function.txt"
#################################function
data<- read.table(paste0(data_subpath,"/ALL_function.txt"), header=TRUE, sep="\t")
data<- as.data.frame(data)
data$Annotation<- paste(data$Source, data$Annotation, sep = "_")

#Extract data to plot

hm<-data[, c("Annotation", "HCC_ALL","HHP_ALL","CHOL_ALL")]
# hm<-data[, c("Annotation", "HCC_POS","HHP_POS","CHOL_POS")]
hm$Annotation<-factor(hm$Annotation,levels = rev(hm$Annotation))
hm_final<-melt(hm,id=c("Annotation"))
colnames(hm_final)[3]<-c("pvalue")
hm_final$pvalue <- -log10(hm_final$pvalue)

y<- hm_final$Annotation
x<- hm_final$variable

pdf(paste0(outs_subpath,"/S9G_HCC_HHP_CHOL_function_ALL.pdf"), width = 10.2, height = 15)
ggplot(hm_final,aes(x=x, y=y)) +
  geom_tile(aes(fill=pvalue))+
  scale_fill_stepsn(colors = c("#8181C3", "#FFFFFF", "#F0DFDF", "#E2BFBF", "#D39F9F", "#C47F7F", "#B65F5F", "#A83F3F", "#991F1F", "#8B0000"),
                    limits=c(0,16),
                    labels = c("0.301", "1.301", "2", "4", "6", "8", "10", "12", "14"), 
                    breaks = c(0.301, 1.301, 2, 4, 6, 8, 10, 12, 14)
  )+
  theme(
    axis.title.x = element_blank(),
    axis.title.y =  element_blank(),
    axis.text.y = element_text(color="black",size=20),
    axis.text.x = element_text(color="black",size=15, angle=30, margin = margin(t = 10)),
    legend.title = element_text(size = 8,face = "bold"),
    legend.text = element_text(color="black", size=6),
    legend.key = element_rect(fill = "NA")
  )+
  geom_vline(xintercept = seq(-0.5,4.5,1),color="grey", size=2)+
  geom_hline(yintercept = seq(-0.5,23.5,1),color="grey", size=2)
dev.off()
