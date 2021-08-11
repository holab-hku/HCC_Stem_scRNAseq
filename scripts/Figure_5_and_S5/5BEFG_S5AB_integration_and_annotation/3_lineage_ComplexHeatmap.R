library(dplyr)
library(Seurat)
library(reticulate)
library(ggplot2)
library(cowplot)
library(gprofiler2)
#library(gProfileR)
library(RCurl)
library(RColorBrewer)

library(readr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(gplots)
library(whisker)
library(amap)
library(dendroextras)
library(scales)

#TD
#setwd("D:/Lena/scRNAseq/6 samples integration/standard workflow f500")
#setwd("D:/Lena/scRNAseq/6 samples integration/standard workflow f500/rename Fig1")

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
data_subpath <- paste0(dir,"/data/",figure)
outs_subpath <- paste0(dir,"/outs/",figure)

integration_dir <- paste0(outs_subpath,"/integration")

day_type.integrated <- readRDS(file = paste0(outs_subpath,"/6_samples_integrated_f500_npcs50_res2_renamed.rds"))

##################################################################################################################################
DefaultAssay(day_type.integrated) <- "RNA"

day_type.integrated$lineage_clusters <- factor(day_type.integrated$lineage_clusters, levels = c("HCC", "Chol", "HHP",  "Endo", "MP", "HSC", "B cell", "NP"))
day_type.integrated$day <- factor(day_type.integrated$day, levels = c("Day3", "Day10", "Day30"))
day_type.integrated$merged_clusters <- factor(day_type.integrated$merged_clusters, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"))

#Find markers
#day_type.integrated.lineage.markers <- FindAllMarkers(day_type.integrated, only.pos = TRUE, min.pct = 0.25)
#write.table(day_type.integrated.lineage.markers, file = paste0(integration_dir,"/All_clusters_marker_genes_pc50_res2_merged_lineage_RNA.txt"), sep = "\t")
#due to differences between softwares and dependencies, it's recommended to use our generated one using Seurat 3 in 2019
day_type.integrated.lineage.markers <- read.table(paste0(integration_dir,"/All_clusters_marker_genes_pc50_res2_merged_lineage_RNA.txt"), header=T, sep = "\t")


Lineage_heatmap <- filter(day_type.integrated.lineage.markers, p_val_adj < 0.05 
                          & (pct.1-pct.2) > 0.24 
                          & duplicated (day_type.integrated.lineage.markers$gene) == 0)

Lineage_heatmap$cluster <- factor(Lineage_heatmap$cluster, levels = c("HCC", "Chol", "HHP",  "Endo", "MP", "HSC", "B cell", "NP"))
Lineage_heatmap <- Lineage_heatmap[order(Lineage_heatmap[, "cluster"]),]

Lineage_heatmap_test<-mutate(Lineage_heatmap,name=factor(Lineage_heatmap$gene, levels=Lineage_heatmap$gene))
#Lineage_heatmap <- as.matrix(day_type.integrated.lineage.markers, rownames=FALSE)[, c("gene", "cluster")]
#Lineage_heatmap.list <- as.list(Lineage_heatmap[,"gene"])

hm <- as.data.frame(day_type.integrated$RNA@data[row.names(day_type.integrated) %in% Lineage_heatmap$gene, ])
#gene <- as.data.frame(Lineage_heatmap, row.names = Lineage_heatmap$gene) [, c("gene", "cluster")]
cell <- as.data.frame(day_type.integrated@meta.data)[, c("lineage_clusters", "merged_clusters", "day", "type")]
hm_test<-hm[match(Lineage_heatmap_test$name,row.names(hm)),]

merged_hm<-merge.data.frame(cell,t(hm_test),by.x=0,by.y=0)
reordered_hm <- merged_hm[order(merged_hm[, "lineage_clusters"],merged_hm[, "merged_clusters"]),]
plot_hm<-reordered_hm[,c(1,6:1397)]
row.names(plot_hm)<-plot_hm[,c(1)]
plot_hm<-plot_hm[,-c(1)]

ha = HeatmapAnnotation(Type=reordered_hm$type,
                       Day=reordered_hm$day,
                       Lineage=reordered_hm$lineage_clusters,
                       Cluster=reordered_hm$merged_clusters,
                       col=list(Type=c("Neg"="#bebebe","Pos"="#fb6a4a"),
                                Lineage=c("MP"="#FC8D62", "HCC"="#8DA0CB", "Endo"="#E78AC3", "Chol"="#A6D854", "HHP"="#FFD92F", "HSC"="#E5C494", "B cell"="#66C2A5", "NP"="#B3B3B3"),
                                Day=c("Day3"="#bdd7e7", "Day10"="#6baed6", "Day30"="#2171b5"),
                                Cluster=c("5"="#F8766D", "2"="#E18A00", "1"= "#BE9C00", "3"=  "#8CAB00", "7"=  "#24B700", "11"=  "#00BE70", "9"=  "#00C1AB", "6"=  "#00BBDA", "4"=  "#00ACFC", "10"=  "#8B93FF", "12"=  "#D575FE", "13"=  "#F962DD", "8"=  "#FF65AC"))
                       #,simple_anno_size = unit(2, "cm")
)

cell$cluster <- factor(cell$lineage_clusters, levels = c("HCC", "Chol", "HHP",  "Endo", "MP", "HSC", "B cell", "NP"))
cell <- cell[order(cell[, "lineage_clusters"]),]

#tiff('heatmap1.tiff', units="px", width=4800, height=3200, res=600) #save file for high resolution
pdf(paste0(outs_subpath,"/5F_Lineage_heatmap.pdf"), width=12, height=18)
Heatmap(t(scale(plot_hm)), name = "Z-score",
        width = unit(21, "cm"), height = unit(30, "cm"),
        col = colorRamp2(c(-2, 0,2), c("#030389", "white", "#8b0000")),
        show_row_dend = FALSE,
        cluster_rows=FALSE,cluster_columns = FALSE,
        clustering_distance_columns="euclidean",clustering_method_columns = "ward.D2",
        clustering_distance_rows="euclidean",clustering_method_rows = "ward.D2",
        column_dend_height = unit(1, "cm"),
        row_dend_width = unit(1, "cm"),
        show_row_names = FALSE, show_column_names = FALSE, row_names_gp = gpar(fontsize = 12),column_names_gp = gpar(fontsize = 12),
        row_split = factor(Lineage_heatmap_test$cluster),
        column_split = factor(cell$cluster),
        top_annotation = ha
) 
dev.off()


