rm(list = ls())
options(stringsAsFactors = F)

library(dplyr)
library(Seurat)
library(reticulate)
library(ggplot2)
library(cowplot)
library(gprofiler2)
library(RCurl)
library(RColorBrewer)

#TD
#setwd("E:/Seq data/scRNAseq/6 samples integration/standard workflow f500")
dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
data_subpath <- paste0(dir,"/data/",figure)
outs_subpath <- paste0(dir,"/outs/",figure)

day_type.integrated <- readRDS(paste0(outs_subpath,"/6_samples_integrated_f500_npcs50_res2.rds"))

integration_dir <- paste0(outs_subpath,"/integration")

Idents(day_type.integrated) <- "seurat_clusters"

pdf(paste0(outs_subpath,"/S5A_Clustering_tsne_pc50_res2_seurat_clusters.pdf"), width=10.7, height=10)
DimPlot(day_type.integrated, reduction = "tsne", pt.size=1)
dev.off()

#RenameIdents
day_type.integrated$Cluster_type<-paste(day_type.integrated$type, day_type.integrated$seurat_clusters,  sep = "_")
Idents(day_type.integrated) <- "Cluster_type"

newcluster.ids <- c("5", "2", "1", "3", "1", "1", "7", "11", "9", "6", "1", "11", "1", "4", "4", "1", "2", "1", "3", "1", "2", "2", "1", "10", "12", "1", "2", "2", "13", "2", "13", "6", "5", "1", "1", "2", "1", "2", "1", "1", "1", "2", "2", "2", "2", "1", "2", "1", "8", "2", "1", "11", "1", "1", "1", "9", "1", "2", "1", "1", "10", "6", "6", "6", "6", "6", "6", "12", "6")
names(newcluster.ids) <- levels(day_type.integrated)
day_type.integrated <- RenameIdents(day_type.integrated, newcluster.ids)
day_type.integrated$merged_clusters <- Idents(day_type.integrated)

##################################################################################################################################
write.table(table(day_type.integrated$merged_clusters, day_type.integrated$day_type), file = paste0(integration_dir,"/1_to_13_fraction"), sep = "\t")

pdf(paste0(outs_subpath,"/5G_Clustering_tsne_pc50_res2_pos_neg.pdf"), width=19.5, height=10)
DimPlot(day_type.integrated, reduction = "tsne", split.by= "type", group.by = "type", cols = c("grey", "#FB6A4A"), pt.size=2)
dev.off()

pdf(paste0(outs_subpath,"/5B_Clustering_tsne_pc50_res2_merged_1_to_13_label.pdf"), width=10.7, height=10)
DimPlot(day_type.integrated, reduction = "tsne", label=TRUE, pt.size=2)
dev.off()

pdf(paste0(integration_dir,"/Clustering_tsne_pc50_res2_merged_1_to_13.pdf"), width=10.7, height=10)
DimPlot(day_type.integrated, reduction = "tsne", pt.size=2)
dev.off()

#Find markers
day_type.integrated.markers_merged_1_to_13 <- FindAllMarkers(day_type.integrated, only.pos = TRUE, min.pct = 0.25)
write.table(day_type.integrated.markers_merged_1_to_13, file=paste0(integration_dir,"/All_clusters_marker_genes_full_pc50_res2_merged_1_to_13.txt"), sep = "\t")

#DoHeatmap
#top20 <- day_type.integrated.markers_merged_1_to_13 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
#pdf(file = "All_clusters_Heatmap_marker_genes_pc50_res2_merged_1_to_13.pdf")
#DoHeatmap(day_type.integrated, features = top20$gene, size = 2, draw.lines = TRUE) + theme(axis.text.y = element_text(size = 1.48)) + NoLegend()
#dev.off()

##################################################################################################################################
#Rename into 8 lineages

write.table(levels(day_type.integrated), file=paste0(integration_dir,"1_to_13_levels.txt"), sep = "\t")

lineage.ids <- c("MP", "HCC", "HCC", "Endo", "Endo", "Chol", "HHP", "Endo", "MP", "HSC", "B cell", "NP", "HCC")
names(lineage.ids) <- levels(day_type.integrated)
day_type.integrated <- RenameIdents(day_type.integrated, lineage.ids)
day_type.integrated$lineage_clusters <- Idents(day_type.integrated)

##################################################################################################################################

pdf(paste0(outs_subpath,"/5E_Clustering_tsne_pc50_res2_merged_lineage_label.pdf"), width=12.3, height=10)
DimPlot(day_type.integrated, reduction = "tsne", label=TRUE, pt.size=2, cols = c("#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#66C2A5", "#B3B3B3"))
dev.off()

pdf(paste0(integration_dir,"/Clustering_tsne_pc50_res2_merged_lineage.pdf"), width=12.3, height=10)
DimPlot(day_type.integrated, reduction = "tsne", pt.size=2, cols = c("#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#66C2A5", "#B3B3B3"))
dev.off()

#Find markers
#day_type.integrated.markers_merged_lineage <- FindAllMarkers(day_type.integrated, only.pos = TRUE, min.pct = 0.25)
#write.table(day_type.integrated.markers_merged_lineage, file="All_clusters_marker_genes_full_pc50_res2_merged_lineage.txt", sep = "\t")

#DoHeatmap
#top20.lineage <- day_type.integrated.markers_merged_lineage %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
#pdf(file = "All_clusters_Heatmap_marker_genes_pc50_res2_merged_lineage.pdf")
#DoHeatmap(day_type.integrated, features = top20.lineage$gene, size = 2, draw.lines = TRUE) + theme(axis.text.y = element_text(size = 1.48)) + NoLegend()
#dev.off()

#Proportion of lineages
write.table(table(day_type.integrated@meta.data$lineage_clusters, day_type.integrated@meta.data$day_type), file = paste0(outs_subpath,"/5H_day_type_in_lineage_table.txt"), sep = "\t")

DefaultAssay(day_type.integrated) <- "RNA"

pdf(paste0(integration_dir,"/Stellate_cells_Fibroblasts_VlnPlot.pt_size_0.pdf"), width=15, height=12)
print(VlnPlot(day_type.integrated, features = c("Col1a1","Col1a2","Acta2","Col3a1","Dcn","Tagln","Pdgfrb","Des"), pt.size = 0))
dev.off()

pdf(paste0(integration_dir,"/Endothelial_cells_VlnPlot.pt_size_0.pdf"), width=15, height=12)
print(VlnPlot(day_type.integrated, features = c("Eng","Pecam1","Fcgr2b","Stab2","Cdh5","Kdr","Erg","Icam2"), pt.size = 0))
dev.off()

pdf(paste0(integration_dir,"/Macrophages_Monocytes_VlnPlot.pt_size_0.pdf"), width=15, height=12)
print(VlnPlot(day_type.integrated, features = c("Clec4f","Cd68","Csf1r","Itgam","Adgre1","Cd5l","Marco","Itgax"), pt.size = 0))
dev.off()

pdf(paste0(integration_dir,"/Cholangiocytes_VlnPlot.pt_size_0.pdf"), width=15, height=12)
print(VlnPlot(day_type.integrated, features = c("Sox9","Krt19","Krt7","Mmp7","Epcam","Cldn4","Fxyd3"), pt.size = 0))
dev.off()

pdf(paste0(integration_dir,"/B_cells_VlnPlot.pt_size_0.pdf"), width=15, height=12)
print(VlnPlot(day_type.integrated, features = c("Cd79a","Cd79b","Jchain","Iglc3","Iglc2","Mzb1","Ms4a1"), pt.size = 0))
dev.off()

pdf(paste0(integration_dir,"/HCC_Hepatocyte_VlnPlot.pt_size_0.pdf"), width=15, height=12)
print(VlnPlot(day_type.integrated, features = c("Trf", "Ttr", "Alb", "Serpina1a","Scd1", "Afp", "Gpc3", "Golm1"), pt.size = 0))
dev.off()

pdf(paste0(integration_dir,"/Neutrophils_VlnPlot.pt_size_0.pdf"), width=15, height=12)
print(VlnPlot(day_type.integrated, features = c("Ngp", "Mpo", "Camp", "Cd44","Cd55", "Itgam", "Fcgr3a", "Fcgr2a"), pt.size = 0))
dev.off()

markers <- c("Alb","Ttr","Ass1","Afp","Gpc3","Clec4f","Csf1r","Cd68","Adgre1","Cd5l","Pecam1","Cdh5","Stab2","Kdr","Eng","Krt19","Krt7","Epcam","Sox9","Spp1","Col1a1","Col1a2","Acta2","Dcn","Pdgfrb","Cd79a","Cd79b","Jchain","Iglc2","Mzb1","Ngp","Camp","Mpo","Cxcl2","Csf3r")
pdf(paste0(outs_subpath,"/S5B_Lineage_marker_genes_featureplot.pdf"), height=13,width=12)
FeaturePlot(day_type.integrated, features=markers, ncol=5)
dev.off()

saveRDS(day_type.integrated, paste0(outs_subpath,"/6_samples_integrated_f500_npcs50_res2_renamed.rds"))