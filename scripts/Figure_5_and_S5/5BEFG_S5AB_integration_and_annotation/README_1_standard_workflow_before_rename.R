library(dplyr)
library(Seurat)
library(reticulate)
library(ggplot2)
library(cowplot)
library(gprofiler2)
library(RCurl)
options(future.globals.maxSize= 1153433600)

#TD
#setwd("D:/Lena/scRNAseq")
dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
data_subpath <- paste0(dir,"/data/",figure)
outs_subpath <- paste0(dir,"/outs/",figure)

dir.create(outs_subpath, recursive=TRUE)

#Read the data for the 6 samples
D3N.data <- Read10X(data.dir = paste0(data_subpath,"/Cell_Ranger_filtered_matrices/Day3/cellranger_v3/Neg"))
D3N.meta <- read.table(paste0(data_subpath,"/Cell_Ranger_filtered_matrices/Day3/cellranger_v3/Neg/metadata"), row.names=1)
D3N <- CreateSeuratObject(counts = D3N.data, project = "Stem_scRNAseq", meta.data = D3N.meta)

D3P.data <- Read10X(data.dir = paste0(data_subpath,"/Cell_Ranger_filtered_matrices/Day3/cellranger_v3/Pos"))
D3P.meta <- read.table(paste0(data_subpath,"/Cell_Ranger_filtered_matrices/Day3/cellranger_v3/Pos/metadata"), row.names=1)
D3P <- CreateSeuratObject(counts = D3P.data, project = "Stem_scRNAseq", meta.data = D3P.meta)

D10N.data <- Read10X(data.dir = paste0(data_subpath,"/Cell_Ranger_filtered_matrices/Day10/cellranger_v3/Neg"))
D10N.meta <- read.table(paste0(data_subpath,"/Cell_Ranger_filtered_matrices/Day10/cellranger_v3/Neg/metadata"), row.names=1)
D10N <- CreateSeuratObject(counts = D10N.data, project = "Stem_scRNAseq", meta.data = D10N.meta)

D10P.data <- Read10X(data.dir = paste0(data_subpath,"/Cell_Ranger_filtered_matrices/Day10/cellranger_v3/Pos"))
D10P.meta <- read.table(paste0(data_subpath,"/Cell_Ranger_filtered_matrices/Day10/cellranger_v3/Pos/metadata"), row.names=1)
D10P <- CreateSeuratObject(counts = D10P.data, project = "Stem_scRNAseq", meta.data = D10P.meta)

D30N.data <- Read10X(data.dir = paste0(data_subpath,"/Cell_Ranger_filtered_matrices/Day30/cellranger_v3/Neg"))
D30N.meta <- read.table(paste0(data_subpath,"/Cell_Ranger_filtered_matrices/Day30/cellranger_v3/Neg/metadata"), row.names=1)
D30N <- CreateSeuratObject(counts = D30N.data, project = "Stem_scRNAseq", meta.data = D30N.meta)

D30P.data <- Read10X(data.dir = paste0(data_subpath,"/Cell_Ranger_filtered_matrices/Day30/cellranger_v3/Pos"))
D30P.meta <- read.table(paste0(data_subpath,"/Cell_Ranger_filtered_matrices/Day30/cellranger_v3/Pos/metadata"), row.names=1)
D30P <- CreateSeuratObject(counts = D30P.data, project = "Stem_scRNAseq", meta.data = D30P.meta)

#TD
#setwd("D:/Lena/scRNAseq/6 samples integration/standard workflow f500")
Neg_data <- merge(x=D3N, y=c(D10N,D30N))
Pos_data <- merge(x=D3P, y=c(D10P,D30P))

#According to the Seurat - Guided Clustering Tutorial https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
#"Standard pre-processing workflow"
#Filter cells with > X_thresh% mitochondria counts, and cells that are outside of Lower_bound_nFeature_RNA and Upper_bound_nFeature_RNA
X_thresh <- 10
Lower_bound_nFeature_RNA <- 200
Upper_bound_nFeature_RNA <- 5000

QC_dir <- paste0(outs_subpath,"/QC")
dir.create(QC_dir)

Neg_data[["percent.mt"]] <- PercentageFeatureSet(Neg_data, pattern = "^mt-")
pdf(paste0(QC_dir,"/Neg_samples_QC_violin_plot.pdf"), width=14, height=9)
VlnPlot(Neg_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
Pos_data[["percent.mt"]] <- PercentageFeatureSet(Pos_data, pattern = "^mt-")
pdf(paste0(QC_dir,"/Pos_samples_QC_violin_plot.pdf"), width=14, height=9)
VlnPlot(Pos_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf(paste0(QC_dir,"/Neg_samples_FeatureScatter_plot.pdf"), width=14, height=9)
plot1 <- FeatureScatter(Neg_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Neg_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

pdf(paste0(QC_dir,"/Pos_samples_FeatureScatter_plot.pdf"), width=14, height=9)
plot1 <- FeatureScatter(Pos_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Pos_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

plot1 <- NULL
plot2 <- NULL

Neg_data
Neg_data <- subset(Neg_data, subset = nFeature_RNA > Lower_bound_nFeature_RNA & nFeature_RNA < Upper_bound_nFeature_RNA  & percent.mt < X_thresh)
Neg_data

Pos_data
Pos_data <- subset(Pos_data, subset = nFeature_RNA > Lower_bound_nFeature_RNA & nFeature_RNA < Upper_bound_nFeature_RNA  & percent.mt < X_thresh)
Pos_data

types <- merge(x=Neg_data, y=c(Pos_data))

#Set the variables no longer in use to NULL to save memory
D3N.data <- NULL
D3N.meta <- NULL
D3N <- NULL

D3P.data <- NULL
D3P.meta <- NULL
D3P <- NULL

D10N.data <- NULL
D10N.meta <- NULL
D10N <- NULL

D10P.data <- NULL
D10P.meta <- NULL
D10P <- NULL

D30N.data <- NULL
D30N.meta <- NULL
D30N <- NULL

D30P.data <- NULL
D30P.meta <- NULL
D30P <- NULL

Neg_data <- NULL
Pos_data <- NULL

##################################################################################################################################
#Proceed according to the Seurat data integration vignette (Standard Workflow): https://satijalab.org/seurat/v3.0/integration.html 
##################################################################################################################################

day_type.list <- SplitObject(types, split.by = "day_type")

types <- NULL

for (i in 1:length(day_type.list)) {
  day_type.list[[i]] <- NormalizeData(day_type.list[[i]], verbose = FALSE)
  day_type.list[[i]] <- FindVariableFeatures(day_type.list[[i]], selection.method = "vst", 
                                             nfeatures = 500, verbose = FALSE)
}

reference.list <- day_type.list[c("Day3_Neg", "Day3_Pos", "Day10_Neg", "Day10_Pos", "Day30_Neg", "Day30_Pos")]

day_type.list <- NULL

day_type.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

reference.list <- NULL

to_integrate <- Reduce(intersect, lapply(day_type.anchors@object.list, rownames))

day_type.integrated <- IntegrateData(anchorset = day_type.anchors, features.to.integrate = to_integrate, dims = 1:30)

day_type.anchors <- NULL
to_integrate <- NULL

DefaultAssay(day_type.integrated) <- "integrated"
##################################################################################################################################
#Cell cycle scoring (no regression)
##################################################################################################################################
s.genes <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Mlf1","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8")
g2m.genes <- c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")

day_type.integrated <- CellCycleScoring(day_type.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

day_type.integrated <- ScaleData(day_type.integrated, features = rownames(day_type.integrated))


integration_dir <- paste0(outs_subpath,"/integration")
dir.create(integration_dir)

#saveRDS
saveRDS(day_type.integrated,file=paste0(integration_dir,"/6_samples_integrated_f500.rds"))

##################################################################################################################################
#day_type.integrated <- RunPCA(day_type.integrated, npcs = 50, verbose = FALSE)
##Jacksaw to determine number of dimension, use the same number of dims as number of PCs in RunPCA
#day_type.integrated <- JackStraw(day_type.integrated,dims = 50)
#day_type.integrated <- ScoreJackStraw(day_type.integrated,dims = 1:50)
#pdf("JackStraw_6_samples_integrated_f500.pdf", width=14, height=9)
#j1<-JackStrawPlot(day_type.integrated, dims=1:50)
#plot_grid(j1)
#dev.off()

##Look up marker genes from 2 to N (N depends on JackStrawPlot)

# for (i in 2:50) {
#   day_type.integrated <- readRDS(file="6_samples_integrated_f500.rds")
#   DefaultAssay(day_type.integrated) <- "integrated"
#   
#   day_type.integrated <- RunPCA(day_type.integrated, npcs = i, verbose = FALSE)
#   day_type.integrated <- RunTSNE(day_type.integrated, reduction = "pca", dims = 1:i)
#   
#   pdf(paste0("Clustering_tsne_pc", i, "_split.pdf"), width=25, height=10)
#   print(DimPlot(day_type.integrated, reduction = "tsne", split.by= "type", pt.size=1))
#   dev.off()
#   
#   pdf(paste0("Clustering_tsne_pc", i, ".pdf"), width=12, height=10)
#   print(DimPlot(day_type.integrated, reduction = "tsne", label=TRUE, pt.size=1))
#   dev.off()
#   
#   DefaultAssay(day_type.integrated) <- "RNA"
#   
#   pdf(paste0("Stellate_Fibroblast_Mesenchyme_FeaturePlot_pc", i, ".pdf"), width=15, height=12)
#   print(FeaturePlot(day_type.integrated, features = c("Col1a1","Col1a2","Acta2","Col3a1","Dcn","Tagln","Pdgfrb","Des"), reduction = "tsne"))
#   dev.off()
#   
#   pdf(paste0("Endothelial_cells_FeaturePlot_pc", i, ".pdf"), width=15, height=12)
#   print(FeaturePlot(day_type.integrated, features = c("Eng","Pecam1","Fcgr2b","Stab2","Cdh5","Kdr","Erg","Icam2"), reduction = "tsne"))
#   dev.off()
#   
#   pdf(paste0("Macrophages_FeaturePlot_pc", i, ".pdf"), width=15, height=12)
#   print(FeaturePlot(day_type.integrated, features = c("Clec4f","Cd68","Csf1r","Itgam","Adgre1","Cd5l","Marco","Itgax"), reduction = "tsne"))
#   dev.off()
#   
#   pdf(paste0("Cholangiocyte_FeaturePlot_pc", i, ".pdf"), width=15, height=12)
#   print(FeaturePlot(day_type.integrated, features = c("Sox9","Krt19","Krt7","Mmp7","Epcam","Cldn4","Fxyd3"), reduction = "tsne"))
#   dev.off()
#   
#   pdf(paste0("B_cells_FeaturePlot_pc", i, ".pdf"), width=15, height=12)
#   print(FeaturePlot(day_type.integrated, features = c("Cd79a","Cd79b","Jchain","Iglc3","Iglc2","Mzb1","Ms4a1"), reduction = "tsne"))
#   dev.off()
#   
#   pdf(paste0("Dividing_cells_FeaturePlot_pc", i, ".pdf"), width=15, height=12)
#   print(FeaturePlot(day_type.integrated, features = c("Mki67","Stmn1","Top2a","Ccna2","Ccnb2"), reduction = "tsne"))
#   dev.off()
#   
#   pdf(paste0("Hepatocyte_HCC_FeaturePlot_pc", i, ".pdf"), width=15, height=12)
#   print(FeaturePlot(day_type.integrated, features = c("Trf", "Ttr", "Alb", "Serpina1a","Scd1", "Afp", "Gpc3", "Golm1"), reduction = "tsne"))
#   dev.off()
#   }

##################################################################################################################################
#Test and define resolution
#Look up resolution from (2, 4, 0.5)
#day_type.integrated <- readRDS(file="6_samples_integrated_f500.rds")
#DefaultAssay(day_type.integrated) <- "integrated"
##################################################################################################################################

# num_pcs <- 50
# res <- 2
# day_type.integrated <- RunPCA(day_type.integrated, npcs = num_pcs, verbose = FALSE)
# day_type.integrated <- FindNeighbors(day_type.integrated, dims = 1:num_pcs)
# day_type.integrated <- FindClusters(day_type.integrated, resolution = res)
# day_type.integrated <- RunTSNE(day_type.integrated, reduction = "pca", dims = 1:num_pcs)
# 
# pdf(paste0("Clustering_tsne_pc", num_pcs, "_res", res, "_split.pdf"), width=25, height=10)
# DimPlot(day_type.integrated, reduction = "tsne", split.by= "type", pt.size=1)
# dev.off()
# 
# pdf(paste0("Clustering_tsne_pc", num_pcs, "_res", res, ".pdf"), width=12, height=10)
# DimPlot(day_type.integrated, reduction = "tsne", label=TRUE, pt.size=1)
# dev.off()
# 
# saveRDS(day_type.integrated,file="6_samples_integrated_f500_npcs50_res2.rds")

day_type.integrated <- readRDS(file=paste0(outs_subpath,"/6_samples_integrated_f500_npcs50_res2.rds"))

# #Find markers
# day_type.integrated.markers <- FindAllMarkers(day_type.integrated, only.pos = TRUE, min.pct = 0.25)
# write.table(day_type.integrated.markers, file="All_clusters_marker_genes_pc50_res2.txt", sep = "\t")
# 
# #DoHeatmap
# top20 <- day_type.integrated.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
# pdf(file = "All_clusters_marker_genes_pc50_res2.pdf")
# DoHeatmap(day_type.integrated, features = top20$gene, size = 3, draw.lines = TRUE) + theme(axis.text.y = element_text(size = 1.48)) + NoLegend()
# dev.off()
# 
# pdf(file = "All_clusters_marker_genes_pc50_res2_Pos.pdf")
# DoHeatmap(subset(day_type.integrated, subset = type =="Pos"), features = top20$gene, size = 3, draw.lines = TRUE) + theme(axis.text.y = element_text(size = 1.48)) + NoLegend()
# dev.off()
# pdf(file = "All_clusters_marker_genes_pc50_res2_Neg.pdf")
# DoHeatmap(subset(day_type.integrated, subset = type =="Neg"), features = top20$gene, size = 3, draw.lines = TRUE) + theme(axis.text.y = element_text(size = 1.48)) + NoLegend()
# dev.off()

#RenameIdents into Cluster_type
day_type.integrated$Cluster_type<-paste(day_type.integrated$type, day_type.integrated$seurat_clusters,  sep = "_")
Idents(day_type.integrated) <- "Cluster_type"
# write.table(levels(day_type.integrated), file="Cluster_type_levels.txt", sep = "\t")

#BuildClusterTree assists clusters merging
day_type.integrated <- BuildClusterTree(day_type.integrated, reorder = F, reorder.numeric = F)
pdf(paste0(integration_dir,"/ClusterTree_raw_clusters.pdf"), width=25, height=10)
PlotClusterTree(day_type.integrated)
dev.off()

#Find markers for Cluster_type
day_type.integrated.markers <- FindAllMarkers(day_type.integrated, only.pos = TRUE, min.pct = 0.25)
write.table(day_type.integrated.markers, file=paste0(integration_dir,"/All_clusters_marker_genes_pc50_res2_Cluster_type.txt"), sep = "\t")

##################################################################################################################################
#Merging clusters according to 
#1)CNV "CNV_levels_Pos_Neg_0-35.pdf" 
#2)clustertree "ClusterTree_raw_clusters.pdf" 
#3)markers "All_clusters_marker_genes_pc50_res2_Cluster_type.txt"
#4)tSNE "Clustering_tsne_pc50_res2.pdf" & "Clustering_tsne_pc50_res2_split.pdf"

#34_Neg 26_Neg 26_Pos = "Cholangiocytes": low CNV, markers "Epcam" "Spp1", same cluster

#33_Neg = "Macrophages/Monocytes" or "Endothelial_cells": low CNV, markers "Csf1r" "Clec4f" "Ptprb" "Kdr" "Cdh5"
#13_Neg = "Endothelial_cells": low CNV, markers "Stab2" "Eng" "Ptprb" "Cdh5"
#29_Pos 33_Pos 24_Pos 7_Pos 14_Pos 13_Pos 11_Pos: no CNV data(too few cells), mixed markers, close to 33_Neg & 13_Neg
#The above clusters = "Endothelial_cells_1"

#19_Neg 19_Pos = "Stellate_cells/Fibroblasts": low CNV or no CNV data(19_Pos too few cells), markers "Acta2" "Dcn" "Col1a1" "Vim"

#14_Neg = "Macrophages/Monocytes": low CNV, markers "Csf1r" "Clec4f" "Cd68"
#29_Neg = "Macrophages/Monocytes": low CNV, markers "Ccl5r", close to 14_Neg 
#The above clusters = "Macrophages/Monocytes_1"

#25_Neg 25_Pos = "B_cells/Plasma_cells": low CNV or no CNV data(25_Pos too few cells), markers "Igxx" "Jchain" "Mzb1"

#24_Neg = "Macrophages/Monocytes": low CNV, markers "H2-Ab1" "H2-Aa" "Lilr4b"
#0_Neg: low CNV, mixed markers, close to 24_Neg
#The above clusters = "Macrophages/Monocytes_2"

#11_Neg = "Endothelial_cells_2": low CNV, markers "Pecam1" 

#35_Neg 32_Neg = "Neutrophil": low CNV, markers "Ngp" "Mpo" 

#7_Neg = "Endothelial_cells": low CNV, markers "Pecam1" "Stab1"
#3_Neg 0_Neg: low CNV, mixed markers, close to 7_Neg 
#The above clusters = "Endothelial_cells_3"

#The rest of cells, separate into 4 groups according to clusters of clustertree
#15_Neg 15_Pos = "HCC/Hepatocytes_1": low CNV, markers "Spp1"
#10_Pos 10_Neg 5_Pos 5_Neg 22_Pos 22_Neg 18_Pos 18_Neg 17_Pos 17_Neg 12_Pos 12_Neg 9_Pos 9_Neg 27_Pos 27_Neg = "HCC/Hepatocytes_2": high CNV, markers "Trf" "Cyp2e1" "Spp1" "Gss" "Scd1" "Cypxx" "Gstm" "Sqstm1"
#4_Pos 4_Neg 28_Pos 28_Neg 21_Pos 8_Pos 8_Neg 30_Pos 30_Neg 20_Pos 23_Pos 31_Pos 2_Pos 2_Neg 16_Pos 23_Neg 6_Pos 6_Neg 1_Pos 1_Neg 16_Neg 0_Pos = "HCC/Hepatocytes_3": high CNV, markers "Orm1" "C3" "Fah" "Gnmt" "Cyp2e1" "Afp" "Ass1" "Scd1" "Trf" G2M/S markers "Alb" Coagulation/Complement factors "Apoa4"
#31_Neg 20_Neg 21_Neg 3_Pos = "HCC/Hepatocytes_4": medium CNV, G2M/S markers "Gstm" "Hp" "Tff3"
##################################################################################################################################

#Merging clusters (raw)
newcluster.ids <- c("Macrophages/Monocytes_1", "HCC/Hepatocytes_2", "HCC/Hepatocytes_3", "Endothelial_cells_3", "HCC/Hepatocytes_3", "HCC/Hepatocytes_4", "Endothelial_cells_2", "Cholangiocytes", "HCC/Hepatocytes_1", "Endothelial_cells_1", "HCC/Hepatocytes_3", "Cholangiocytes", "HCC/Hepatocytes_3", "Macrophages/Monocytes_2", "Macrophages/Monocytes_2", "HCC/Hepatocytes_3", "HCC/Hepatocytes_2", "HCC/Hepatocytes_4", "Endothelial_cells_3", "HCC/Hepatocytes_3", "HCC/Hepatocytes_2", "HCC/Hepatocytes_2", "HCC/Hepatocytes_3", "Stellate_cells/Fibroblasts", "B_cells/Plasma_cells", "HCC/Hepatocytes_3", "HCC/Hepatocytes_2", "HCC/Hepatocytes_2", "Neutrophil", "HCC/Hepatocytes_2", "Neutrophil", "Endothelial_cells_1", "Macrophages/Monocytes_1", "HCC/Hepatocytes_4", "HCC/Hepatocytes_3", "HCC/Hepatocytes_2", "HCC/Hepatocytes_3", "HCC/Hepatocytes_2", "HCC/Hepatocytes_3", "HCC/Hepatocytes_3", "HCC/Hepatocytes_3", "HCC/Hepatocytes_2", "HCC/Hepatocytes_2", "HCC/Hepatocytes_2", "HCC/Hepatocytes_2", "HCC/Hepatocytes_3", "HCC/Hepatocytes_2", "HCC/Hepatocytes_3", "HCC/Hepatocytes_4", "HCC/Hepatocytes_2", "HCC/Hepatocytes_3", "Cholangiocytes", "HCC/Hepatocytes_3", "HCC/Hepatocytes_3", "HCC/Hepatocytes_3", "HCC/Hepatocytes_1", "HCC/Hepatocytes_3", "HCC/Hepatocytes_2", "HCC/Hepatocytes_3", "HCC/Hepatocytes_3", "Stellate_cells/Fibroblasts", "Endothelial_cells_1", "Endothelial_cells_1", "Endothelial_cells_1", "Endothelial_cells_1", "Endothelial_cells_1", "Endothelial_cells_1", "B_cells/Plasma_cells", "Endothelial_cells_1")
names(newcluster.ids) <- levels(day_type.integrated)
day_type.integrated <- RenameIdents(day_type.integrated, newcluster.ids)
day_type.integrated$merged_clusters <- Idents(day_type.integrated)

pdf(paste0(integration_dir,"/Clustering_tsne_pc50_res2_merged_split.pdf"), width=25, height=10)
DimPlot(day_type.integrated, reduction = "tsne", split.by= "type", pt.size=1)
dev.off()

pdf(paste0(integration_dir,"/Clustering_tsne_pc50_res2_merged.pdf"), width=12, height=10)
DimPlot(day_type.integrated, reduction = "tsne", label=TRUE, pt.size=1)
dev.off()

table(day_type.integrated@meta.data$merged_clusters, day_type.integrated@meta.data$type)

#Find markers
day_type.integrated.markers_merged <- FindAllMarkers(day_type.integrated, only.pos = TRUE, min.pct = 0.25)
write.table(day_type.integrated.markers_merged, file=paste0(integration_dir,"/All_clusters_marker_genes_full_pc50_res2_merged.txt"), sep = "\t")

# #DoHeatmap
# top20 <- day_type.integrated.markers_merged %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
# pdf(file = "All_clusters_Heatmap_marker_genes_pc50_res2_merged.pdf")
# DoHeatmap(day_type.integrated, features = top20$gene, size = 2, draw.lines = TRUE) + theme(axis.text.y = element_text(size = 1.48)) + NoLegend()
# dev.off()

DefaultAssay(day_type.integrated) <- "RNA"

pdf(paste0(integration_dir,"/Stellate_cells_Fibroblasts_VlnPlot.pdf"), width=15, height=12)
print(VlnPlot(day_type.integrated, features = c("Col1a1","Col1a2","Acta2","Col3a1","Dcn","Tagln","Pdgfrb","Des"), pt.size = 0.1))
dev.off()

pdf(paste0(integration_dir,"/Endothelial_cells_VlnPlot.pdf"), width=15, height=12)
print(VlnPlot(day_type.integrated, features = c("Eng","Pecam1","Fcgr2b","Stab2","Cdh5","Kdr","Erg","Icam2"), pt.size = 0.1))
dev.off()

pdf(paste0(integration_dir,"/Macrophages_Monocytes_VlnPlot.pdf"), width=15, height=12)
print(VlnPlot(day_type.integrated, features = c("Clec4f","Cd68","Csf1r","Itgam","Adgre1","Cd5l","Marco","Itgax"), pt.size = 0.1))
dev.off()

pdf(paste0(integration_dir,"/Cholangiocytes_VlnPlot.pdf"), width=15, height=12)
print(VlnPlot(day_type.integrated, features = c("Sox9","Krt19","Krt7","Mmp7","Epcam","Cldn4","Fxyd3"), pt.size = 0.1))
dev.off()

pdf(paste0(integration_dir,"/B_cells_VlnPlot.pdf"), width=15, height=12)
print(VlnPlot(day_type.integrated, features = c("Cd79a","Cd79b","Jchain","Iglc3","Iglc2","Mzb1","Ms4a1"), pt.size = 0.1))
dev.off()

pdf(paste0(integration_dir,"/HCC_Hepatocyte_VlnPlot.pdf"), width=15, height=12)
print(VlnPlot(day_type.integrated, features = c("Trf", "Ttr", "Alb", "Serpina1a","Scd1", "Afp", "Gpc3", "Golm1"), pt.size = 0.1))
dev.off()

pdf(paste0(integration_dir,"/Neutrophils_VlnPlot.pdf"), width=15, height=12)
print(VlnPlot(day_type.integrated, features = c("Ngp", "Mpo", "Camp", "Cd44","Cd55", "Itgam", "Fcgr3a", "Fcgr2a"), pt.size = 0.1))
dev.off()

##################################################################################################################################