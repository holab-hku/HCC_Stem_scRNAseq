rm(list = ls())
options(stringsAsFactors = F)

library(dplyr)
library(Seurat)
library(reticulate)
library(ggplot2)
library(cowplot)
library(gprofiler2)
#library(gProfileR)
library(RCurl)
library(RColorBrewer)
options(future.globals.maxSize = 4500 * 1024^2)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
data_subpath <- paste0(dir,"/data/",figure)
outs_subpath <- paste0(dir,"/outs/",figure)

integration_dir <- paste0(outs_subpath,"/integration")

day_type.integrated <- readRDS(paste0(outs_subpath,"/6_samples_integrated_f500_npcs50_res2_renamed.rds"))

##################################################################################################################################
#subset HCC cells
day_type.integrated.HCC <- subset(day_type.integrated, idents="HCC")
                                  
metadata <- select(day_type.integrated.HCC@meta.data, type, day, day_type)
                                  
#get the raw counts data and create seurat object
HCC.raw_data <- GetAssayData(day_type.integrated.HCC, slot = "counts", assay = "RNA")
types.HCC <- CreateSeuratObject(counts = HCC.raw_data, project = "Stem_scRNAseq", meta.data = metadata)
                                  
day_type.list <- SplitObject(types.HCC, split.by = "day_type")
                                  
day_type.integrated <- NULL
day_type.integrated.HCC <- NULL
types.HCC <- NULL
HCC.raw_data <- NULL
metadata <- NULL

##################################################################################################################################
#Proceed according to the Seurat data integration vignette (Standard Workflow): https://satijalab.org/seurat/v3.0/integration.html 
##################################################################################################################################

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

HCC.standard.integrated <- IntegrateData(anchorset = day_type.anchors, features.to.integrate = to_integrate, dims = 1:30)

day_type.anchors <- NULL
to_integrate <- NULL

HCC.standard.integrated <- ScaleData(HCC.standard.integrated, features = rownames(HCC.standard.integrated))


num_pcs <- 50
res <- 1
HCC.standard.integrated <- RunPCA(HCC.standard.integrated, npcs = num_pcs, verbose = FALSE)
HCC.standard.integrated <- FindNeighbors(HCC.standard.integrated, dims = 1:num_pcs)
HCC.standard.integrated <- FindClusters(HCC.standard.integrated, resolution = res)
HCC.standard.integrated <- RunTSNE(HCC.standard.integrated, reduction = "pca", dims = 1:num_pcs)

pdf(paste0(integration_dir,"/HCC.standard_f500_raw_tsne_pc", num_pcs, "_res", res, "_split.pdf"), width=25, height=10)
DimPlot(HCC.standard.integrated, reduction = "tsne", split.by= "type", pt.size=1)
dev.off()

pdf(paste0(integration_dir,"/HCC.standard_f500_raw_tsne_pc", num_pcs, "_res", res, ".pdf"), width=12, height=10)
DimPlot(HCC.standard.integrated, reduction = "tsne", label=TRUE, pt.size=1)
dev.off()

DefaultAssay(HCC.standard.integrated) <- "RNA"

#Find markers
HCC.standard.integrated.markers <- FindAllMarkers(HCC.standard.integrated, only.pos = TRUE, min.pct = 0.25)
write.table(HCC.standard.integrated.markers, file=paste0(integration_dir,"/HCC.standard_f500_raw_pc50_res1.txt"), sep = "\t")

#saveRDS
saveRDS(HCC.standard.integrated,file=paste0(outs_subpath,"/HCC.standard_f500_raw_npcs50_res1.rds"))
#HCC.standard.integrated <- readRDS(file="HCC.standard_f500_raw_npcs50_res1.rds")

##################################################################################################################################
#Exclude cells of no interest (6/12 = dividing HCC, 14 = non-HCC, 16 = virus-responding)
day_type.integrated.HCC.final <- subset(HCC.standard.integrated, idents=c("0", "1", "2", "3", "4", "5", "7", "8", "9", "10", "11", "13", "15", "17", "18"))

metadata <- select(day_type.integrated.HCC.final@meta.data, type, day, day_type)

#get the raw counts data and create seurat object
HCC.raw_data <- GetAssayData(day_type.integrated.HCC.final, slot = "counts", assay = "RNA")
types.HCC <- CreateSeuratObject(counts = HCC.raw_data, project = "Stem_scRNAseq", meta.data = metadata)

day_type.list <- SplitObject(types.HCC, split.by = "day_type")
#this step is: day_type.list <- SplitObject(types, split.by = "day_type") in all the scripts so far

saveRDS(day_type.list,file=paste0(outs_subpath,"/HCC.standard_f500_final_raw.rds"))
#day_type.list <- readRDS("HCC.standard_f500_final_raw.rds")

HCC.standard.integrated <- NULL
day_type.integrated.HCC.final <- NULL
types.HCC <- NULL
HCC.raw_data <- NULL
metadata <- NULL

##################################################################################################################################
#Proceed according to the Seurat data integration vignette (SCT): https://satijalab.org/seurat/v3.0/integration.html 
##################################################################################################################################

for (i in 1:length(day_type.list)) {
  day_type.list[[i]] <- SCTransform(day_type.list[[i]], verbose = FALSE, variable.features.n = 10000, return.only.var.genes = FALSE)
}


day_type.features <- SelectIntegrationFeatures(object.list = day_type.list, nfeatures = 10000)
day_type.list <- PrepSCTIntegration(object.list = day_type.list, anchor.features = day_type.features, verbose = FALSE)
day_type.anchors <- FindIntegrationAnchors(object.list = day_type.list, normalization.method = "SCT", anchor.features = day_type.features)

day_type.features <- NULL
day_type.list <- NULL

HCC.final.SCT <- IntegrateData(anchorset = day_type.anchors, normalization.method = "SCT")

day_type.anchors <- NULL

HCC.final.SCT <- ScaleData(HCC.final.SCT, features = rownames(HCC.final.SCT))

#saveRDS
saveRDS(HCC.final.SCT,file=paste0(outs_subpath,"/HCC.standard.final.integrated_f500_SCT.rds"))
#HCC.final.SCT <- readRDS ("HCC.standard.final.integrated_f500_SCT.rds")

#HCC.final.SCT <- RunPCA(HCC.final.SCT, npcs = 50, verbose = FALSE)
#Jacksaw to determine number of dimension, use the same number of dims as number of PCs in RunPCA
#HCC.final.SCT <- JackStraw(HCC.final.SCT,dims = 50)
#HCC.final.SCT <- ScoreJackStraw(HCC.final.SCT,dims = 1:50)
#pdf("JackStraw_HCC.final.SCT.pdf", width=14, height=9)
#j1<-JackStrawPlot(HCC.final.SCT, dims=1:50)
#plot_grid(j1)
#dev.off()

HCC.final.SCT$day_type<-factor(HCC.final.SCT$day_type, levels=c("Day3_Pos","Day10_Pos","Day30_Pos","Day3_Neg","Day10_Neg","Day30_Neg")) 

res <- 0.1
#num_pcs <- 20
num_pcs <- 10

HCC.final.SCT <- RunPCA(HCC.final.SCT, npcs = num_pcs, verbose = FALSE)
HCC.final.SCT <- FindNeighbors(HCC.final.SCT, dims = 1:num_pcs)
HCC.final.SCT <- FindClusters(HCC.final.SCT, resolution = res)
HCC.final.SCT <- RunTSNE(HCC.final.SCT, reduction = "pca", dims = 1:num_pcs)
HCC.final.SCT <- RunUMAP(HCC.final.SCT, reduction = "pca", dims = 1:num_pcs)

#saveRDS
saveRDS(HCC.final.SCT,file=paste0(outs_subpath,"/HCC.final.SCT_pc", num_pcs, "_res", res, ".rds"))

pdf(paste0(integration_dir,"/HCC.final.SCT_umap_pc", num_pcs, "_res", res, "_split.pdf"), width=20, height=10)
DimPlot(HCC.final.SCT, reduction = "umap", split.by= "day_type", pt.size=1, ncol = 3)
dev.off()

pdf(paste0(integration_dir,"/HCC.final.SCT_umap_pc", num_pcs, "_res", res, ".pdf"), width=12, height=10)
DimPlot(HCC.final.SCT, reduction = "umap", label=TRUE, pt.size=1)
dev.off()

#Proportion of clusters
#write.table(table(HCC.final.SCT@meta.data$seurat_clusters, HCC.final.SCT@meta.data$day_type), file = "day_type_in_HCC_clusters", sep = "\t")
