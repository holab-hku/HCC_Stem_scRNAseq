library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(tidyr)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_6_and_S6/6ABC_S6AC_HCC_subcluster_UMAP_and_DEG_GO"
outs_subpath <- paste0(dir,"/outs/",figure)
dir.create(outs_subpath, recursive=TRUE)

other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
other_outs_subpath <- paste0(dir,"/outs/",other_figure)

res <- 0.1
num_pcs <- 10

HCC.final.SCT <- readRDS(paste0(other_outs_subpath,"/HCC.final.SCT_pc10_res0.1.rds"))
DefaultAssay(HCC.final.SCT) <- "SCT"


for (set in c("All", "Pos")) {
	if (set == "All") {
		plot_obj <- HCC.final.SCT
		plot_out_path <- paste0(outs_subpath,"/6A_HCC_subclusters_UMAP_",set,".pdf")
	} else if (set == "Pos") {
		plot_obj <- HCC.final.SCT
		Idents(plot_obj) <- "type"
		plot_obj <- subset(plot_obj, idents="Pos")
		Idents(plot_obj) <- "seurat_clusters"
		plot_out_path <- paste0(outs_subpath,"/6C_HCC_subclusters_UMAP_",set,".pdf")
	}
	print(plot_obj)
	
	pdf(plot_out_path, width=13.7, height=10)
	print(DimPlot(plot_obj, reduction = "umap", pt.size=2, cols = c("#377eb8", "#ff7f00", "#4daf4a", "#e41a1c"))+
	  theme(axis.title = element_text(size = 50, color = "black", face = "bold"),
	        axis.text  = element_text(size = 40, color = "black", face = "bold"),
	        legend.position = "none",
	        legend.title = element_text(size = 30,face = "bold")))
	dev.off()
}

Idents(HCC.final.SCT) <- "type"
HCC.final.SCT$type<-factor(HCC.final.SCT$type, levels=c("Pos","Neg")) 
pdf(paste0(outs_subpath,"/6B_HCC_subclusters_UMAP_type_split.pdf"), width=19.4, height=8)
DimPlot(HCC.final.SCT, reduction = "umap", pt.size=2, cols = c("#bebebe", "#fb6a4a"), split.by = "type")+
  theme(axis.title = element_text(size = 50, color = "black", face = "bold"),
        axis.text  = element_text(size = 40, color = "black", face = "bold"),
        legend.position = "none",
        legend.title = element_text(size = 30,face = "bold"))
dev.off()

proportion <- table(HCC.final.SCT@meta.data$seurat_clusters, HCC.final.SCT@meta.data$day_type)
write.table(proportion, paste0(outs_subpath,"/S6A_HCC_subcluster_sample_proportions.txt"),quote=F)

stemness.genes <-  c("Epcam", "Cd24a", "Cd44", "Aldh1a1", "Anpep", "Cd47","Sox9", "Icam1", "Prom1","Lgr5", "Gpc3", "Ly6d")

HCC.final.SCT.POS <- subset(HCC.final.SCT, idents= "Pos")

stemness.genes_featureplot_dir <- paste0(outs_subpath,"/S6C_stemness_genes_FeaturePlot")
dir.create(stemness.genes_featureplot_dir)

DefaultAssay(HCC.final.SCT.POS) <- "SCT"

for (i in stemness.genes){
  pdf(paste0(stemness.genes_featureplot_dir,"/Featureplot_", i, ".pdf"), width=5, height=4)
  f <- FeaturePlot(HCC.final.SCT.POS, features = i, slot = "data", pt.size = 4, order = TRUE, cols = c("#F2E2E2", "#8b0000")) +
   # scale_color_gradientn( colours = c("#F2E2E2", "#8b0000"),  limits = c(0, 2.5)) +
    scale_color_gradientn( colours = c("#F2E2E2", "#8b0000")) +
    theme(plot.title = element_text(color="black", size=50, face = "plain" ),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks  = element_blank(),
          legend.title = element_text(size = 8,face = "bold"),
          legend.text = element_text(color="black", size=20),
          legend.key = element_rect(fill = "NA"))
  print(f)                                                                                                                                                              
  dev.off()}



find_markers_dir <- paste0(outs_subpath,"/FindMarkers")
dir.create(find_markers_dir)

##################################################################################################################################
#DEG for annotation and clinical correlation
##################################################################################################################################

#Find markers 
Idents(HCC.final.SCT) <- "seurat_clusters"
HCC.final.SCT.markers <- FindAllMarkers(HCC.final.SCT, slot = "data", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
write.table(HCC.final.SCT.markers, file = paste0(find_markers_dir,"/HCC.final.SCT.cluster_ALL_min.pct0.1_logfc0.2.txt"), sep = "\t")

#Find markers for Pos/Neg only
Idents(HCC.final.SCT) <- "type"
HCC.final.SCT.Pos <- subset(HCC.final.SCT, idents = "Pos")
HCC.final.SCT.Neg <- subset(HCC.final.SCT, idents = "Neg")
Idents(HCC.final.SCT.Pos) <- "seurat_clusters"
Idents(HCC.final.SCT.Neg) <- "seurat_clusters"

HCC.final.SCT.Pos.markers <- FindAllMarkers(HCC.final.SCT.Pos, slot = "data", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
write.table(HCC.final.SCT.Pos.markers, file = paste0(find_markers_dir,"/HCC.final.SCT.cluster_Pos_min.pct0.1_logfc0.2.txt"), sep = "\t")
HCC.final.SCT.Neg.markers <- FindAllMarkers(HCC.final.SCT.Neg, slot = "data", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
write.table(HCC.final.SCT.Neg.markers, file = paste0(find_markers_dir,"/HCC.final.SCT.cluster_Neg_min.pct0.1_logfc0.2.txt"), sep = "\t")

##################################################################################################################################
#details of "gProfiler_selected.txt"
#GO annotation
#Submit the above three DEG list to gProfiler
#Get a multiquery "gProfiler_mmusculus_17-06-2020_14-09-31__multiquery.txt"
#Highlight liver function related annotations
#export to "gProfiler_selected.txt"
##################################################################################################################################

data<- read.table(paste0(find_markers_dir,"/gProfiler_selected.txt"), header=TRUE, sep="\t")
data<- as.data.frame(data [, 1:14])
data$Annotation<- paste(data$Source, data$Annotation, sep = "_")

#Extract data to plot
#hm<-data[, c("Annotation", "ALL_3","ALL_0","ALL_2", "ALL_1")]
hm<-data[, c("Annotation", "POS_3","POS_0","POS_2", "POS_1")]
#hm<-data[, c("Annotation", "NEG_3","NEG_0","NEG_2", "NEG_1")]
hm$Annotation<-factor(hm$Annotation,levels = rev(hm$Annotation))
hm_final<-melt(hm,id=c("Annotation"))
colnames(hm_final)[3]<-c("pvalue")
hm_final$pvalue <- -log10(hm_final$pvalue)

y<- hm_final$Annotation
x<- hm_final$variable

pdf(paste0(outs_subpath,"/6C_GO_liver_function_Pos.pdf"), width = 10.2, height = 10)
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

##################################################################################################################################
