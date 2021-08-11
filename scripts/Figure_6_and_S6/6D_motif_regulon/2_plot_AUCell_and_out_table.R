library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggsignif)
library(ggpubr)

rm(list=ls())

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_6_and_S6/6D_motif_regulon"
outs_subpath <- paste0(dir,"/outs/",figure)
dir.create(outs_subpath, recursive=TRUE)

other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
other_outs_subpath <- paste0(dir,"/outs/",other_figure)

param <- "pc10_res0.1"

scRNA <- readRDS(paste0(other_outs_subpath,"/HCC.final.SCT_",param,".rds"))

##read regulonAUCell
Regulon <- readRDS(paste0(outs_subpath,"/3.4_regulonAUC.Rds"))
Regulon <- Regulon@assays@data@listData$AUC
dim(Regulon)

# Regulon_high <- Regulon[!grepl("extended", rownames(Regulon)),] #Regulon_high or Regulon_low
# row.names(Regulon_high) <- str_split(row.names(Regulon_high), "_", simplify = T)[,1]
Regulon_low <- Regulon[grepl("extended", rownames(Regulon)),]
row.names(Regulon_low) <- str_split(row.names(Regulon_low), "_", simplify = T)[,1]

#filter out TF of low AUC
# Regulon_high <- Regulon_high[rowSums(Regulon_high>0)>0.5*length(colnames(Regulon_high)),]
# dim(Regulon_high)
Regulon_low <- Regulon_low[rowSums(Regulon_low>0)>0.5*length(colnames(Regulon_low)),]
dim(Regulon_low)

# scRNA[["Regulon"]] <- CreateAssayObject(counts = Regulon_high) 
scRNA[["Regulon"]] <- CreateAssayObject(counts = Regulon_low) 

DefaultAssay(scRNA) <- "Regulon"
scRNA$seurat_clusters <- factor(scRNA$seurat_clusters, levels = c("3", "0", "2", "1"))
scRNA$day <- factor(scRNA$day, levels = c("Day3", "Day10", "Day30"))
dim(scRNA)

################## Pos subset
Idents(scRNA) <- "type"
Pos <- subset(scRNA, idents= "Pos")
DefaultAssay(Pos) <- "Regulon"

for (TF in c("Hnf4a", "Cebpa", "Foxa3")){
  VlnPlot(Pos, features = TF, 
          group.by = "seurat_clusters", 
          pt.size=0,
          cols = c("#e41a1c", "#377eb8", "#4daf4a", "#ff7f00"))+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA, show.legend = FALSE)+
    labs(x="Cluster", y = "AUCell", title = paste0(TF, " motif regulon")) +
    scale_x_discrete(labels=c("3" = "C3", "0" = "C0",
                              "2" = "C2", "1" = "C1")) +
    theme(axis.text.x = element_text(face="plain", color="black", 
                                     size=10, angle=0, hjust = 0.5),
          axis.text.y = element_text(face="plain", color="black", 
                                     size=10, angle=0),
          axis.title.x = element_text(face="bold", color="black", 
                                      size=12, angle=0),
          axis.title.y = element_text(face="bold", color="black", 
                                      size=12, angle=90),
          plot.title = element_text(color="black", size=14, face = "bold", hjust = 1), 
          legend.position = "none")+
    geom_signif(comparisons = list(c("C3", "C0"),
                                   c("C2", "C1")))
  
  ggsave(paste0(outs_subpath,"/6D_",TF, "_AUCell_VlnPlot.pdf"), width = 2.3, height = 4, device="pdf")
}

# write.table(row.names(Pos[["Regulon"]]), "TF_list.csv", row.names = F, col.names = F)


#DotPlot for TF
##############################
DefaultAssay(Pos) <- "SCT"
Idents(Pos) <- "seurat_clusters"

DotPlot(Pos, features = c("Hnf4a", "Cebpa", "Foxa3"), 
               scale.by = "size") + 
       scale_colour_gradient(low = "#F2E2E2", high = "#8b0000")+
       coord_flip()
ggsave(paste0(outs_subpath,"/6D_TF_expr_Dotplot.pdf"), width = 5, height = 2.5, device="pdf")

################## Regulon full list (Table S4)
DefaultAssay(Pos) <- "Regulon"
Idents(Pos) <- "seurat_clusters"

pos_regulon = FindAllMarkers(Pos, slot = 'counts', logfc.threshold = 0, only.pos = T, return.thresh = 0.05)
write.table(pos_regulon, paste0(outs_subpath,"/Regulon_full_list.txt"))
