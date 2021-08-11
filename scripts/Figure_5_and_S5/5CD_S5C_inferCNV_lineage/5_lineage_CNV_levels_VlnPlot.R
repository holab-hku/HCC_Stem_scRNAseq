library(Seurat)
library(ggplot2)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_5_and_S5/5CD_S5C_inferCNV_lineage"
other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"

outs_subpath <- paste0(dir,"/outs/",figure)

other_outs_subpath <- paste0(dir,"/outs/",other_figure)

day_type.integrated <- readRDS(paste0(other_outs_subpath,"/6_samples_integrated_f500_npcs50_res2.rds"))

prefix <- "S5C_inferCNV"

CNV.meta <- read.table(paste0(outs_subpath,"/",prefix,"/",prefix,".hotspot_region_mean_CNV.cutoff0.05.txt"), row.names=1)

cnv_levels <- CNV.meta$Mean_Absolute_Gain_or_Loss_of_Copy_Number
names(cnv_levels) <- rownames(x = CNV.meta)
day_type.integrated <- AddMetaData(
  object = day_type.integrated,
  metadata = cnv_levels,
  col.name = 'Mean_Absolute_Gain_or_Loss_of_Copy_Number'
)

day_type.integrated$Cluster_type<-paste(day_type.integrated$type, day_type.integrated$seurat_clusters,  sep = "_")
Idents(day_type.integrated) <- "Cluster_type"
newcluster.ids <- c("5", "2", "1", "3", "1", "1", "7", "11", "9", "6", "1", "11", "1", "4", "4", "1", "2", "1", "3", "1", "2", "2", "1", "10", "12", "1", "2", "2", "13", "2", "13", "6", "5", "1", "1", "2", "1", "2", "1", "1", "1", "2", "2", "2", "2", "1", "2", "1", "8", "2", "1", "11", "1", "1", "1", "9", "1", "2", "1", "1", "6", "6", "6", "6", "6", "6", "6", "12", "6")
names(newcluster.ids) <- levels(day_type.integrated)
day_type.integrated <- RenameIdents(day_type.integrated, newcluster.ids)
day_type.integrated$merged_clusters <- Idents(day_type.integrated)

pdf(paste0(outs_subpath,"/5C_CNV_levels_merged_1_13.pdf"),height=10,width=25)
VlnPlot(object=day_type.integrated,features="Mean_Absolute_Gain_or_Loss_of_Copy_Number", sort=TRUE, pt.size = 0.00, cols = c("#BE9C00","#E18A00","#FF65AC","#00C1AB", "#00ACFC","#00BE70","#00BBDA","#F962DD","#F8766D","#24B700", "#8CAB00","#8B93FF","#D575FE")) + theme(axis.text.x = element_text(size = 40, angle = 0, hjust=0.5), axis.text.y = element_text(size = 30), axis.title.x = element_text(size = 33), axis.title.y = element_text(size = 28), legend.position = 'none') + xlab("Cluster") + ylab("Mean Absolute Gain or Loss of Copy Number") + ggtitle(" ")
dev.off()

lineage.ids <- c("MP", "HCC", "HCC", "Endo", "Endo", "Chol", "HHP", "Endo", "MP", "HSC", "B cell", "NP", "HCC")
names(lineage.ids) <- levels(day_type.integrated)
day_type.integrated <- RenameIdents(day_type.integrated, lineage.ids)
day_type.integrated$lineage_clusters <- Idents(day_type.integrated)

pdf(paste0(outs_subpath,"/5D_CNV_levels_lineage.pdf"),height=10,width=25)
VlnPlot(object=day_type.integrated,features="Mean_Absolute_Gain_or_Loss_of_Copy_Number", sort=TRUE, pt.size = 0.00, cols = c("#8DA0CB","#FFD92F","#A6D854","#FC8D62","#B3B3B3","#E78AC3","#E5C494","#66C2A5")) + theme(axis.text.x = element_text(size = 40, angle = 0, hjust=0.5), axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 28), legend.position = 'none') + xlab(" ") + ylab("Mean Absolute Gain or Loss of Copy Number") + ggtitle(" ")
dev.off()
