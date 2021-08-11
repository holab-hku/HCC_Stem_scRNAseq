#devtools::install_github("hms-dbmi/UpSetR")
library(UpSetR)
library(ComplexHeatmap)
library(Seurat)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_7_and_S7/7AC_HCC_upreg_UpSet_plot_and_GO"
outs_subpath <- paste0(dir,"/outs/",figure)
dir.create(outs_subpath, recursive=TRUE)

other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
other_outs_subpath <- paste0(dir,"/outs/",other_figure)

param <- "pc10_res0.1"
HCC.final.SCT <- readRDS(paste0(other_outs_subpath,"/HCC.final.SCT_",param,".rds"))

HCC.final.SCT$cluster_type<-paste(HCC.final.SCT$seurat_clusters, HCC.final.SCT$type, sep = "_")

Idents(HCC.final.SCT) <- "cluster_type"
DefaultAssay(HCC.final.SCT) <- "SCT"

#DEG identification
#cluster0.pos.markers.SCT <- FindMarkers(HCC.final.SCT,ident.1 = "0_Pos",ident.2 = "0_Neg", min.pct = 0.1, base )
#cluster1.pos.markers.SCT <- FindMarkers(HCC.final.SCT,ident.1 = "1_Pos",ident.2 = "1_Neg", min.pct = 0.1, base)
#cluster2.pos.markers.SCT <- FindMarkers(HCC.final.SCT,ident.1 = "2_Pos",ident.2 = "2_Neg", min.pct = 0.1, base)
#cluster3.pos.markers.SCT <- FindMarkers(HCC.final.SCT,ident.1 = "3_Pos",ident.2 = "3_Neg", min.pct = 0.1, base)

#Gene list for each cluster
#ListInput<-list(C0=row.names(cluster0.pos.markers.SCT[which(cluster0.pos.markers.SCT$p_val_adj<=0.05 & cluster0.pos.markers.SCT$avg_log2FC>0),]),
#                C1=row.names(cluster1.pos.markers.SCT[which(cluster1.pos.markers.SCT$p_val_adj<=0.05 & cluster1.pos.markers.SCT$avg_log2FC>0),]),
#                C2=row.names(cluster2.pos.markers.SCT[which(cluster2.pos.markers.SCT$p_val_adj<=0.05 & cluster2.pos.markers.SCT$avg_log2FC>0),]),
#                C3=row.names(cluster3.pos.markers.SCT[which(cluster3.pos.markers.SCT$p_val_adj<=0.05 & cluster3.pos.markers.SCT$avg_log2FC>0),]))

cluster0.pos.markers.SCT <- FindMarkers(HCC.final.SCT,ident.1 = "0_Pos",ident.2 = "0_Neg", min.pct = 0.1, base = 2.718281828459)
cluster1.pos.markers.SCT <- FindMarkers(HCC.final.SCT,ident.1 = "1_Pos",ident.2 = "1_Neg", min.pct = 0.1, base = 2.718281828459)
cluster2.pos.markers.SCT <- FindMarkers(HCC.final.SCT,ident.1 = "2_Pos",ident.2 = "2_Neg", min.pct = 0.1, base = 2.718281828459)
cluster3.pos.markers.SCT <- FindMarkers(HCC.final.SCT,ident.1 = "3_Pos",ident.2 = "3_Neg", min.pct = 0.1, base = 2.718281828459)

#Gene list for each cluster
ListInput<-list(C0=row.names(cluster0.pos.markers.SCT[which(cluster0.pos.markers.SCT$p_val_adj<=0.05 & cluster0.pos.markers.SCT$avg_log2.718281828459FC>0),]),
                C1=row.names(cluster1.pos.markers.SCT[which(cluster1.pos.markers.SCT$p_val_adj<=0.05 & cluster1.pos.markers.SCT$avg_log2.718281828459FC>0),]),
                C2=row.names(cluster2.pos.markers.SCT[which(cluster2.pos.markers.SCT$p_val_adj<=0.05 & cluster2.pos.markers.SCT$avg_log2.718281828459FC>0),]),
                C3=row.names(cluster3.pos.markers.SCT[which(cluster3.pos.markers.SCT$p_val_adj<=0.05 & cluster3.pos.markers.SCT$avg_log2.718281828459FC>0),]))


#To check specific intersection set
m = make_comb_mat(ListInput)
conserved_genes<-extract_comb(m, "1111")

#upset plot
pdf(paste0(outs_subpath,"/HCC_subclusters_UpSet_plot.pdf"),height=4)
upset(fromList(ListInput),
                  sets = c("C1","C2","C0","C3"),
                  keep.order = TRUE,
                  order.by = "freq",
                  empty.intersections = "on",
                  sets.x.label = "Upregulated genes" 
)
dev.off()
