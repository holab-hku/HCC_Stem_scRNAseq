library(Seurat)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_5_and_S5/5CD_S5C_inferCNV_lineage"
other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"

data_subpath <- paste0(dir,"/data/",figure)
outs_subpath <- paste0(dir,"/outs/",figure)

dir.create(outs_subpath)
counts_and_meta_dir <- paste0(outs_subpath,"/lineage_counts_matrix_and_meta")
dir.create(counts_and_meta_dir)

other_data_subpath <- paste0(dir,"/data/",other_figure)
other_outs_subpath <- paste0(dir,"/outs/",other_figure)

day_type.integrated <- readRDS(paste0(other_outs_subpath,"/6_samples_integrated_f500_npcs50_res2.rds"))

day_type.integrated$Cluster_type<-paste(day_type.integrated$type, day_type.integrated$seurat_clusters,  sep = "_")
Idents(day_type.integrated) <- "Cluster_type"

newcluster.ids <- c("5", "2", "1", "3", "1", "1", "7", "11", "9", "6", "1", "11", "1", "4", "4", "1", "2", "1", "3", "1", "2", "2", "1", "10", "12", "1", "2", "2", "13", "2", "13", "6", "5", "1", "1", "2", "1", "2", "1", "1", "1", "2", "2", "2", "2", "1", "2", "1", "8", "2", "1", "11", "1", "1", "1", "9", "1", "2", "1", "1", "10", "6", "6", "6", "6", "6", "6", "12", "6")
names(newcluster.ids) <- levels(day_type.integrated)
day_type.integrated <- RenameIdents(day_type.integrated, newcluster.ids)
day_type.integrated$merged_clusters <- Idents(day_type.integrated)

day_type.integrated.meta <- day_type.integrated@meta.data
day_type.integrated.meta <- day_type.integrated.meta[c("nCount_RNA", "nFeature_RNA","merged_clusters")]
write.table(day_type.integrated.meta, paste0(counts_and_meta_dir,"/meta.tsv"), sep = "\t", col.names = NA, quote=F)

raw_data_all <- GetAssayData(day_type.integrated, slot = "counts", assay = "RNA")
raw_data_all_df <- tibble::rownames_to_column(as.data.frame(raw_data_all), "gene")
write.table(raw_data_all_df, paste0(counts_and_meta_dir,"/exprMatrix.tsv"), row.names=FALSE, sep = "\t", quote=F)