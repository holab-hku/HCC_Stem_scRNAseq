library(Seurat)
dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_6_and_S6/S6B_HCC_samples_CNV_levels"
outs_subpath <- paste0(dir,"/outs/",figure)
dir.create(outs_subpath, recursive=TRUE)

other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
other_figure_2 <- "Figure_5_and_S5/5CD_S5C_inferCNV_lineage/S5C_inferCNV"

param <- "pc10_res0.1"

samples <- readRDS(paste0(dir,"/outs/",other_figure,"/HCC.final.SCT_",param,".rds"))
samples.bcs <- colnames(samples)

CNV_obs_mat <- read.table(paste0(dir,"/outs/",other_figure_2,"/infercnv.observations.txt"),header=T,check.names=F)

CNV_obs_mat <- CNV_obs_mat[,samples.bcs]

CNV_obs_HCC <- paste0(outs_subpath,"/CNV_obs_HCC")
dir.create(CNV_obs_HCC)

write.table(CNV_obs_mat, paste0(CNV_obs_HCC,"/infercnv.observations.txt"))