ibrary(Seurat)
library(loomR)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_6_and_S6/6E_S6D_scVelo_RNA_velocity"

data_subpath <- paste0(dir,"/data/",figure)

other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
other_outs_subpath <- paste0(dir,"/outs/",other_figure)

param <- "pc10_res0.1"
HCC_integrated <- readRDS(paste0(other_outs_subpath,"/HCC.final.",mode,"_",param,".rds"))

Idents(HCC_integrated) <- "type"
HCC_integrated.pos <- subset(HCC_integrated, idents="Pos")

saveRDS(HCC_integrated.pos, paste0(output_dir,"/HCC.final.",mode,"_",param,".Pos.rds"))
HCC_integrated.pos.loom <- as.loom(HCC_integrated.pos, filename =  paste0(output_dir,"/HCC.final.",mode,"_",param,".Pos.loom"), verbose = FALSE)
HCC_integrated.pos.loom
HCC_integrated.pos.loom$close_all() 

HCC_integrated.neg <- subset(HCC_integrated, idents="Neg")

saveRDS(HCC_integrated.neg, paste0(output_dir,"/HCC.final.",mode,"_",param,".Neg.rds"))
HCC_integrated.neg.loom <- as.loom(HCC_integrated.neg, filename =  paste0(output_dir,"/HCC.final.",mode,"_",param,".Neg.loom"), verbose = FALSE)
HCC_integrated.neg.loom
HCC_integrated.neg.loom$close_all() 
