library(Seurat)
library(ggplot2)
library(ggsignif)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_6_and_S6/S6B_HCC_samples_CNV_levels"
outs_subpath <- paste0(dir,"/outs/",figure)

prefix <- "CNV_obs_HCC"

other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"

param <- "pc10_res0.1"

samples <- readRDS(paste0(dir,"/outs/",other_figure,"/HCC.final.SCT_",param,".rds"))

Y_axis <- "Mean_Absolute_Gain_or_Loss_of_Copy_Number"

CNV.meta <- read.table(paste0(outs_subpath,"/",prefix,"/",prefix,".hotspot_region_mean_CNV.cutoff0.05.txt"))

barcodes <- rownames(samples@meta.data)

CNV.meta <- subset(CNV.meta, rownames(CNV.meta) %in% barcodes)

CNV_column <- ""
Y_axis_name <- ""
if (Y_axis == "Mean_Absolute_Gain_or_Loss_of_Copy_Number") {
	CNV_column <- "Mean_Absolute_Gain_or_Loss_of_Copy_Number"
	Y_axis_name <- "Mean Absolute Gain or Loss of Copy Number"
} else if (Y_axis == "Number_of_high_CNV_genes") {
	CNV_column <- "Num_high_CNV_genes"
	Y_axis_name <- "Number of high CNV genes"
}

cnv_levels <- CNV.meta[[CNV_column]]
names(cnv_levels) <- rownames(x = CNV.meta)
samples <- AddMetaData(
  object = samples,
  metadata = cnv_levels,
  col.name = Y_axis
)

if (Y_axis == "Mean_Absolute_Gain_or_Loss_of_Copy_Number") {
	y_max_2 <- 0.164
} else if (Y_axis == "Number_of_high_CNV_genes") {
	y_max <- 900
}

#samples$day <- factor(samples$day, levels = c("Day3", "Day10", "Day30"))
#samples$type <- factor(samples$type, levels = c("Pos", "Neg"))

samples$day_type <- factor(samples$day_type, levels = c("Day3_Pos", "Day10_Pos", "Day30_Pos", "Day3_Neg", "Day10_Neg", "Day30_Neg"))
comparisons <- list(c("Day3_Pos", "Day10_Pos"), c("Day3_Pos", "Day30_Pos"), c("Day10_Pos", "Day30_Pos"), c("Day3_Neg", "Day10_Neg"), c("Day3_Neg", "Day30_Neg"), c("Day10_Neg", "Day30_Neg"), c("Day3_Pos", "Day3_Neg"), c("Day10_Pos", "Day10_Neg"), c("Day30_Pos", "Day30_Neg"))
Idents(samples) <- "day_type"
#for (msl in c(TRUE,FALSE)) {
for (msl in c(TRUE)) {	
	pdf(paste0(outs_subpath,"/S6B_HCC_sample_CNV_levels.map_signif_level_",msl,".pdf"),height=15,width=15)
	print(VlnPlot(object=samples,features=Y_axis, pt.size = 0.00, y.max = y_max_2, cols = c("#FB6A4A", "#FB6A4A", "#FB6A4A", "#D3D3D3", "#D3D3D3", "#D3D3D3")) + theme(axis.text.x = element_text(size = 45, angle = 45), axis.text.y = element_text(size = 45), axis.title.x = element_blank(), axis.title.y = element_text(size = 40), plot.title = element_blank(), legend.position = 'none') + ylab(Y_axis_name) + geom_signif(comparisons = comparisons, map_signif_level=msl, y_position = c(0.134, 0.146, 0.14, 0.134, 0.146, 0.14, 0.152, 0.158, 0.164)))
	#y_position = c(0.134, 0.14, 0.146)
	dev.off()
}
