library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_6_and_S6/6F_S6EFG_monocle3_pseudotime_and_gene_modules"
outs_subpath <- paste0(dir,"/outs/",figure)
dir.create(outs_subpath, recursive=TRUE)

other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
other_outs_subpath <- paste0(dir,"/outs/",other_figure)

################################### Variables ######################################################
font_size = 18

#run <- "All"

# slot <- "data"
#slot <- "scale.data"

##standard
#mode <- "standard"
#param <- "pc30_res0.1"
#HCC_integrated <- readRDS(paste0(input_dir,"/",mode,"/HCC.final.",mode,"_",param,".rds"))

# ##SCT
# mode <- "SCT"
# param <- "pc10_res0.1"
# #param <- "pc20_res0.1"

##modified by using RNA raw count
##keep SCT dimension reduction
mode <- "RNA"
param <- "pc10_res0.1"
slot <- "counts"
#param <- "pc20_res0.1"

HCC_integrated <- readRDS(paste0(other_outs_subpath,"/HCC.final.SCT_",param,".rds"))

cluster_ident <- "seurat_clusters"
####################################################################################################


###################### Use this for Pos cells only, otherwise comment out ##########################
run <- "Pos"
Idents(HCC_integrated) <- "type"
HCC_integrated <- subset(HCC_integrated, idents = run)
####################################################################################################

###################### Use this for Neg cells only, otherwise comment out ##########################
#run <- "Neg"
#Idents(HCC_integrated) <- "type"
#HCC_integrated <- subset(HCC_integrated, idents = run)
####################################################################################################


reduction <- "UMAP"
seurat_reduction <- "umap"

# assay <- ""
# if (mode == "SCT") {
# 	assay <- "SCT"
# } else if (mode == "standard") {
# 	assay <- "RNA"
# }

assay <- mode

Idents(HCC_integrated) <- cluster_ident
DefaultAssay(HCC_integrated) <- assay


monocle3_output_dir <- paste0(outs_subpath,"/monocle3_trajectory_pseudotime.",mode,"_",param,"_",slot,"_modified")
# mkdir_monocle3_output_dir <- paste0("mkdir ",monocle3_output_dir)
# system(mkdir_monocle3_output_dir)
dir.create(monocle3_output_dir)

saveRDS(HCC_integrated@meta.data, paste0(monocle3_output_dir,"/meta.",run,".rds"))

run_output_dir <- paste0(monocle3_output_dir,"/",run)
# mkdir_run_output_dir <- paste0("mkdir ", run_output_dir)
# system(mkdir_run_output_dir)
dir.create(run_output_dir)

##Following the Nature fibrotic niche paper: "We removed mitochondrial and ribosomal genes from the gene set for the purposes of trajectory analysis." 
##Omit genes with "mt-", "Rpl", "Rps", "Mrpl", "Mrps"
genenames <- rownames(HCC_integrated)
omit_genes <- grep("^Rpl", genenames, value = TRUE)
omit_genes <- c(omit_genes, grep("^Rps", genenames, value = TRUE))
omit_genes <- c(omit_genes, grep("^mt-", genenames, value = TRUE))
omit_genes <- c(omit_genes, grep("Mrpl", genenames, value = TRUE))
omit_genes <- c(omit_genes, grep("Mrps", genenames, value = TRUE))

omit_genes

omit_data <- GetAssayData(HCC_integrated, assay = assay, slot=slot)
omit_data <- omit_data[-(which(rownames(omit_data) %in% omit_genes)),]
HCC_integrated <- subset(HCC_integrated, features = rownames(omit_data))

expression_matrix <- GetAssayData(HCC_integrated, assay=assay, slot=slot)
cell_metadata <- HCC_integrated@meta.data
gene_metadata <- data.frame("gene_short_name"=as.character(rownames(HCC_integrated)), row.names = as.character(rownames(HCC_integrated)))
gene_metadata$gene_short_name <- as.character(gene_metadata$gene_short_name)

# Make CDS object
cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_metadata)

# Assign Reduction coordinate
Reduction <- HCC_integrated@reductions[[seurat_reduction]]@cell.embeddings
cds@int_colData@listData[["reducedDims"]]@listData[[reduction]] <- Reduction

# Process and align cds
# the data importantly should not be normalized again, hence norm_method = "none"
cds <- preprocess_cds(cds, num_dim = 100, norm_method = "none")

# Cluster the cells
# The cells really don't need to be clustered again at all, but running this step is required before learn_graph (which is essential) can be run
# This finds another set of cluster numbers for the data, but we will replace them back with the seurat_clusters immediately afterwards
cds <- cluster_cells(cds, reduction_method = reduction)

# Assign seurat clusters to monocle
list_cluster <- HCC_integrated@meta.data[[cluster_ident]]
names(list_cluster) <- HCC_integrated@assays[[assay]]@data@Dimnames[[2]]
cds@clusters@listData[[reduction]][["clusters"]] <- list_cluster

# Learn a graph
cds <- learn_graph(cds, use_partition = F, learn_graph_control = list(ncenter=200, minimal_branch_len = 5))

tr_plot_output_file <- paste0(run_output_dir,"/",cluster_ident,"_",reduction,"_",run,".",assay,".",slot,".trajectory.pdf")
pdf(tr_plot_output_file, width=10, height=8)
plot <- plot_cells(cds, color_cells_by = cluster_ident, label_cell_groups=F, graph_label_size = 5.5, cell_size = 0.65)
plot <- plot + FontSize(x.text = font_size, y.text = font_size, x.title = font_size, y.title = font_size, main = font_size+10) + theme(legend.position = "none") + ggtitle(run) + theme(plot.title = element_text(hjust=0.5))
print(plot)
dev.off()

# Choose the root cell(s) interactively
cds <- order_cells(cds)

pt_plot_output_file <- paste0(outs_subpath,"/S6E_",cluster_ident,"_",reduction,"_",run,".",assay,".",slot,".pseudotime.pdf")
pdf(pt_plot_output_file, width=11, height=8)
plot <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=F, graph_label_size = 5.5, cell_size = 0.65)
plot <- plot + FontSize(x.text = font_size, y.text = font_size, x.title = font_size, y.title = font_size, main = font_size+10) + ggtitle(run) + theme(legend.text=element_text(size=16), legend.title=element_text(size=16),plot.title = element_text(hjust=0.5))
print(plot)
dev.off()


saveRDS(cds, paste0(run_output_dir,"/cds_",run,".",assay,".",slot,".rds"))
#cds <- readRDS(paste0(run_output_dir,"/cds_",run,".",assay,".",slot,".rds"))


#######################################################
#Finding genes that change as a function of pseudotime
#######################################################
cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)

saveRDS(cds_pr_test_res, paste0(run_output_dir,"/cds_pr_test_res_",run,".",assay,".",slot,".rds"))
#cds_pr_test_res <- readRDS(paste0(run_output_dir,"/cds_pr_test_res_",run,".",assay,".",slot,".rds"))

# Add pseudotime bin to metadata
pt <- pseudotime(cds)
colData(cds)$pseudotime <- as.numeric(pt)

pseudotime_bin <- c()
pt_bin_number_of_cells <- list()
for (t in pt) {
  bin <- ""
  if (t >= 0 & t < 2.5) {
    bin <- "0.0-2.5"
  } else if (t >= 2.5 & t < 5.0) {
    bin <- "2.5-5.0"
  } else if (t >= 5.0 & t < 7.5) {
    bin <- "5.0-7.5"
  } else if (t >= 7.5 & t < 10.0) {
    bin <- "7.5-10.0"
  } else if (t >= 10.0 & t < 12.5) {
    bin <- "10.0-12.5"
  } else if (t >= 12.5 & t < 15.0) {
    bin <- "12.5-15.0"
  } else if (t >= 15.0 & t < 17.5) {
    bin <- "15.0-17.5"
  } else if (t >= 17.5 & t < 20.0) {
    bin <- "17.5-20.0"
  } else if (t >= 20.0 & t < 22.5) {
    bin <- "20.0-22.5"
  } else if (t >= 22.5 & t < 25.0) {
    bin <- "22.5-25.0"
  } else if (t >= 25.0 & t < 27.5) {
    bin <- "25.0-27.5"
  } else if (t >= 27.5) {
    bin <- ">= 27.5"
  }
  
  pseudotime_bin <- c(pseudotime_bin, bin)
  
  if (bin %in% names(pt_bin_number_of_cells) == "FALSE") {
    pt_bin_number_of_cells[bin] <- 0
  } else {
    pt_bin_number_of_cells[bin] <- as.numeric(pt_bin_number_of_cells[bin]) + 1
  }
}
#number of cells in each pseudotime_bin
pt_bin_number_of_cells

colData(cds)$pseudotime_bin <- pseudotime_bin

colData(cds)$pseudotime_bin <- factor(x = pseudotime_bin, levels = c("0.0-2.5", "2.5-5.0", "5.0-7.5", "7.5-10.0", "10.0-12.5", "12.5-15.0", "15.0-17.5", "17.5-20.0", "20.0-22.5","22.5-25.0", "25.0-27.5", ">= 27.5"))

# Add the pseudotime_bin to seurat object for plotting later, if necessary
HCC_integrated$pseudotime_bin <- colData(cds)$pseudotime_bin


head(colData(cds))
tail(colData(cds))


saveRDS(cds, paste0(run_output_dir,"/cds_",run,".",assay,".",slot,".rds"))
#cds <- readRDS(paste0(run_output_dir,"/cds_",run,".",assay,".",slot,".rds"))


# Check what kinds of q-values are obtained
quantile_out <- quantile(cds_pr_test_res$q_value, probs = c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1), na.rm=TRUE)
quantile_out

# Use q-value cutoff used in several publications
q_val_cutoff <- 0.05

pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < q_val_cutoff))

######## Determines the number of gene modules #############
resolution <- 1e-03
############################################################

gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=resolution)

# number of genes in each module
table(gene_module_df$module)

saveRDS(gene_module_df, paste0(run_output_dir,"/gene_module_df_",run,".",assay,".",slot,".rds"))
#gene_module_df <- readRDS(paste0(run_output_dir,"/gene_module_df_",run,".",assay,".",slot,".rds"))

# gene module UMAP plot
pdf(paste0(run_output_dir,"/Pseudotime_gene_modules_",reduction,"_",run,".",assay,".",slot,".pdf"), width=12)
plot_cells(cds, genes=gene_module_df, label_cell_groups=FALSE)
dev.off()

module_dir <- paste0(run_output_dir,"/Pseudotime_gene_modules_",run,".",assay,".",slot)
# mkdir_module_dir <- paste0("mkdir ", module_dir)
# system(mkdir_module_dir)
dir.create(module_dir)

# output module gene lists
for (i in 1:length(unique(gene_module_df$module))) {
	mod <- gene_module_df %>% filter(module %in% c(unique(gene_module_df$module)[i]))
	
	# also obtain the graph_test info of the module genes
	cds_pr_test_res_module <- cds_pr_test_res[mod$id,]
	cds_pr_test_res_module$module <- mod$module
	cds_pr_test_res_module$supermodule <- mod$supermodule
	cds_pr_test_res_module$dim_1 <- mod$dim_1
	cds_pr_test_res_module$dim_2 <- mod$dim_2
	
	write.csv(cds_pr_test_res_module,paste0(module_dir, "/Module",unique(gene_module_df$module)[i],".csv"))
}

# seurat_clusters module heatmap
#cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
#                                cell_group=colData(cds)$seurat_clusters)
#agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
#row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
#pdf(paste0(run_output_dir,"/seurat_clusters_Pseudotime_gene_modules_heatmap_",run,".",assay,".",slot,".pdf"))
#pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2", angle_col = 0)
#dev.off()

# pseudotime_bin module heatmap
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$pseudotime_bin)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

for (pheatmap_scale in c("row", "column")) {
	pdf(paste0(run_output_dir,"/pseudotime_bin_Pseudotime_gene_modules_heatmap_",run,".",assay,".",slot,"_",pheatmap_scale,".pdf"))
	pheatmap::pheatmap(agg_mat, scale=pheatmap_scale, clustering_method="ward.D2", cluster_cols = FALSE, angle_col = 45)
	dev.off()
}

####################
# Combined Modules
####################

combined_module_heatmap <- function(gene_module_df, cds, run, assay, slot, run_output_dir, pheatmap_scale) {
	# output combined module gene lists
	combined_module_dir <- paste0(run_output_dir,"/Pseudotime_gene_modules_",run,".",assay,".",slot,"_",pheatmap_scale,"_combined_module")
	# mkdir_combined_module_dir <- paste0("mkdir ", combined_module_dir)
	# system(mkdir_combined_module_dir)
	dir.create(combined_module_dir)
	
	for (i in 1:length(unique(gene_module_df$module))) {
		mod <- gene_module_df %>% filter(module %in% c(unique(gene_module_df$module)[i]))
	
		# also obtain the graph_test info of the module genes
		cds_pr_test_res_module <- cds_pr_test_res[mod$id,]
		cds_pr_test_res_module$combined_module <- mod$module
		cds_pr_test_res_module$supermodule <- mod$supermodule
		cds_pr_test_res_module$dim_1 <- mod$dim_1
		cds_pr_test_res_module$dim_2 <- mod$dim_2
	
		write.csv(cds_pr_test_res_module,paste0(combined_module_dir, "/Combined_Module",unique(gene_module_df$module)[i],".csv"))
	}
	
	cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                   cell_group=colData(cds)$pseudotime_bin)
	agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
	row.names(agg_mat) <- stringr::str_c("Combined Module ", row.names(agg_mat))
	pdf(paste0(run_output_dir,"/pseudotime_bin_Pseudotime_gene_modules_heatmap_",run,".",assay,".",slot,"_",pheatmap_scale,"_combined_module.pdf"), height=3.5)
	pheatmap::pheatmap(agg_mat, scale=pheatmap_scale, clustering_method="ward.D2", cluster_cols = FALSE, angle_col = 45)
	dev.off()
}

# #In pseudotime_bin module heatmap, scale="row" then cut the dendrogram:
# #r_combined_module1: 13
# #r_combined_module2: 9, 10, 8
# #r_combined_module3: 5, 3, 4, 1, 7
# #r_combined_module4: 11
# #r_combined_module5: 6, 2, 12
# 
# gene_module_df_row_comb <- gene_module_df
# 
# row_comb_modules <- c()
# 
# for (module in as.vector(gene_module_df_row_comb$module)) {
# 	if (module == 13) {
# 		row_comb_modules <- c(row_comb_modules, 1)
# 	} else if (module == 9 | module == 8 | module == 10) {
# 		row_comb_modules <- c(row_comb_modules, 2)
# 	} else if (module == 5 | module == 3 | module == 4 | module == 1 | module == 7) {
# 		row_comb_modules <- c(row_comb_modules, 3)
# 	} else if (module == 11) {
# 		row_comb_modules <- c(row_comb_modules, 4)
# 	} else if (module == 6 | module == 2 | module == 12) {
# 		row_comb_modules <- c(row_comb_modules, 5)
# 	} else {
# 		print("Error, module not found")
# 	}
# }

#In pseudotime_bin module heatmap, scale="row" then cut the dendrogram:
#r_combined_module3: 13
#r_combined_module2: 9, 10, 8
#r_combined_module4: 5, 3, 4, 1, 7
#r_combined_module1: 11
#r_combined_module5: 6, 2, 12


gene_module_df_row_comb <- gene_module_df

row_comb_modules <- c()

for (module in as.vector(gene_module_df_row_comb$module)) {
  if (module == 13) {
    row_comb_modules <- c(row_comb_modules, 3)
  } else if (module == 9 | module == 8 | module == 10) {
    row_comb_modules <- c(row_comb_modules, 2)
  } else if (module == 5 | module == 3 | module == 4 | module == 1 | module == 7) {
    #row_comb_modules <- c(row_comb_modules, 1)
	row_comb_modules <- c(row_comb_modules, 5)
  } else if (module == 11) {
    row_comb_modules <- c(row_comb_modules, 4)
  } else if (module == 6 | module == 2 | module == 12) {
    #row_comb_modules <- c(row_comb_modules, 5)
	row_comb_modules <- c(row_comb_modules, 1)
  } else {
    print("Error, module not found")
  }
}

gene_module_df_row_comb$module <- factor(row_comb_modules)

combined_module_heatmap(gene_module_df_row_comb, cds, run, assay, slot, run_output_dir, "row")

######
#Output files for the pseudotime heatmap
#q-value cutoff
q_val_cutoff <- 0.05
######

cds_pr_test_res_in_modules <- subset(cds_pr_test_res, rownames(cds_pr_test_res) %in% gene_module_df_row_comb$id)
cds_pr_test_res_in_modules <- cds_pr_test_res_in_modules[cds_pr_test_res_in_modules$q_value < q_val_cutoff,]

# Add the pseudotime_bin to seurat object for plotting later, if necessary
HCC_integrated$pseudotime_bin <- colData(cds)$pseudotime_bin
HCC_integrated$pseudotime <- colData(cds)$pseudotime

# Export sorted column (pseudotime) and row (combined module) metadata, and sorted expression matrix rds for input to ComplexHeatmap
pseudotime_cells <- data.frame(row.names = colnames(HCC_integrated), "pseudotime" = HCC_integrated$pseudotime)
pseudotime_cells <- pseudotime_cells[order(pseudotime_cells$pseudotime), , drop=FALSE]
write.csv(pseudotime_cells, paste0(run_output_dir,"/pseudotime_cells_",run,".",assay,".",slot,".csv"))

row_combined_module <- data.frame(row.names = gene_module_df_row_comb$id, "Combined_Module" = gene_module_df_row_comb$module)
row_combined_module <- row_combined_module[rownames(cds_pr_test_res_in_modules), , drop=FALSE]
row_combined_module <- row_combined_module[order(row_combined_module$Combined_Module), , drop=FALSE]
write.csv(row_combined_module, paste0(run_output_dir,"/row_combined_module_",run,".",assay,".",slot,".csv"))

# col_combined_module <- data.frame(row.names = gene_module_df_col_comb$id, "Combined_Module" = gene_module_df_col_comb$module)
# col_combined_module <- col_combined_module[rownames(cds_pr_test_res_in_modules), , drop=FALSE]
# col_combined_module <- col_combined_module[order(col_combined_module$Combined_Module), , drop=FALSE]
# write.csv(col_combined_module, paste0(run_output_dir,"/col_combined_module_",run,".",assay,".",slot,".csv"))

#expr_matrix <- GetAssayData(HCC_integrated, assay = assay, slot = slot)
expr_matrix = SingleCellExperiment::counts(cds)
expr_matrix@x = expr_matrix@x / rep.int(size_factors(cds), diff(expr_matrix@p))

expr_matrix <- expr_matrix[rownames(row_combined_module),]
expr_matrix <- expr_matrix[,rownames(pseudotime_cells)]

saveRDS(expr_matrix, paste0(run_output_dir,"/expr_matrix_",run,".",assay,".",slot,".rds"))

####################
# plot_genes_in_pseudotime
####################
genes <- c("Tgfb1", "Tgfb2", "Met", "Zeb1", "Ccna2", "Mki67", "Ube2c", "Cdk1", "E2f1", "Ect2", "Hells", "Top2a")
#genes <- c("Tgfb1")
genes_cds <- cds[rowData(cds)$gene_short_name %in% genes,]

plot_genes_in_pseudotime(genes_cds,
                         color_cells_by="seurat_clusters",
                         min_expr=0.5,
                         cell_size = 0.5, 
                         panel_order = genes,
                         ncol=2) +
  scale_color_manual(values=c("#377eb8", "#ff7f00",  "#4daf4a", "#e41a1c")) +
  geom_hline(yintercept=0) +
  theme(legend.title = element_blank(), 
        legend.position = "bottom",
        strip.text.x = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))
ggsave(paste0(outs_subpath,"/S6G_HCC_subcluster_along_pseudotime_selected_genes.pdf"), width = 5.7, height = 13, units="cm")
