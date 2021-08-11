library(infercnv)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_5_and_S5/5CD_S5C_inferCNV_lineage"

outs_subpath <- paste0(dir,"/outs/",figure)

input_dir <- paste0(outs_subpath,"/lineage_counts_matrix_and_meta")

counts_matrix <- paste0(input_dir,"/exprMatrix.genes_out.tsv")
annot_file <- paste0(input_dir,"/meta.noheader.tsv")
gene_ordering <- paste0(input_dir,"/gene_ordering.tsv")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=annot_file,
                                    delim="\t",
                                    gene_order_file=gene_ordering,
									ref_group_names=c("5","4"))

									

infercnv_output_dir <- paste0(outs_subpath,"/S5C_inferCNV")
dir.create(infercnv_output_dir)

threads <- 12

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=infercnv_output_dir,  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             num_threads = threads,
							 denoise=T,
							 output_format="pdf",
                             HMM=F
							 )