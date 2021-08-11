library(AUCell)
library(GSEABase)
library(Seurat)
library(ggpubr)
library(tidyverse)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_7_and_S7/S7D_MTORC1_and_NFE2L2.V2_gene_sets_AUCell"
outs_subpath <- paste0(dir,"/outs/",figure)
dir.create(outs_subpath, recursive=TRUE)

data_subpath <- paste0(dir,"/data/",figure)
data_subpath_0 <- paste0(dir,"/data/Figure_7_and_S7")

other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
other_outs_subpath <- paste0(dir,"/outs/",other_figure)

param <- "pc10_res0.1"
samples <- readRDS(paste0(other_outs_subpath,"/HCC.final.SCT_",param,".rds"))

Idents(samples) <- "seurat_clusters"
Idents(samples) <- factor(Idents(samples), levels = c("3", "0", "2", "1"))

Idents(samples) <- "day_type"
Idents(samples) <- factor(Idents(samples), levels = c("Day3_Neg", "Day3_Pos", "Day10_Neg", "Day10_Pos", "Day30_Neg", "Day30_Pos"))

Idents(samples) <- "day"
Idents(samples) <- factor(Idents(samples), levels = c("Day3", "Day10", "Day30"))

Idents(samples) <- "type"
Idents(samples) <- factor(Idents(samples), levels = c("Pos", "Neg"))

# ###POS ONLY####
# Idents(samples) <- "type"
# samples <- subset(samples, idents="Pos")
# ###############

# Idents(samples) <- "seurat_clusters"
Idents(samples) <- "type"

assay <- "RNA"
slot <- "counts"

#assay <- "SCT"
#slot <- "data"

exprMatrix <- GetAssayData(samples, assay = assay, slot = slot)

human_to_mouse_gene_names <- function(human_genes_list) {
  human_to_mouse_gene <- read.table(paste0(data_subpath_0,"/human_to_mouse_gene.txt"), header=T, sep="\t")
  rownames(human_to_mouse_gene) <- human_to_mouse_gene$human_gene
  
  intersected_genes <- intersect(human_genes_list, rownames(human_to_mouse_gene))
  
  human_to_mouse_gene_intersect <- human_to_mouse_gene[intersected_genes,]
  final_list <- human_to_mouse_gene_intersect$mouse_gene
  
  no_intersect <- setdiff(human_genes_list, intersected_genes)
  if (length(no_intersect) > 0) {
    cat("The following genes can't be found in the MGI table, will convert them anyway to the lowercase format:\n")
    print(no_intersect)
    
    no_intersect_lc <- str_to_sentence(no_intersect)
    print(no_intersect_lc)
    
    final_list <- c(final_list, no_intersect_lc)
  }
  return(final_list)
  
}

############################
#create gene set collection#
############################

#1. HALLMARK_MTORC1_SIGNALING
Hallmark_MTORC1_genes <- scan(paste0(data_subpath,"/HALLMARK_MTORC1_SIGNALING.txt"),"")
Hallmark_MTORC1_genes
Hallmark_MTORC1_genes <- Hallmark_MTORC1_genes[10:length(Hallmark_MTORC1_genes)]
Hallmark_MTORC1_genes <- human_to_mouse_gene_names(Hallmark_MTORC1_genes)
Hallmark_MTORC1_geneSet <- GeneSet(Hallmark_MTORC1_genes, setName="Hallmark_MTORC1")

#2. NFE2L2.V2_SIGNALING
NFE2L2.V2_genes <- scan(paste0(data_subpath,"/NFE2L2.V2.txt"),"")
NFE2L2.V2_genes
NFE2L2.V2_genes <- NFE2L2.V2_genes[17:length(NFE2L2.V2_genes)]
NFE2L2.V2_genes <- human_to_mouse_gene_names(NFE2L2.V2_genes)
NFE2L2.V2_geneSet <- GeneSet(unique(NFE2L2.V2_genes), setName="NFE2L2.V2_genes")

##tutorial
#add two random sets and a hs-like set
#set.seed(321)
#random_50_geneSet <- GeneSet(sample(rownames(exprMatrix), 100), setName="Random")
#random_500_geneSet <- GeneSet(sample(rownames(exprMatrix), 500), setName="Random_2")

#countsPerGene <- apply(exprMatrix, 1, function(x) sum(x>0))
# "Housekeeping-like" - actually it is the 100 most commonly expressed genes in the dataset, so not necessarily hk genes
#hklike_geneSet <-  GeneSet(sample(names(countsPerGene)[which(countsPerGene>quantile(countsPerGene, probs=.95))], 100), setName="HK_like")

#geneSets <- GeneSetCollection(c(Wang_EMT_Up_geneSet, Wang_EMT_Down_geneSet, Hallmark_EMT_geneSet, MGI_EMT_table_up_geneSet, MGI_EMT_table_down_geneSet, Jechlinger_Mouse_EMT_Up_geneSet, Jechlinger_Mouse_EMT_Down_geneSet, random_50_geneSet, random_500_geneSet, hklike_geneSet))
###

geneSets <- GeneSetCollection(c(Hallmark_MTORC1_geneSet, NFE2L2.V2_geneSet))

names(geneSets)
###############################################################################################################
#AUCell pipeline: https://www.bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html
###############################################################################################################
#out_dir <- paste0(getwd(),"/outs_paper")
#dir.create(out_dir)

set_of_cells <- paste(unique(samples$type), collapse = "_")
run_pref <- paste0(assay,"_",slot,"_",set_of_cells,"_cells")

run_dir <- paste0(outs_subpath,"/AUCell_run")
dir.create(run_dir,recursive=TRUE)

cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE)
save(cells_rankings, file=paste0(run_dir,"/",run_pref,"_cells_rankings.RData"))

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
save(cells_AUC, file=paste0(run_dir,"/",run_pref,"_cells_AUC.RData"))

# #set.seed(123)
# par(mfrow=c(3,6)) 
# #pdf(paste0(run_dir,"/",run_pref,"_AUCell_explore_Thresholds.pdf"), width=20, height = 13)
# pdf(paste0(run_dir,"/",run_pref,"_AUCell_explore_Thresholds.pdf"), width=8, height=6)
# cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
# dev.off()

##Make the matrix of cell x gene_set AUC values
out_AUC <- as.data.frame(matrix(nrow = length(cells_AUC@colData@rownames)))
for (i in 1:length(cells_AUC@NAMES)) {
  out_AUC[[cells_AUC@NAMES[i]]] <- as.vector(cells_AUC@assays@data@listData$AUC[i,])
}
out_AUC$V1 <- NULL
rownames(out_AUC) <- cells_AUC@colData@rownames

##Add AUC back to seurat metadata
for (i in 1:dim(out_AUC)[2]) {
  to_add <- out_AUC[[i]]
  names(to_add) <- rownames(x=out_AUC)
  
  samples <- AddMetaData(
    object = samples,
    metadata = to_add,
    col.name = colnames(out_AUC)[i]
  )
}

##plot better graphs
Idents(samples) <- "type"
Idents(samples) <- factor(Idents(samples), levels = c("Pos", "Neg"))

Idents(samples) <- "seurat_clusters"
Idents(samples) <- factor(Idents(samples), levels = c("3", "0", "2", "1"))

for (i in names(geneSets)){
png(paste0(outs_subpath,"/S7D_", i ,"_cluster_type.png"), units = "cm", res=600, width = 12, height = 6)
print(VlnPlot(samples, features = i, pt.size=0, split.by = "type", cols = c("#bebebe", "#fb6a4a"))+
  # geom_boxplot(width=0.1, fill="white", outlier.shape = NA, show.legend = FALSE)+
  labs(x="Cluster", y = "AUCell", title = i) +
  scale_x_discrete(labels=c("3" = "C3", "0" = "C0",
                            "2" = "C2", "1" = "C1")) +
  scale_color_manual(labels = c("Neg"="tdT-", "Pos"="tdT+"), values = c("#bebebe", "#fb6a4a"))+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA, show.legend = FALSE)+ 
  theme(axis.text.x = element_text(face="plain", color="black", 
                                   size=12, angle=0, hjust = 0.5),
        axis.text.y = element_text(face="plain", color="black", 
                                   size=12, angle=0),
        axis.title.x = element_text(face="bold", color="black", 
                                    size=12, angle=0),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=12, angle=90),
        plot.title = element_text(color="black", size=14, face = "bold"), 
        legend.position = "right"))
dev.off()
}

###############################################################################################################
#calculate p value
###############################################################################################################
calc_within_subcluster_type_p_value <- function(samples, identity_class, column_name) {
  Idents(samples) <- identity_class
  
  comparisons <- c()
  p_values <- c()
  
  for (subcluster in as.numeric(names(table(Idents(samples))))) {
    samples.c <- subset(samples, idents=subcluster)
    Idents(samples.c) <- "type"
    samples.c.Pos <- subset(samples.c, idents="Pos")
    samples.c.Neg <- subset(samples.c, idents="Neg")
    
    test <- wilcox.test(samples.c.Pos[[column_name]][,1], samples.c.Neg[[column_name]][,1])
    
    comp <- paste0(subcluster,"_Pos vs ",subcluster,"_Neg")
    cat(comp,"\n")
    
    comparisons <- c(comparisons, comp)
    p_values <- c(p_values, test$p.value)
    
  }
  
  out_df <- data.frame("col1" = comparisons, "p.value" = p_values)
  names(out_df)[1] <- column_name
  
  print(out_df)
  
  write.csv(out_df, paste0(run_dir,"/",column_name,".within_subcluster_type.p_values.csv"), row.names=F, quote=F)	
}

calc_within_subcluster_type_p_value(samples, "seurat_clusters", "Hallmark_MTORC1")
calc_within_subcluster_type_p_value(samples, "seurat_clusters", "NFE2L2.V2_genes")
