library(AUCell)
library(GSEABase)
library(Seurat)
library(ggpubr)
library(tidyverse)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_6_and_S6/6G_AUCell_gene_sets"
outs_subpath <- paste0(dir,"/outs/",figure)
dir.create(outs_subpath, recursive=TRUE)

data_subpath <- paste0(dir,"/data/",figure)

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

###POS ONLY####
Idents(samples) <- "type"
samples <- subset(samples, idents="Pos")
###############

Idents(samples) <- "seurat_clusters"

assay <- "RNA"
slot <- "counts"

#assay <- "SCT"
#slot <- "data"

exprMatrix <- GetAssayData(samples, assay = assay, slot = slot)

#dir <- getwd()

human_to_mouse_gene_names <- function(human_genes_list) {
  human_to_mouse_gene <- read.table(paste0(data_subpath,"/human_to_mouse_gene.txt"), header=T, sep="\t")
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
#MGI (Mouse Genome Informatics) EMT: http://www.informatics.jax.org/go/term/GO:0001837
#Can split into Pos (Up), Neg (Down)
MGI_EMT_table <- read.table(paste0(data_subpath,"/MGI_EMT/GO_term_summary_20210427_002714.txt"), sep="\t", row.names=NULL, header = TRUE)

MGI_EMT_table_up <- MGI_EMT_table[grep("positive", MGI_EMT_table$Qualifier), ]
MGI_EMT_table_down <- MGI_EMT_table[grep("negative", MGI_EMT_table$Qualifier), ]

MGI_EMT_table_up_genes <- unique(MGI_EMT_table_up$MGI.Gene.Marker.ID)
MGI_EMT_table_down_genes <- unique(MGI_EMT_table_down$MGI.Gene.Marker.ID)

MGI_EMT_table_up_geneSet <- GeneSet(MGI_EMT_table_up_genes, setName="EMT_Up")
MGI_EMT_table_down_geneSet <- GeneSet(MGI_EMT_table_down_genes, setName="EMT_Down")

#Chiang LC SUB PROL Up https://www.gsea-msigdb.org/gsea/msigdb/cards/CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_UP and Down https://www.gsea-msigdb.org/gsea/msigdb/cards/CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_DN
#CHIANG
Chiang_LC_SUB_PROL_Up_genes <- scan(paste0(data_subpath,"/CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_UP/geneset.txt"),"")
Chiang_LC_SUB_PROL_Down_genes <- scan(paste0(data_subpath,"/CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_DN/geneset.txt"),"")

Chiang_LC_SUB_PROL_Up_genes <- Chiang_LC_SUB_PROL_Up_genes[29:length(Chiang_LC_SUB_PROL_Up_genes)]
Chiang_LC_SUB_PROL_Down_genes <- Chiang_LC_SUB_PROL_Down_genes[29:length(Chiang_LC_SUB_PROL_Down_genes)]

Chiang_LC_SUB_PROL_Up_genes <- human_to_mouse_gene_names(Chiang_LC_SUB_PROL_Up_genes)
Chiang_LC_SUB_PROL_Down_genes <- human_to_mouse_gene_names(Chiang_LC_SUB_PROL_Down_genes)

Chiang_LC_SUB_PROL_Up_geneSet <- GeneSet(Chiang_LC_SUB_PROL_Up_genes, setName="Chiang_Proliferation_Up")
Chiang_LC_SUB_PROL_Down_geneSet <- GeneSet(Chiang_LC_SUB_PROL_Down_genes, setName="Chiang_Proliferation_Down")

#Liu_ES_LIKE
Liu_ES_LIKE <- c("E2f1",  "Uhrf1",  "Top2a",   "Hells",  "Ect2",   "Foxm1",  "Ccnb1",  "Ccna2",  "Pold1",  "Orc1",   "Cdc45",  "Mcm10",  "Cdc25a", "Mcm8", "Gli1", "Gli2",  "Zic2",   "Gli3",   "Hey2",   "Heyl",  "Notch3",  "Notch2",  "Rbpj",   "Notch1", "Fosl1",  "Fzd3",   "Tcf7l1", "Tcf7l2",   "Axin2",  "Myc",   "Nanog",  "Sox2", "Pou5f1", "Smad6",  "Tgfbr1")
Liu_ES_LIKE_geneSet <- GeneSet(Liu_ES_LIKE, setName="Liu_ES_Like")

##tutorial
#add two random sets and a hs-like set
set.seed(321)
#random_50_geneSet <- GeneSet(sample(rownames(exprMatrix), 100), setName="Random")
random_500_geneSet <- GeneSet(sample(rownames(exprMatrix), 500), setName="Random_500")

#countsPerGene <- apply(exprMatrix, 1, function(x) sum(x>0))
# "Housekeeping-like" - actually it is the 100 most commonly expressed genes in the dataset, so not necessarily hk genes
#hklike_geneSet <-  GeneSet(sample(names(countsPerGene)[which(countsPerGene>quantile(countsPerGene, probs=.95))], 100), setName="HK_like")

#geneSets <- GeneSetCollection(c(Wang_EMT_Up_geneSet, Wang_EMT_Down_geneSet, Hallmark_EMT_geneSet, MGI_EMT_table_up_geneSet, MGI_EMT_table_down_geneSet, Jechlinger_Mouse_EMT_Up_geneSet, Jechlinger_Mouse_EMT_Down_geneSet, random_50_geneSet, random_500_geneSet, hklike_geneSet))
###

geneSets <- GeneSetCollection(c(MGI_EMT_table_up_geneSet, 
                                MGI_EMT_table_down_geneSet, 
                                Chiang_LC_SUB_PROL_Up_geneSet, 
                                Chiang_LC_SUB_PROL_Down_geneSet,
                                Liu_ES_LIKE_geneSet,
                                random_500_geneSet))

names(geneSets)
###############################################################################################################
#AUCell pipeline: https://www.bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html
###############################################################################################################
out_dir <- paste0(outs_subpath,"/AUCell_run")
dir.create(out_dir)

set_of_cells <- paste(unique(samples$type), collapse = "_")
run_pref <- paste0(assay,"_",slot,"_",set_of_cells,"_cells")

run_dir <- paste0(out_dir,"/",run_pref)
dir.create(run_dir)

cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE)
save(cells_rankings, file=paste0(run_dir,"/",run_pref,"_cells_rankings.RData"))

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
save(cells_AUC, file=paste0(run_dir,"/",run_pref,"_cells_AUC.RData"))

#set.seed(123)
par(mfrow=c(3,6)) 
#pdf(paste0(run_dir,"/",run_pref,"_AUCell_explore_Thresholds.pdf"), width=20, height = 13)
pdf(paste0(run_dir,"/",run_pref,"_AUCell_explore_Thresholds.pdf"), width=8, height=6)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
dev.off()

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


#for (identity in c("seurat_clusters", "day", "type")) {
for (identity in c("seurat_clusters")) {
  Idents(samples) <- identity
  
  if (identity == "seurat_clusters") {
    Idents(samples) <- factor(Idents(samples), levels = c("3", "0", "2", "1"))
  }
  
  if (identity == "type") {
    
    if (length(unique(samples$type)) == 1) {
      next
    }
    pdf(paste0(run_dir,"/VlnPlot.AUCell.type.",run_pref,".pdf"), height=10, width=15)
    print(VlnPlot(samples, features = names(geneSets), pt.size=0))
    dev.off()
    
    
  } else {
    x_label_size <- 18
    for (pt_size in c(0.0001, 0)) {
      
      if (length(unique(samples$type)) == 1) {
        plots <- VlnPlot(samples, features = names(geneSets), ncol=3, pt.size=pt_size, combine=FALSE, cols = c("#e41a1c", "#377eb8", "#4daf4a", "#ff7f00"))
      } else {
        plots <- VlnPlot(samples, features = names(geneSets), ncol=3, pt.size=pt_size, split.by = "type", combine=FALSE)
      }
      
      for (i in 1:length(plots)) {
        plots[[i]] <- plots[[i]] +
          geom_boxplot(width=0.1, fill="white", outlier.shape = NA, show.legend = FALSE)+ 
          theme(legend.position = 'none', axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust=0.5, size = x_label_size)) + ylab("AUCell")
      }
      
	  out_plot_name <- paste0("VlnPlot.AUCell.",identity,".pt_size",pt_size,".",run_pref,".pdf")
	  out_plot_dir <- run_dir
	  
	  if (pt_size ==0) {
		  out_plot_name <- paste0("6G_",out_plot_name)
		  out_plot_dir <- outs_subpath
	  }
	  
      pdf(paste0(out_plot_dir,"/",out_plot_name), height=18, width=4.8)
      print(ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], ncol=2, nrow = 5))
      dev.off()
      
    }
  }
}
