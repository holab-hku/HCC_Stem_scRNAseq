##############Install packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::version()
# If your bioconductor version is previous to 4.0, see the section bellow

## Required
BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"))
# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC")
packageVersion("SCENIC")

# # For R version 3.6 and Bioconductor 3.9:  
# devtools::install_github("aertslab/SCENIC@v1.1.2")
##############Load library
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(tidyr)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(SCENIC)
library(tidyverse)
library(patchwork)
library(arrow)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_6_and_S6/6D_motif_regulon"
outs_subpath <- paste0(dir,"/outs/",figure)
data_subpath <- paste0(dir,"/data/",figure)
dir.create(outs_subpath, recursive=TRUE)
dir.create(data_subpath, recursive=TRUE)

other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
other_outs_subpath <- paste0(dir,"/outs/",other_figure)

SCENIC_dir <- paste0(outs_subpath,"/SCENIC")
dir.create(SCENIC_dir)

param <- "pc10_res0.1"

HCC.final.SCT <- readRDS(paste0(other_outs_subpath,"/HCC.final.SCT_",param,".rds"))

exprMat <- as.matrix(HCC.final.SCT@assays$SCT@data)

##prepare cell meta data
cellInfo <- data.frame(HCC.final.SCT@meta.data)
#colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="type")] <- "celltype"
cellInfo$cluster_type<-paste0(cellInfo$cluster,"_",cellInfo$celltype)
cellInfo <- cellInfo[,c("cluster","celltype","cluster_type")]
saveRDS(cellInfo, file=paste0(SCENIC_dir,"/cellInfo.Rds"))

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=c("microglia"="forestgreen", 
                           "endothelial-mural"="darkorange", 
                           "astrocytes_ependymal"="magenta4", 
                           "oligodendrocytes"="hotpink", 
                           "interneurons"="red3", 
                           "pyramidal CA1"="skyblue", 
                           "pyramidal SS"="darkblue"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file=paste0(SCENIC_dir,"/colVars.Rds"))
#plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))


#https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/
#download:
#mm9-tss-centered-10kb-7species.mc9nr.feather
#mm9-500bp-upstream-7species.mc9nr.feather
#mm9-tss-centered-10kb-7species.mc9nr.feather.zsync
#mm9-tss-centered-10kb-7species.mc9nr.feather.gosync
#mm9-500bp-upstream-7species.mc9nr.feather.zsync
#mm9-500bp-upstream-7species.mc9nr.feather.gosync
#and place files in data_subpath/cisTarget_databases
dbDir <- paste0(data_subpath,"/cisTarget_databases") # RcisTarget databases location



#CANNOT COMPLETE
#dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
#              "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
#dir.create("cisTarget_databases")
# setwd("cisTarget_databases") # if needed
#for(featherURL in dbFiles) {
#	download.file(featherURL, destfile= paste0(dbDir, "/", basename(featherURL) )) # saved in current dir
#}




##############Initialize settings (check sha256sum before this step)

# mydbs <- c("mm9-500bp-upstream-7species.mc9nr.feather",
#            "mm9-tss-centered-10kb-7species.mc9nr.feather")
# names(mydbs) <- c("500bp", "10kb")
# scenicOptions <- initializeScenic(org="mgi", 
#                                   nCores=1,
#                                   dbDir=cisTarget_databases, 
#                                   dbs = mydbs,
#                                   datasetTitle = "HCC.final.SCT_SCENIC")

#setwd(SCENIC_dir)

org <- "mgi" # or hgnc, or dmel

myDatasetTitle <- "HCC.final.SCT_SCENIC" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 

saveRDS(scenicOptions, paste0(SCENIC_dir,"/scenicOptions.rds"))

##############Co-expression network
## Gene filter/selection
genesKept <- geneFiltering(exprMat, 
                           scenicOptions, 
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

# check whether any known relevant genes are filtered-out 
interestingGenes <- c("Nfe2l2", "Fos")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

## Correlation
runCorrelation(exprMat_filtered, scenicOptions)
#add log if necessary
#exprMat_filtered_log <- log2(exprMat_filtered+1) 

# Run GENIE3 (increase "nParts" to save memory)
runGenie3(exprMat_filtered, scenicOptions)

# Check result
#test <- readRDS("./int/1.4_GENIE3_linkList.rds")
#head(test[test$method,])

##############Build and score the GRN

#scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
#scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

#Possible failure here#
# BiocManager::install("BiocParallel")
library(BiocParallel)
library(parallel)
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 1

scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)

#Possible failure here#
# install.packages("doSNOW")
# install.packages("doParallel")
# install.packages("doMPI")
library("doSNOW")
library("doParallel")
library("doMPI")

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat)
saveRDS(scenicOptions, file=paste0(SCENIC_dir,"/scenicOptions.Rds")) # To save status

# Optional: Binarize activity
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
# savedSelections <- shiny::runApp(aucellApp)
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- paste0(SCENIC_dir,"/newThresholds.Rds")
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

# Export:
# saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
export2loom(scenicOptions, exprMat)

# To save the current status, or any changes in settings, save the object again:
saveRDS(scenicOptions, file=paste0(SCENIC_dir,"/scenicOptions.Rds") )

### Exploring output 
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
viewMotifs(tableSubset) 

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 

# Cell-type specific regulators (RSS): 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)






#Extract cell information
HCC.final.SCT$cluster_type <- paste0("C", HCC.final.SCT$seurat_clusters, "_", HCC.final.SCT$type)
Idents(HCC.final.SCT) <- "cluster_type"
cellInfo <- data.frame(seuratCluster=Idents(HCC.final.SCT))

#Extract expression matrix
expr <- GetAssayData(object = HCC.final.SCT, assay= "SCT", slot = "data") 
expr <- as(Class = 'matrix', object = expr)
