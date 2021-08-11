##############################
#Analysis of TCGA-LIHC RNASeq#
##############################
library(DESeq2)
library(apeglm)
library(dplyr)
library(tidyverse)
library(corrplot)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_6_and_S6/6H_clinical_correlation"
outs_subpath <- paste0(dir,"/outs/",figure)
data_subpath <- paste0(dir,"/data/",figure)
dir.create(outs_subpath, recursive=TRUE)

#Read in count data from TCGA-LIHC
myfile=paste0(data_subpath,"/TCGA-LIHC_count.txt")
countsTable <- read.delim(myfile, header=TRUE, row.names=1)

#input clinical data for sample
labels <- read.table(paste0(data_subpath,"/TCGA-LIHC_clinical.txt"),header=TRUE,row.names = 1)

#Rearrange sample order
countsTable <- countsTable[, rownames(labels)]

#Remove genes without expression and patients without tumor grade information
countsTable<-countsTable[rowSums(countsTable)>0, ]
countsTable<-countsTable[,row.names(labels[!is.na(labels$Tumor_grade),])]
labels<-labels[row.names(labels[!is.na(labels$Tumor_grade),]),]

#Differential gene expression analysis of different tumor grade by DESeq2
cds <- DESeqDataSetFromMatrix(countData=countsTable, colData=labels, design = ~ Tumor_grade)
featureData <- data.frame(gene=rownames(countsTable))
mcols(cds) <- DataFrame(mcols(cds), featureData)
cds$Tumor_grade<-factor(cds$Tumor_grade, levels = c("Normal","G1","G2","G3","G4"))
cds<-DESeq(cds)

#DEG analysis
res_NTVsG1 <- results(cds,contrast = c("Tumor_grade","Normal","G1"),independentFiltering=FALSE)
res_NTVsG2 <- results(cds,contrast = c("Tumor_grade","Normal","G2"),independentFiltering=FALSE)
res_NTVsG3 <- results(cds,contrast = c("Tumor_grade","Normal","G3"),independentFiltering=FALSE)
res_NTVsG4 <- results(cds,contrast = c("Tumor_grade","Normal","G4"),independentFiltering=FALSE)
res_G1VsG2 <- results(cds,contrast = c("Tumor_grade","G1","G2"),independentFiltering=FALSE)
res_G1VsG3 <- results(cds,contrast = c("Tumor_grade","G1","G3"),independentFiltering=FALSE)
res_G1VsG4 <- results(cds,contrast = c("Tumor_grade","G1","G4"),independentFiltering=FALSE)
res_G2VsG3 <- results(cds,contrast = c("Tumor_grade","G2","G3"),independentFiltering=FALSE)
res_G2VsG4 <- results(cds,contrast = c("Tumor_grade","G2","G4"),independentFiltering=FALSE)
res_G3VsG4 <- results(cds,contrast = c("Tumor_grade","G3","G4"),independentFiltering=FALSE)

#Merge DEG analysis results
merge<-merge.data.frame(res_NTVsG1[,c(2,3,6)],res_NTVsG2[,c(2,3,6)],by.x=0, by.y=0,suffixes = c("NTvsG1","NTvsG2"))
merge<-merge.data.frame(merge,res_NTVsG3[,c(2,3,6)],by.x=1, by.y=0,suffixes = c("","NTvsG3"))
merge<-merge.data.frame(merge,res_NTVsG4[,c(2,3,6)],by.x=1, by.y=0,suffixes = c("","NTvsG4"))
merge<-merge.data.frame(merge,res_G1VsG2[,c(2,3,6)],by.x=1, by.y=0,suffixes = c("","G1VsG2"))
merge<-merge.data.frame(merge,res_G1VsG3[,c(2,3,6)],by.x=1, by.y=0,suffixes = c("","G1VsG3"))
merge<-merge.data.frame(merge,res_G1VsG4[,c(2,3,6)],by.x=1, by.y=0,suffixes = c("","G1VsG4"))
merge<-merge.data.frame(merge,res_G2VsG3[,c(2,3,6)],by.x=1, by.y=0,suffixes = c("","G2VsG3"))
merge<-merge.data.frame(merge,res_G2VsG4[,c(2,3,6)],by.x=1, by.y=0,suffixes = c("","G2VsG4"))
merge<-merge.data.frame(merge,res_G3VsG4[,c(2,3,6)],by.x=1, by.y=0,suffixes = c("","G3VsG4"))
row.names(merge)<-merge[,c(1)]
DESeq_Res_combined<-merge[,-1]
colnames(DESeq_Res_combined)[7:9]<-c("log2FoldChangeNTvsG3","lfcSENTvsG3","padjNTvsG3")

#Calculate variance-stabilized matrix of expression values
vsd <- vst(cds, blind=FALSE)

####################################
#Mapping ensemble ID to gene symbol#
####################################
library(stringi)
library(biomaRt)
human_ensembl = useEnsembl(biomart = "ensembl", dataset="hsapiens_gene_ensembl",mirror="asia")
human_ensembl_id <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), mart = human_ensembl)
vsd_merged<-merge.data.frame(human_ensembl_id,assay(vsd),by.x=1, by.y=0)
row.names(vsd_merged)<-vsd_merged$ensembl_gene_id
Res_mapped<-merge.data.frame(human_ensembl_id,DESeq_Res_combined,by.x=1, by.y=0)

#Define signatures for each tumor grade for TCGA-LIHC
Sig_NT<-Res_mapped[(Res_mapped$padjNTvsG1 < 0.05 &
                      Res_mapped$padjNTvsG2 < 0.05 &
                      Res_mapped$padjNTvsG3 < 0.05 &
                      Res_mapped$padjNTvsG4 < 0.05),]
Sig_G1<-Res_mapped[(Res_mapped$padjNTvsG1 < 0.05 &
                      Res_mapped$padjG1VsG2 < 0.05 &
                      Res_mapped$padjG1VsG3 < 0.05 &
                      Res_mapped$padjG1VsG4 < 0.05),]
Sig_G2<-Res_mapped[(Res_mapped$padjNTvsG2 < 0.05 &
                      Res_mapped$padjG1VsG2 < 0.05 &
                      Res_mapped$padjG2VsG3 < 0.05 &
                      Res_mapped$padjG2VsG4 < 0.05),]
Sig_G3<-Res_mapped[(Res_mapped$padjNTvsG3 < 0.05 &
                      Res_mapped$padjG1VsG3 < 0.05 &
                      Res_mapped$padjG2VsG3 < 0.05 &
                      Res_mapped$padjG3VsG4 < 0.05),]
Sig_G4<-Res_mapped[(Res_mapped$padjNTvsG4 < 0.05 &
                      Res_mapped$padjG1VsG4 < 0.05 &
                      Res_mapped$padjG2VsG4 < 0.05 &
                      Res_mapped$padjG3VsG4 < 0.05),]

#extract average expression of signature genes for each tumor grade
avg_NT<-rowMeans(vsd_merged[,row.names(labels[labels$Tumor_grade == "Normal",])])
avg_G1<-rowMeans(vsd_merged[,row.names(labels[labels$Tumor_grade == "G1",])])
avg_G2<-rowMeans(vsd_merged[,row.names(labels[labels$Tumor_grade == "G2",])])
avg_G3<-rowMeans(vsd_merged[,row.names(labels[labels$Tumor_grade == "G3",])])
avg_G4<-rowMeans(vsd_merged[,row.names(labels[labels$Tumor_grade == "G4",])])

avg_tumorgrade<-merge.data.frame(avg_NT,avg_G1,by.x=0,by.y=0,suffixes = c("NT","G1"))
avg_tumorgrade<-merge.data.frame(avg_tumorgrade,avg_G2,by.x=1,by.y=0,suffixes = c("","G2"))
avg_tumorgrade<-merge.data.frame(avg_tumorgrade,avg_G3,by.x=1,by.y=0,suffixes = c("","G3"))
avg_tumorgrade<-merge.data.frame(avg_tumorgrade,avg_G4,by.x=1,by.y=0,suffixes = c("","G4"))
colnames(avg_tumorgrade)<-c("EnsembleID","NT","G1","G2","G3","G4")
row.names(avg_tumorgrade)<-avg_tumorgrade[,1]
avg_tumorgrade<-avg_tumorgrade[,-1]

#calculate log fold change of signature genes for each tumor grade
logFC_NT<-rowMeans(vsd_merged[,row.names(labels[labels$Tumor_grade == "Normal",])])-rowMeans(vsd_merged[,row.names(labels[!(labels$Tumor_grade == "Normal"),])])
logFC_G1<-rowMeans(vsd_merged[,row.names(labels[labels$Tumor_grade == "G1",])])-rowMeans(vsd_merged[,row.names(labels[!(labels$Tumor_grade == "G1"),])])
logFC_G2<-rowMeans(vsd_merged[,row.names(labels[labels$Tumor_grade == "G2",])])-rowMeans(vsd_merged[,row.names(labels[!(labels$Tumor_grade == "G2"),])])
logFC_G3<-rowMeans(vsd_merged[,row.names(labels[labels$Tumor_grade == "G3",])])-rowMeans(vsd_merged[,row.names(labels[!(labels$Tumor_grade == "G3"),])])
logFC_G4<-rowMeans(vsd_merged[,row.names(labels[labels$Tumor_grade == "G4",])])-rowMeans(vsd_merged[,row.names(labels[!(labels$Tumor_grade == "G4"),])])

logFC_tumorgrade<-merge.data.frame(logFC_NT,logFC_G1,by.x=0,by.y=0,suffixes = c("NT","G1"))
logFC_tumorgrade<-merge.data.frame(logFC_tumorgrade,logFC_G2,by.x=1,by.y=0,suffixes = c("","G2"))
logFC_tumorgrade<-merge.data.frame(logFC_tumorgrade,logFC_G3,by.x=1,by.y=0,suffixes = c("","G3"))
logFC_tumorgrade<-merge.data.frame(logFC_tumorgrade,logFC_G4,by.x=1,by.y=0,suffixes = c("","G4"))
colnames(logFC_tumorgrade)<-c("EnsembleID","NT","G1","G2","G3","G4")
row.names(logFC_tumorgrade)<-logFC_tumorgrade[,1]
logFC_tumorgrade<-logFC_tumorgrade[,-1]

######################################
#read in data from single cell RNAseq#
######################################
DEGlist <- read.delim(paste0(data_subpath,"/scRNA_DEG_list.txt"), header=TRUE)
DEGlist1<- DEGlist %>% pivot_wider(id_cols=gene
                                  ,names_from=c(cluster)
                                  ,values_from = c(avg_logFC,p_val_adj)
                                  )
DEGlist1<-unnest(DEGlist1,cols = c('avg_logFC_SCT3+0', avg_logFC_SCT1, avg_logFC_SCT2, 
                                   'p_val_adj_SCT3+0', p_val_adj_SCT1, p_val_adj_SCT2))
DEGlist1[is.na(DEGlist1)]<-0

###############################################################
#convert map mouse gene ensemble ID to human human ensemble ID#
###############################################################
ensembl = useEnsembl(biomart="ensembl")
mouse_ensembl = useEnsembl(biomart = "ensembl", dataset="mmusculus_gene_ensembl")
mouse_ensembl_id <- getBM(attributes=c('ensembl_gene_id','mgi_symbol','external_gene_name','entrezgene_id','description','chromosome_name','start_position','end_position','strand','gene_biotype'), mart = mouse_ensembl) #extract id and symbol
mouse_ensembl_id_human_ortholog<-getBM(attributes=c('ensembl_gene_id','hsapiens_homolog_ensembl_gene','hsapiens_homolog_associated_gene_name'), mart = mouse_ensembl)
mouse_ensemble<-merge.data.frame(mouse_ensembl_id,mouse_ensembl_id_human_ortholog,by.x=1,by.y=1)
head(mouse_ensemble)

#Convert mouse ensemble id to human ortholog
DEGlist_merged_human<-merge.data.frame(DEGlist1,mouse_ensemble,by.x=1,by.y=3)
DEGlist_merged_human<-distinct(DEGlist_merged_human,mgi_symbol,.keep_all = TRUE)
row.names(DEGlist_merged_human)<-DEGlist_merged_human[,1]
DEGlist_merged_human<-DEGlist_merged_human[,-1]
data_corr_FC_merged<-merge.data.frame(DEGlist_merged_human,logFC_tumorgrade,by.x=16,by.y=0)
data_corr_FC_merged<-data_corr_FC_merged[,c(1:4,18:22)]
data_corr_FC_merged<-data_corr_FC_merged[!duplicated(data_corr_FC_merged$hsapiens_homolog_ensembl_gene),]
row.names(data_corr_FC_merged)<-data_corr_FC_merged[,1]
colnames(data_corr_FC_merged)[2:4]<-c("SCT3+0","SCT1","SCT2")
data_corr_FC_merged<-data_corr_FC_merged[,c(1,2,4,3,5:9)]

#remove low correlation genes
data_corr_FC_merged[(data_corr_FC_merged$'SCT3+0' < 0.2 & data_corr_FC_merged$'SCT3+0' > -0.2),c(2)]<-NA
data_corr_FC_merged[(data_corr_FC_merged$SCT2 < 0.2 & data_corr_FC_merged$SCT2 > -0.2),c(3)]<-NA
data_corr_FC_merged[(data_corr_FC_merged$SCT1 < 0.2 & data_corr_FC_merged$SCT1 > -0.2),c(4)]<-NA

#Calculate pearson correlation
corr_FC_merged_SCT3 <- cor(data_corr_FC_merged[!is.na(data_corr_FC_merged["SCT3+0"]),c(2,5:9)])
corr_FC_merged_SCT1 <- cor(data_corr_FC_merged[!is.na(data_corr_FC_merged["SCT1"]),c(4,5:9)])
corr_FC_merged_SCT2 <- cor(data_corr_FC_merged[!is.na(data_corr_FC_merged["SCT2"]),c(3,5:9)])
corr_merged<-data.matrix(data.frame(corr_FC_merged_SCT3[c(2:nrow(corr_FC_merged_SCT3)),1],
                                  corr_FC_merged_SCT2[c(2:nrow(corr_FC_merged_SCT2)),1],
                                  corr_FC_merged_SCT1[c(2:nrow(corr_FC_merged_SCT1)),1]
))
colnames(corr_merged)<-c("SCT3+0","SCT2","SCT1")

pdf(paste0(outs_subpath,"/6H_HCC_subclusters_clinical_correlation.pdf"), height=5)
corrplot(t(corr_merged)
         ,order=c("original")
         ,cl.lim = c(-1, 1)
         ,col = colorRampPalette(c("#030389","white","#8b0000"))(200)
)
dev.off()
