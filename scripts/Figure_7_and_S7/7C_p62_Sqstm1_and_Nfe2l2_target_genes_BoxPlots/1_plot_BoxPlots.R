library(AUCell)
library(GSEABase)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(reshape2)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_7_and_S7/7C_p62_Sqstm1_and_Nfe2l2_target_genes_BoxPlots"
outs_subpath <- paste0(dir,"/outs/",figure)
dir.create(outs_subpath, recursive=TRUE)

other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
other_outs_subpath <- paste0(dir,"/outs/",other_figure)

param <- "pc10_res0.1"
HCC.final.SCT <- readRDS(paste0(other_outs_subpath,"/HCC.final.SCT_",param,".rds"))

Keap1_Nrf3 <- c("Sqstm1", "Gpx1", "Gstm1", "Gstp1", "Gsta3", "Nqo1", "Idh1", "Ftl1", "Gclm")

gene_cell <- as.data.frame(t(as.data.frame(HCC.final.SCT@assays$SCT@counts)))
gene_cell <- log2(gene_cell+1)
gene_cell$type <- HCC.final.SCT$type
gene_cell$seurat_clusters <- HCC.final.SCT$seurat_clusters

gene_cell_stacked<-melt(gene_cell[,c(Keap1_Nrf3,"seurat_clusters","type")],id=c("seurat_clusters","type"))
colnames(gene_cell_stacked)[3]<-c("Gene")
gene_cell_stacked$seurat_clusters<-factor(gene_cell_stacked$seurat_clusters,levels=c("3","0","2","1"), labels =c("C3","C0","C2","C1"))
gene_cell_stacked$type<-factor(gene_cell_stacked$type,levels=c("Pos","Neg"), labels =c("tdT+","tdT-"))
pdf(paste0(outs_subpath,"/7C_Boxplot_Keap1_Nrf3.pdf"), width = 8, height = 4)
ggplot (gene_cell_stacked, aes(x=type, y=value, fill=seurat_clusters)) + 
  #geom_col()+
  geom_boxplot(outlier.shape = NA) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=0.02, aes(color = seurat_clusters)) +
  scale_color_manual(values=c("#e41a1c","#377eb8", "#4daf4a", "#ff7f00")) +
  scale_fill_manual(values=c("#e41a1c","#377eb8", "#4daf4a", "#ff7f00")) +
  facet_grid(row=vars(seurat_clusters),cols = vars(Gene)) +
  labs(y="Expression [Log2(count+1)]") +
  theme(plot.title = element_text(color="black", size=20, face = "plain" ),
        panel.grid =  element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        axis.title.x = element_blank(),
        #axis.line = element_blank(),
        axis.text.x  = element_text(size = 12, color = "black", face = "bold", angle = 30, vjust=1),
        #axis.ticks  = element_blank(),
        legend.position = "none",
        legend.title = element_text(size = 12,face = "bold"),
        #legend.text = element_text(color="black", size=20),
        #legend.key = element_rect(fill = "NA"),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold"),
        strip.placement = "outside")
dev.off()

#Calculate LogFC for BoxPlot
##############################
gene_cell_stacked_mean <- aggregate(gene_cell_stacked$value,list(Cluster=gene_cell_stacked$seurat_clusters,Type=gene_cell_stacked$type,Gene=gene_cell_stacked$Gene),mean) 
gene_cell_stacked_mean_cast <- dcast(gene_cell_stacked_mean,Cluster+Gene~Type, value.var = "x")
gene_cell_stacked_mean_cast$LogFC <- gene_cell_stacked_mean_cast$`tdT+` - gene_cell_stacked_mean_cast$`tdT-`

#################################show log2FC(COUNT+1) for each cluster

result <- data.frame(matrix(nrow = 1, ncol = 3))
colnames(result) <- c("Cluster", "Gene", "Log2FC")

for (i in c("0","1","2","3")){
object = subset(HCC.final.SCT, idents= i)

Idents(object)= "type"
object.pos = subset(object, idents = "Pos")
object.neg = subset(object, idents = "Neg")

for (target in c(Keap1_Nrf3)) {
  a <- mean(object.pos@assays$SCT@counts[row.names(object.pos@assays$SCT@counts) == target,])
  b <- mean(object.neg@assays$SCT@counts[row.names(object.neg@assays$SCT@counts) == target,])
  x <- log2((a+1)/(b+1))
  y <- data.frame(matrix(nrow = 1, ncol = 3))
  colnames(y) <- c("Cluster", "Gene", "Log2FC")
  y[1,1] <- i
  y[1,2] <- target
  y[1,3] <- x
  result <- rbind(result, y)
}
write.table(result,paste0(outs_subpath,"/7C_log_fold_changes.txt"), sep = "\t")}

