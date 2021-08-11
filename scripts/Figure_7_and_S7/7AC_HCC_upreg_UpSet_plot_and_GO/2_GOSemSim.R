library(AnnotationHub)
library(GOSemSim)
library(rt)
library(ggfittext)
library(org.Mm.eg.db)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_7_and_S7/7AC_HCC_upreg_UpSet_plot_and_GO"
outs_subpath <- paste0(dir,"/outs/",figure)
data_subpath <- paste0(dir,"/data/",figure)
dir.create(outs_subpath, recursive=TRUE)

#Build annotation
hub <- AnnotationHub()
q <- query(hub, "Mus musculus")
#id <- q$ah_id[length(q)]
#Mmus <- hub[[id]]

#Prepare GO data
MmGO <- godata('org.Mm.eg.db', ont="BP")

#Read in GO list
go_list <- read.csv( file = paste0(data_subpath,"/gProfiler_results.csv"), header = TRUE)
go_bp = as.character(go_list[which(go_list$source == "GO:BP" &
                                     go_list$term_size < 1000),3])


goSim_res<-0
for (i in go_bp) {
  goSim_res<-cbind(goSim_res,mgoSim(go_bp, i, semData=MmGO, measure="Jiang", combine=NULL))
}
goSim_res<-goSim_res[,-1]
goSim_res<-as.data.frame(goSim_res)
goSim_res<-goSim_res[which(rowSums(is.na(goSim_res)) != ncol(goSim_res)),which(colSums(is.na(goSim_res)) != nrow(goSim_res))]

library("Hmisc")
library(corrplot)
library(hyperSpec)
library(psych)
library(factoextra)

palette = colorRampPalette(c("white", "red")) (100)
tiff(paste0(outs_subpath,"/corrmap_GOSemSim_tdT+_heatmap.tiff"), units="px", width=3200, height=3200, res=600) #save file for high resolution
heatmap(x = t(goSim_res), col = palette, symm = TRUE)
dev.off()

#get cluster group
fviz_nbclust(goSim_res, FUN = hcut, method = "wss")
h.row <- hclust(dist((t(goSim_res)),method="euclidean"),method="ward.D2")# for calculating distance with pearson,use pearson.dist()
plot(h.row,cex=1)
rect.hclust(h.row, h=3,border=1:20)
h.row.cut <- cutree(h.row, h=3)# use kmean will not yield the same results as clustering in heatmap
write.table(h.row.cut,paste0(outs_subpath,"/gProfiler_GOterm_Clusters.txt"),sep="\t")
