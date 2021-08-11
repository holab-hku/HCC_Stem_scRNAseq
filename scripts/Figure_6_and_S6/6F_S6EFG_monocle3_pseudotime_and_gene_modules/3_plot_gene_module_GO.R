library(ggplot2)
library(tidyverse)
library(reshape2)
library(cowplot)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_6_and_S6/6F_S6EFG_monocle3_pseudotime_and_gene_modules"
outs_subpath <- paste0(dir,"/outs/",figure)

run <- "Pos"
slot <- "counts"
mode <- "RNA"
param <- "pc10_res0.1"

monocle3_output_dir <-  paste0(outs_subpath,"/monocle3_trajectory_pseudotime.",mode,"_",param,"_",slot,"_modified")
run_output_dir <- paste0(monocle3_output_dir,"/",run)


#details of GO_to_plot.txt
#
data_plot<- read.table(paste0(run_output_dir,"/GO_to_plot.txt"), header=TRUE,sep="\t")
# modules_palette <- hcl.colors(length(unique(data_plot$Module)), palette = 'Zissou1')
#modules_palette_new <- c("#e41a1c", "#4daf4a", "#ff7f00", "#984ea3", "#377eb8")
modules_palette_new <- c("#377eb8", "#4daf4a", "#ff7f00", "#e41a1c", "#984ea3")

data_plot$term <- paste0(data_plot$source, "_", data_plot$term_name)

plots <- list()

for (i in unique(data_plot$Module)) {

	barchart<-data_plot[data_plot$Module == i & data_plot$Print == "Y",]

	y<-barchart$term #assign value for y axis
	x<- -log10(barchart$adjusted_p_value) #assign value for x asix

	pdf(paste0(outs_subpath,"/6F_S6F_Module_", i, "_GO.pdf"), height=4, width=8)
	print(ggplot(barchart,aes(x,y=reorder(y,-adjusted_p_value)))+
	  coord_cartesian(xlim =c(0, as.integer(max(x)+1)), ylim = c(1, 5))+
	  geom_col(fill=modules_palette_new[i], alpha = 0.5)+
	  scale_x_continuous(breaks = seq(0, as.integer(max(x)+1.5), by = as.integer((max(x)+1)/6)))+
	  labs(x = "-log10(p-value)",
	       title = paste0("Module ", i))+
	  theme(
	    plot.title = element_text(color="black", size=18, face="bold"),
	    axis.title.x = element_text(color="black", size=18, face="bold"),
	    axis.title.y =  element_blank(),
	    axis.text.y = element_blank(),
	    axis.text.x = element_text(color="black",size=18),
	    legend.title = element_text(size = 16,face = "bold"),
	    legend.text = element_text(color="black", size=6),
	    legend.key = element_rect(fill = "NA"),
	    plot.background = element_rect(color = "white"),
	    panel.background = element_rect(fill = NA),
	    plot.margin = unit(c(10, 10, 10, 10), "pt"),
	    axis.line.x = element_line(colour = "black"),
	    axis.line.y = element_line(colour = "black"),
	  )+
	  annotate(geom="text",size=6,x=0.1,y=y,label=paste0(y, " (", barchart$intersection_size, ")"), hjust=0))
	dev.off()
}
