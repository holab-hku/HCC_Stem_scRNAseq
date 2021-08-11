library(ComplexHeatmap)
library(dendsort)
library(dendextend)
library(circlize)
library(Matrix)
library(Hmisc)
library(dplyr)
library(EnrichedHeatmap)
library(ggpubr)
library(tidyverse)
library(cowplot)

dir <- "HCC_Stem_scRNAseq"
##########################

figure <- "Figure_6_and_S6/6F_S6EFG_monocle3_pseudotime_and_gene_modules"
outs_subpath <- paste0(dir,"/outs/",figure)
dir.create(outs_subpath, recursive=TRUE)

#other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
#other_outs_subpath <- paste0(dir,"/outs/",other_figure)

################################### Variables ######################################################

combined_module_scale <- "row"
#combined_module_scale <- "col"

#run <- "All"
run <- "Pos"
#run <- "Neg"

slot <- "counts"
#slot <- "scale.data"

##standard
mode <- "RNA"
#param <- "pc30_res0.1"

##SCT
# mode <- "SCT"
param <- "pc10_res0.1"
#param <- "pc20_res0.1"

assay <- ""
if (mode == "SCT") {
  assay <- "SCT"
} else if (mode == "RNA") {
  assay <- "RNA"
}

monocle3_output_dir <-  paste0(outs_subpath,"/monocle3_trajectory_pseudotime.",mode,"_",param,"_",slot,"_modified")
run_output_dir <- paste0(monocle3_output_dir,"/",run)

####################################################################################################

pseudotime_cells <- ""
if (mode == "SCT") {
  pseudotime_cells <- read.csv(paste0(run_output_dir,"/pseudotime_cells_",run,".SCT",".",slot,".csv"), header=T, row.names=1)
} else {
  pseudotime_cells <- read.csv(paste0(run_output_dir, "/pseudotime_cells_Pos.RNA.counts.csv"), header=T, row.names=1)
}

#meta <- readRDS(paste0(monocle3_output_dir,"/meta.",run,".rds"))

#meta <- meta[rownames(pseudotime_cells),]

expr_matrix <- readRDS(paste0(run_output_dir,"/expr_matrix_",run,".",assay,".",slot,".rds"))
expr_matrix <- as.matrix(expr_matrix)

combined_module <- ""

combined_module <- read.csv(paste0(run_output_dir,"/",combined_module_scale,"_combined_module_",run,".",assay,".",slot,".csv"), header=T, row.names=1)
expr_matrix <- expr_matrix[rownames(combined_module),]

#################################
#Module Mean Expression Curves function
#################################
x_y_axis_text_size = 17
line_size = 1.5
point_size = 0.2
modules <- rownames(table(combined_module$Combined_Module))
y_axis_name <- "Expression (z-scored)"

#modules_palette <- hcl.colors(length(modules), palette = 'Zissou1')
#modules_palette <- c("#e41a1c", "#4daf4a", "#ff7f00", "#984ea3", "#377eb8")
modules_palette <- c("#377eb8", "#4daf4a", "#ff7f00", "#e41a1c", "#984ea3")

#out_file_dir <- paste0(run_output_dir,"/Expression_curve")
#dir.create(out_file_dir)

plots <- list()

for (i in as.integer(modules)) {
combined_module_specific <- combined_module[combined_module$Combined_Module == i, , drop = FALSE]
expr_matrix_specific <- expr_matrix[rownames(combined_module_specific),]
    
combined_module_mean_expr <- colMeans(t(scale(t(expr_matrix_specific))))
  
df <- data.frame("expr" = as.vector(combined_module_mean_expr), 
                 "pseudotime" = pseudotime_cells$pseudotime, 
                 "set" = str_sub(names(combined_module_mean_expr),-3,-1))
  
theme_set(theme_bw())
  

plots[[i]] <- ggplot(df, aes(x = pseudotime, y = expr, color=set)) + 
      geom_smooth(aes(fill=set), size=line_size, show.legend = FALSE, method = "loess") + #"gam", "loess"
      ylab(paste0(y_axis_name)) + 
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"), 
            axis.text = element_text(size = x_y_axis_text_size), 
            axis.title = element_text(size = x_y_axis_text_size), 
            plot.title = element_text(color = modules_palette[i:i], size = x_y_axis_text_size+3, hjust=0.5), 
            legend.title = element_blank(), 
            legend.text=element_blank()) + 
            xlim(0, 30)+
            ggtitle(paste0("Combined Module ",i)) + 
            scale_color_manual(values=paste0(modules_palette[i])) + 
            scale_fill_manual(values=paste0(modules_palette[i]))
}


pdf(paste0(outs_subpath, "/6F_S6F_gene_module_expr_curves.pdf"), height=20, width=5)
print(plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], ncol=1, align = "v"))
dev.off()
