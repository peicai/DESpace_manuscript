rm(list = ls())
library(SingleCellExperiment)
library(ggplot2)
library(DESpace)
library(dplyr)
library(SpatialExperiment)
library(scales)
library(assertthat)
`%notin%` <- Negate(`%in%`)
## bottom pattern
load("./Simulation/Output/151673/bottom_patch/probs_0.5_0.9FALSE/probs_0.5_0.9_final_object.rda")
load("./Simulation/Output/151673/bottom_patch/probs_0.5_0.9FALSE/result_SV_edgeR_counts.rda")
genes <- read.table("./Simulation/Output/151673/bottom_patch/probs_0.5_0.9FALSE/SV_genes.txt", header = TRUE)
results_sub <- spatial_gene_results %>% filter(adj.p.value < 0.01) 
sum(results_sub$genes %in% genes$gene_ids)
FeaturePlot(final_object, feature = results_sub$genes[1:10], coordinates = c("row", "col"),
            ncol = 5, title = TRUE)


## bottom pattern
load("./Simulation/Output/151673/bottom_patch/probs_0.5_0.9FALSE/probs_0.9_final_object.rda")
load("./Simulation/Output/151673/bottom_patch/probs_0.5_0.9FALSE/result_SV_edgeR_counts.rda")
genes <- read.table("./Simulation/Output/151673/bottom_patch/probs_0.5_0.9FALSE/SV_genes.txt", header = TRUE)
results_sub <- spatial_gene_results %>% filter(adj.p.value < 0.01) 
sum(results_sub$genes %in% genes$gene_ids)
gg_A = gg
results_sub_A = results_sub
(LIBD_A <- FeaturePlot(gg_A, feature = results_sub_A$genes[1], coordinates = c("row", "col"),
            ncol = 1, title = FALSE))

## circular pattern
load("./Simulation/Output/151673/circle_patch/probs_0.5_0.9FALSE/probs_0.9_final_object.rda")
load("./Simulation/Output/151673/circle_patch/probs_0.5_0.9FALSE/result_SV_edgeR_counts.rda")
genes <- read.table("./Simulation/Output/151673/circle_patch/probs_0.5_0.9FALSE/SV_genes.txt", header = TRUE)
results_sub <- spatial_gene_results %>% filter(adj.p.value < 0.01) 
sum(results_sub$genes %in% genes$gene_ids)
gg_B = gg
results_sub_B = results_sub
(LIBD_B <- FeaturePlot(gg_B, feature = results_sub_B$genes[1], coordinates = c("row", "col"),
                       ncol = 1, title = FALSE))
## Annotation pattern
load("./Simulation/Output/151673/Manual_clusters_patch/probs_0.5_0.9FALSE/probs_0.9_final_object.rda")
load("./Simulation/Output/151673/Manual_clusters_patch/probs_0.5_0.9FALSE/result_SV_edgeR_counts.rda")
genes <- read.table("./Simulation/Output/151673/Manual_clusters_patch/probs_0.5_0.9FALSE/SV_genes.txt", header = TRUE)
results_sub <- spatial_gene_results %>% filter(adj.p.value < 0.01) 
sum(results_sub$genes %in% genes$gene_ids)
gg_C = gg
results_sub_C = results_sub
(LIBD_C <- FeaturePlot(gg_C, feature = results_sub_C$genes[1], coordinates = c("row", "col"),
                       ncol = 1, title = FALSE))

(AA = egg::ggarrange( plots = list(LIBD_A + labs(title = "Bottom/Right") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)),
                                   LIBD_B + labs(title = "Circular")+ theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)),
                                   LIBD_C + labs(title = "Annotations")+ theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12))),ncol = 3, nrow = 1))
path_save = "./Figures/Figures/Simulated/"
png(filename = paste0(path_save,"Expression_plot.png"), 
    width = 10, height = 4.5, units = "in",
    pointsize = 12, bg = "white", res = 300)   

AA

dev.off()
