rm(list = ls())
library(SingleCellExperiment)
library(ggplot2)
library(DESpace)
library(dplyr)
library(SpatialExperiment)
library(scales)
library(assertthat)
library(scater)
`%notin%` <- Negate(`%in%`)

#################### LIBD ###############################
## mixture patch
load("./Simulation/Output/151673/mixture_patch/probs_0.5_0.9FALSE/probs_0.5_0.9_final_object.rda")
load("./Simulation/Output/151673/mixture_patch/probs_0.5_0.9FALSE/result_SV_edgeR_counts.rda")
genes <- read.table("./Simulation/Output/151673/mixture_patch/probs_0.5_0.9FALSE/SV_genes.txt", header = TRUE)
results_sub <- spatial_gene_results %>% filter(adj.p.value < 0.01) 
sum(results_sub$genes %in% genes$V1)
combined_object <- logNormCounts(combined_object)

## Replace gene ids by gene names

# # Connect to ExperimentHub
# ehub <- ExperimentHub::ExperimentHub()
# # Download the full real data (about 2.1 GB in RAM) use:
# spe_all <- spatialLIBD::fetch_data(type = "spe", eh = ehub)
# df <- rowData(spe_all)[c("gene_name","gene_id")]
# rm(spe_all); rm(ehub); gc()
# ## Replace rownames of sce from gene ids to gene names
# rownames(combined_object) <- df$gene_name[match(rownames(combined_object),df$gene_id)]
# genes$V1 <- df$gene_name[match(genes$V1,df$gene_id)]
p1 <- FeaturePlot(combined_object, feature = genes$V1[c(2,482,969,1473,2477)], coordinates = c("row", "col"),
                  ncol = 5, title = FALSE)
# Break index for each sub pattern: 1,479,969,1471, 1968
## inverted mixture patch
load("./Simulation/Output/151673/mixture_reverse_patch/probs_0.5_0.9FALSE/probs_0.5_0.9_final_object.rda")
load("./Simulation/Output/151673/mixture_reverse_patch/probs_0.5_0.9FALSE/result_SV_edgeR_counts.rda")
genes <- read.table("./Simulation/Output/151673/mixture_reverse_patch/probs_0.5_0.9FALSE/SV_genes.txt", header = TRUE)
results_sub <- spatial_gene_results %>% filter(adj.p.value < 0.01) 
sum(results_sub$genes %in% genes$V1)
combined_object <- logNormCounts(combined_object)
p2 <- FeaturePlot(combined_object, feature = genes$V1[c(2,482,969,1473,2477)], coordinates = c("row", "col"),
                  ncol = 5, title = FALSE)
(AA = cowplot::plot_grid(plotlist = c(list(p1, p2)), nrow = 2))
# do.call("grid.arrange", c(list(p1,p2), nrow=2))
path_save = "./Figures/Figures/Supplementary/"
ggsave(filename = paste0(path_save,"Mix_Patterns_LIBD.pdf"),
       plot = AA,
       device = "pdf",
       width = 10,
       height = 5.5,
       units = "in",
       dpi = 300,
       limitsize = TRUE)


#################### melanoma ###############################
rm(list = ls())
library(SingleCellExperiment)
library(ggplot2)
library(DESpace)
library(dplyr)
library(SpatialExperiment)
library(scales)
library(assertthat)
library(scater)
`%notin%` <- Negate(`%in%`)
## mixture patch
load("./Simulation/Output/mel1_rep1/mixture_patch/probs_0.5_0.9FALSE/probs_0.5_0.9_final_object.rda")
load("./Simulation/Output/mel1_rep1/mixture_patch/probs_0.5_0.9FALSE/result_SV_edgeR_counts.rda")
genes <- read.table("./Simulation/Output/mel1_rep1/mixture_patch/probs_0.5_0.9FALSE/SV_genes.txt", header = TRUE)
results_sub <- spatial_gene_results %>% filter(adj.p.value < 0.01) 
sum(results_sub$genes %in% genes$V1)
combined_object <- logNormCounts(combined_object)
genes$num_pattern <- as.factor(substring(sub("^[^_]*_", "", rownames(genes)),1,1))
table(genes$num_pattern)
colData(combined_object)$row <- as.numeric(as.character(colData(combined_object)$row))
colData(combined_object)$col <- as.numeric(as.character(colData(combined_object)$col))
(p1 <- FeaturePlot(combined_object, feature = genes$V1[8], 
                  coordinates = c("row", "col"),
                  platform = "ST",
                  ncol = 1, title = FALSE) + scale_y_reverse())
(p2 <- FeaturePlot(combined_object, feature = genes$V1[1290], 
                   coordinates = c("row", "col"),
                   platform = "ST",
                   ncol = 1, title = FALSE) + scale_y_reverse())
(p3 <- FeaturePlot(combined_object, feature = genes$V1[2613], 
                   coordinates = c("row", "col"),
                   platform = "ST",
                   ncol = 1, title = FALSE) + scale_y_reverse())

# Break index for each sub pattern: 1, 1289, 2604
## inverted mixture patch
load("./Simulation/Output/mel1_rep1/mixture_reverse_patch/probs_0.5_0.9FALSE/probs_0.5_0.9_final_object.rda")
load("./Simulation/Output/mel1_rep1/mixture_reverse_patch/probs_0.5_0.9FALSE/result_SV_edgeR_counts.rda")
genes <- read.table("./Simulation/Output/mel1_rep1/mixture_reverse_patch/probs_0.5_0.9FALSE/SV_genes.txt", header = TRUE)
results_sub <- spatial_gene_results %>% filter(adj.p.value < 0.01) 
sum(results_sub$genes %in% genes$V1)
combined_object <- logNormCounts(combined_object)
colData(combined_object)$row <- as.numeric(as.character(colData(combined_object)$row))
colData(combined_object)$col <- as.numeric(as.character(colData(combined_object)$col))

(q1 <- FeaturePlot(combined_object, feature = genes$V1[8], 
                   coordinates = c("row", "col"),
                   platform = "ST",
                   ncol = 1, title = FALSE) + scale_y_reverse())
(q2 <- FeaturePlot(combined_object, feature = genes$V1[1290], 
                   coordinates = c("row", "col"),
                   platform = "ST",
                   ncol = 1, title = FALSE) + scale_y_reverse())
(q3 <- FeaturePlot(combined_object, feature = genes$V1[2613], 
                   coordinates = c("row", "col"),
                   platform = "ST",
                   ncol = 1, title = FALSE) + scale_y_reverse())

(AA = cowplot::plot_grid(plotlist = c(list(p1, p2, p3, q1, q2, q3)), nrow = 2)) 
# do.call("grid.arrange", c(list(p1,p2), nrow=2))
path_save = "./Figures/Figures/Supplementary/"
ggsave(filename = paste0(path_save,"Mix_Patterns_melanoma.pdf"),
       plot = AA,
       device = "pdf",
       width = 10,
       height = 5.5,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

#################### slideseq2 ###############################
rm(list = ls())
library(SingleCellExperiment)
library(ggplot2)
library(DESpace)
library(dplyr)
library(SpatialExperiment)
library(scales)
library(assertthat)
library(scater)
`%notin%` <- Negate(`%in%`)
low <-  "#3A3A98"
mid <-  "#F0F0F0"
high <- "#832424"
## mixture patch
load("./Simulation/Output/SlideSeq2/mixture_patch/probs_0.5_0.9FALSE/probs_0.5_0.9_final_object.rda")
load("./Simulation/Output/SlideSeq2/mixture_patch/probs_0.5_0.9FALSE/result_SV_edgeR_counts.rda")
genes <- read.table("./Simulation/Output/SlideSeq2/mixture_patch/probs_0.5_0.9FALSE/SV_genes.txt", header = TRUE)
load("./Simulation/Output/SlideSeq2/mixture_patch/probs_0.5_0.9FALSE/result_SV_edgeR_stLearncounts.rda")
results_sub <- spatial_gene_results %>% filter(adj.p.value < 0.01) 
sum(results_sub$genes %in% genes$V1)
combined_object <- logNormCounts(combined_object)
genes$num_pattern <- as.factor(substring(sub("^[^_]*_", "", rownames(genes)),1,1))
table(genes$num_pattern)

## help to find the most representative expr gene plot for each pattern
# pat2 <- genes$V1[genes$num_pattern %in% 2]
# pat4 <- genes$V1[genes$num_pattern %in% 4]
# which(spatial_gene_results$genes %in% pat4)

sub_sce <- combined_object[rownames(combined_object) == "Kcnip4"]
CD <- as.data.frame(colData(sub_sce))
CD$expr <- as.vector(logcounts(sub_sce))
(p1 <- ggplot(as.data.frame(CD), 
              aes(x=row, y=col, color=expr)) +
    geom_point(size = 0.1) + labs(color = "") + scale_y_reverse() + 
    coord_equal()+ guides(color = FALSE, size = FALSE) +
    theme_void()+ 
    labs(title = "", fill = "")+ 
    scale_color_gradient2(low=low, mid=mid,
                          high=high) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15) ))

sub_sce <- combined_object[rownames(combined_object) == spatial_gene_results$genes[98]]
CD <- as.data.frame(colData(sub_sce))
CD$expr <- as.vector(logcounts(sub_sce))
(p2 <- ggplot(as.data.frame(CD), 
              aes(x=row, y=col, color=expr)) +
    geom_point(size = 0.1) + labs(color = "") + scale_y_reverse() + 
    coord_equal()+ guides(color = FALSE, size = FALSE) +
    theme_void()+ 
    labs(title = "", fill = "")+ 
    scale_color_gradient2(low=low, mid=mid,
                          high=high) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15) ))


sub_sce <- combined_object[rownames(combined_object) == "Gdi1"]
CD <- as.data.frame(colData(sub_sce))
CD$expr <- as.vector(logcounts(sub_sce))
(p3 <- ggplot(as.data.frame(CD), 
              aes(x=row, y=col, color=expr)) +
    geom_point(size = 0.1) + labs(color = "") + scale_y_reverse() + 
    coord_equal()+ guides(color = FALSE, size = FALSE) +
    theme_void()+ 
    labs(title = "", fill = "")+ 
    scale_color_gradient2(low=low, mid=mid,
                          high=high) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15) ))


sub_sce <- combined_object[rownames(combined_object) == "Kcnd2"]
CD <- as.data.frame(colData(sub_sce))
CD$expr <- as.vector(logcounts(sub_sce))
(p4 <- ggplot(as.data.frame(CD), 
                            aes(x=row, y=col, color=expr)) +
          geom_point(size = 0.1) + labs(color = "") + scale_y_reverse() + 
          coord_equal()+ guides(color = FALSE, size = FALSE) +
          theme_void()+ 
          labs(title = "", fill = "")+ 
    scale_color_gradient2(low=low, mid=mid,
                           high=high) + 
          theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15) ))
  
# Break index for each sub pattern: 1, 1289, 2604
## inverted mixture patch
load("./Simulation/Output/SlideSeq2/mixture_reverse_patch/probs_0.5_0.9FALSE/probs_0.5_0.9_final_object.rda")
load("./Simulation/Output/SlideSeq2/mixture_reverse_patch/probs_0.5_0.9FALSE/result_SV_edgeR_counts.rda")
genes <- read.table("./Simulation/Output/Simulated/SlideSeq2/mixture_reverse_patch/probs_0.5_0.9FALSE/SV_genes.txt", header = TRUE)
results_sub <- spatial_gene_results %>% filter(adj.p.value < 0.01) 
sum(results_sub$genes %in% genes$V1)
combined_object <- logNormCounts(combined_object)
sub_sce <- combined_object[rownames(combined_object) == "Kcnip4"]
CD <- as.data.frame(colData(sub_sce))
CD$expr <- as.vector(logcounts(sub_sce))
(q1 <- ggplot(as.data.frame(CD), 
              aes(x=row, y=col, color=expr)) +
    geom_point(size = 0.1) + labs(color = "") + scale_y_reverse() + 
    coord_equal()+ guides(color = FALSE, size = FALSE) +
    theme_void()+ 
    labs(title = "", fill = "")+ 
    scale_color_gradient2(low=low, mid=mid,
                          high=high) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15) ))

sub_sce <- combined_object[rownames(combined_object) == spatial_gene_results$genes[98]]
CD <- as.data.frame(colData(sub_sce))
CD$expr <- as.vector(logcounts(sub_sce))
(q2 <- ggplot(as.data.frame(CD), 
              aes(x=row, y=col, color=expr)) +
    geom_point(size = 0.1) + labs(color = "") + scale_y_reverse() + 
    coord_equal()+ guides(color = FALSE, size = FALSE) +
    theme_void()+ 
    labs(title = "", fill = "")+ 
    scale_color_gradient2(low=low, mid=mid,
                          high=high) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15) ))


sub_sce <- combined_object[rownames(combined_object) == "Gdi1"]
CD <- as.data.frame(colData(sub_sce))
CD$expr <- as.vector(logcounts(sub_sce))
(q3 <- ggplot(as.data.frame(CD), 
              aes(x=row, y=col, color=expr)) +
    geom_point(size = 0.1) + labs(color = "") + scale_y_reverse() + 
    coord_equal()+ guides(color = FALSE, size = FALSE) +
    theme_void()+ 
    labs(title = "", fill = "")+ 
    scale_color_gradient2(low=low, mid=mid,
                          high=high) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15) ))


sub_sce <- combined_object[rownames(combined_object) == "Kcnd2"]
CD <- as.data.frame(colData(sub_sce))
CD$expr <- as.vector(logcounts(sub_sce))
(q4 <- ggplot(as.data.frame(CD), 
              aes(x=row, y=col, color=expr)) +
    geom_point(size = 0.1) + labs(color = "") + scale_y_reverse() + 
    coord_equal()+ guides(color = FALSE, size = FALSE) +
    theme_void()+ 
    labs(title = "", fill = "")+ 
    scale_color_gradient2(low=low, mid=mid,
                          high=high) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15) ))
(AA = cowplot::plot_grid(plotlist = c(list(p1, p2, p3, p4, q1, q2, q3, q4)), nrow = 2)) 
# do.call("grid.arrange", c(list(p1,p2), nrow=2))
path_save = "./Figures/Figures/Supplementary/"
ggsave(filename = paste0(path_save,"Mix_Patterns_slideseq2.pdf"),
       plot = AA,
       device = "pdf",
       width = 15,
       height = 10,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
