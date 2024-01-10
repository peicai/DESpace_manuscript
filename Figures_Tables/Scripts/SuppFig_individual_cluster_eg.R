rm(list = ls())

load("./Data/LIBD/151673_sce_LIBD.rda")
library(DESpace)
library(ggplot2)
library(ggforce)
library(patchwork)
`%notin%` <- Negate(`%in%`)
# Specify column names of spatial coordinates in colData(spe) 
coordinates <- c("array_row", "array_col")
# Specify column names of spatial clusters in colData(spe) 
spatial_cluster <- 'layer_guess_reordered'

# Connect to ExperimentHub
ehub <- ExperimentHub::ExperimentHub()
# Download the full real data (about 2.1 GB in RAM) use:
sce_all <- spatialLIBD::fetch_data(type = "sce", eh = ehub)
set.seed(123)
results <- DESpace_test(spe = sce_one,
                        spatial_cluster = spatial_cluster,
                        verbose = TRUE)
set.seed(123)
cluster_results <- individual_test(sce_one,
                                   spatial_cluster = spatial_cluster,
                                   edgeR_y = results$estimated_y)
results_WM_both <- top_results(cluster_results = cluster_results,
                               cluster = "WM",
                               high_low = "both")

# "MOBP"
(AA1 <- FeaturePlot(sce_one, feature = "ENSG00000168314", spatial_cluster = spatial_cluster, 
            coordinates = coordinates, cluster = 'WM', 
            legend_cluster = FALSE, Annotated_cluster = FALSE, 
            linewidth = 1.4, title = FALSE))

(feature <- rownames(results_WM_both$low_genes)[5]) # "ENSG00000171617"

(AA2 <- FeaturePlot(sce_one, feature, spatial_cluster = spatial_cluster, 
                    coordinates = coordinates, cluster = 'WM', 
                    legend_cluster = FALSE, Annotated_cluster = FALSE, 
                    linewidth = 1.4, title = FALSE))

results_3_both <- top_results(cluster_results = cluster_results,
                              cluster = "Layer3", 
                              high_low = "both")
(feature <- results_3_both$high_genes$gene_id[5])# ENSG00000182601

(AA3 <- FeaturePlot(sce_one, feature, spatial_cluster = spatial_cluster, 
                    coordinates = coordinates, cluster = 'Layer3', 
                    legend_cluster = FALSE, Annotated_cluster = FALSE, 
                    linewidth = 1.4, title = FALSE))

# SVGs with lower than average abundance in WM
(feature <- results_3_both$low_genes$gene_id[10]) # "ENSG00000171476"

(AA4 <- FeaturePlot(sce_one, feature, spatial_cluster = spatial_cluster, 
                    coordinates = coordinates, ncol = 3,cluster = 'Layer3', 
                    legend_cluster = FALSE, Annotated_cluster = FALSE, 
                    linewidth = 1.4, title = FALSE))


AA = ggpubr::ggarrange( AA1 + labs(title = "MOBP")+ theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15)),
                         AA2 + labs(title = "ENC1")+ theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15)),
                         AA4 + labs(title = "HOPX")+ theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15)),
                        AA3 + labs(title = "HS3ST4")+ theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15)),
                         ncol = 4, nrow = 1)

path_save = "./Figures/Figures/Supplementary/"
ggsave(filename = paste0('Individual_cluster_examples.pdf'),
       plot = AA,
       device = "pdf",
       path = path_save,
       width = 15,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

