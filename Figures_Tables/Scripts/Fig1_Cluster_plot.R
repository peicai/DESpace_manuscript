library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
#Connect to ExperimentHub
ehub <- ExperimentHub::ExperimentHub()
#Download the example spe data
sce_all <- spatialLIBD::fetch_data(type = "sce", eh = ehub)
#Select one sample only:
LIBD_subset <- sce_all[, colData(sce_all)$sample_name == '151673']

(LIBD_cluster_plot <- clusterPlot(
  LIBD_subset,
  label = "layer_guess_reordered",
  color = NA,
  platform = "Visium",
  is.enhanced = FALSE))

### Melanoma
load("./Data/melanoma/mel1_rep1_melanoma.rda")
melanoma_sce = sce_one
colData(melanoma_sce)$row <- as.numeric(as.character(colData(melanoma_sce)$row))
colData(melanoma_sce)$col <- as.numeric(as.character(colData(melanoma_sce)$col))

(melanoma_cluster_plot <- clusterPlot(
  melanoma_sce,
  label = "spatial.cluster",
  platform = "ST",
  is.enhanced = FALSE)+ scale_y_reverse()+ scale_x_reverse()+ 
    labs(fill = "")+  labs(title = "melanoma",labs(color = "")) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15))+ theme(legend.position="bottom"))

### mouse cerebellum
load("./Data/mouse_cerebellum/cerebellum_filtered_100.rda")
CD <- colData(sce_one)
cell_ID = rownames(colData(sce_one))
colData(sce_one) <- cbind(CD, cell_ID)
library(readr)
stLearn_results <- read_csv("./Data/mouse_cerebellum/100_stLearn_results.csv")
colnames(stLearn_results) <- c("cell_ID","col","row","stLearn_pca_kmeans")

CD <- colData(sce_one)
colData(sce_one) <- cbind(CD, as.factor(stLearn_results$stLearn_pca_kmeans+1))
colnames(colData(sce_one))[16] <- "stLearn_pca_kmeans"
min(colData(sce_one)$row);max(colData(sce_one)$row)
min(colData(sce_one)$col);max(colData(sce_one)$col)

mouse_sce <- sce_one
colData(mouse_sce)$row = colData(sce_one)$row/900
colData(mouse_sce)$col = colData(sce_one)$col/900

(mouse_cluster_plot <- ggplot(as.data.frame(colData(mouse_sce)), 
              aes(x=row, y=col, color=factor(stLearn_pca_kmeans))) +
    geom_point(size = 0.1) + labs(color = "") + scale_y_reverse() + 
    coord_equal()+ guides(color = FALSE, size = FALSE) +
    theme_void()+ 
    labs(title = "mouse cerebellum", fill = "")+ 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15) ))

# function to extract legend from plot
get_only_legend <- function(plot) {
  plot_table <- ggplot_gtable(ggplot_build(plot))
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  legend <- plot_table$grobs[[legend_plot]]
  return(legend)
}

# extract legend from plot1 using above function
legend <- get_only_legend(melanoma_cluster_plot)   

# combine both plots using grid.arrange()
combined_plot <- grid.arrange(melanoma_cluster_plot+ theme(legend.position="none"), mouse_cluster_plot, ncol = 2)

# final combined plot with shared legend
plot1 <- ggplotGrob(LIBD_cluster_plot + labs(title = "LIBD",fill = "",color = "") + 
                      theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)) + 
                      theme(legend.position="bottom"))
plot2 <- grid.arrange(melanoma_cluster_plot+
                      theme(plot.title = element_text(vjust = 15,size = 12), legend.position="none"), 
                      legend,nrow =2,heights = c(7, 1))
plot3 <- grid.arrange(mouse_cluster_plot + labs(title = "mouse cerebellum") + 
                        theme(plot.title = element_text(hjust = 0.5, vjust = 15,face = "bold", size=12)), 
                      legend, nrow = 2, heights = c(7, 1))

path_save = "./Figures/Figures/Real/"
png(filename = paste0(path_save,"Cluster_plot.pdf"), 
    width = 10, height = 4.5, units = "in",
    #pointsize = 12, 
    bg = "white", res = 300)   

grid.arrange(plot1,plot2,plot3, ncol = 3, widths = c(1,1,1))

dev.off()
