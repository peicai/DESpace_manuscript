rm(list = ls())
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Functions.R")
library(data.table)
library(dplyr)
library(iCOBRA)
path_save = "./Figures/Figures/Simulated"
raw_path = "./Simulation/Output/MultiSample/"
spatial_probs = c(0.5,0.8)
# remove the first "Manual_DESpace" color
#colours = c(all_colours[c(2,3,12,14)],"burlywood3", "darkgoldenrod1", "gold", "white")


## For LIBD
dataset ="LIBD/"
gg_roc_LIBD = gg_fdr_LIBD = list()
#sample_names = c("151507", "151669", "151673")
patterns = c("bottom_pattern","circle_pattern",
             "Manual_clusters_pattern",
             "MixClusters_pattern","MixClusters_reverse_pattern")
for(i in seq_along(patterns)){
  gg = overall_roc_fdr_plot_multi_samples(
    pattern_names = patterns[i],
    path = "./multi_samples/",
    spatial_probs = c(0.5,0.8), all_colours = all_colours,
    shape_border = shape_border, shape_fill = shape_fill,
    dataset = dataset)
  gg_roc_LIBD[[i]] = gg[[1]]
  gg_fdr_LIBD[[i]] = gg[[2]]
}

dataset ="melanoma/"
gg_roc_melanoma = gg_fdr_melanoma = list()
sample_names = c("mel1_rep1", "mel2_rep1", "mel3_rep1", "mel4_rep1")
patterns = c("right_pattern","circle_pattern",
             "BayesSpace_clusters_pattern",
             "MixClusters_pattern","MixClusters_reverse_pattern")
for(i in seq_along(patterns)){
  gg = overall_roc_fdr_plot_multi_samples(
    pattern_names = patterns[i],
    path ="./multi_samples/",
    spatial_probs = c(0.5,0.8), all_colours = all_colours,
    shape_border = shape_border, shape_fill = shape_fill,
    dataset = dataset)
  gg_roc_melanoma[[i]] = gg[[1]]
  gg_fdr_melanoma[[i]] = gg[[2]]
}

colours = c(all_colours[c(3,12,14)],"burlywood3", "darkgoldenrod1", "gold", "white")
shape_border = c(0,2,5,6,3,1)
# points, fill:
shape_fill = c(15,17,23,25,3,1)
legend <- ggpubr::get_legend(gg_fdr_LIBD[[5]] +
                               guides(colour = guide_legend(ncol = 4, byrow = FALSE,
                                                            override.aes = list(shape = shape_fill,
                                                                                fill = colours[-length(colours)]) ) ) +
                               theme(legend.key.height=unit(3,"line"),
                                     legend.text=element_text(size=50),
                                     legend.key.width=unit(6,"line")))

ggpubr::as_ggplot(legend)

my_theme = theme(legend.position = "none", 
                 text = element_text(size = 16),
                 aspect.ratio = 0.9,
                 plot.title = element_text(hjust = 0, face = "bold", size=50),
                 plot.subtitle = element_text(hjust = 0.5, face = "bold", size=45),
                 axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(3)),
                 axis.text.y=element_text(size=rel(3)),
                 axis.title.y = element_text(size=rel(3)),
                 axis.title.x = element_text(size=rel(3)))


AA1 = ggpubr::ggarrange( gg_fdr_LIBD[[1]] + xlab("") + labs(title = "LIBD", subtitle = "Bottom/Right") + my_theme,
                         gg_fdr_LIBD[[2]] + xlab("") + ylab("") + labs(subtitle = "Circular") + my_theme,
                         gg_fdr_LIBD[[3]] + xlab("") + ylab("") + labs(subtitle = "Annotations") + my_theme,
                         gg_fdr_LIBD[[4]] + xlab("") + ylab("") + labs(subtitle = "Mixture") + my_theme,
                         gg_fdr_LIBD[[5]] + xlab("") + ylab("") + labs(subtitle = "Inverted mixture") + my_theme,
                         ncol = 5, nrow = 1)
AA2 = ggpubr::ggarrange( gg_fdr_melanoma[[1]] + labs(title = "melanoma", subtitle = "") + my_theme,
                         gg_fdr_melanoma[[2]] + ylab("")  + labs(subtitle = "") + my_theme,
                         gg_fdr_melanoma[[3]] + ylab("")  + labs(subtitle = "") + my_theme,
                         gg_fdr_melanoma[[4]] + ylab("")  + labs(subtitle = "") + my_theme,
                         gg_fdr_melanoma[[5]] + ylab("")  + labs(subtitle = "") + my_theme,
                         ncol = 5, nrow = 1)

AA =  ggpubr::ggarrange( AA1, AA2, legend, heights=c(3,3,1),
                         ncol = 1, nrow = 3 )
ggsave(filename = paste0('MultipleSamples_',spatial_probs[1],'_',
                         spatial_probs[2],"_FDR.pdf"),
       plot = AA,
       device = "pdf",
       #path = path_save,
       width = 45,
       height = 22,
       units = "in",
       dpi = 300,
       limitsize = TRUE)