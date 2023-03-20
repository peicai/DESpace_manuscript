## Figure 4: FDR - multiple samples
rm(list = ls())
path = "~/Desktop/master_thesis/Paper/Data/Simulated/MultiSample/"
path_save = "~/Desktop/master_thesis/Paper/Figures/Simulated/"
source("~/Desktop/master_thesis/Paper/Scripts/Plots_Function.R")
source("~/Desktop/master_thesis/Paper/Scripts/All_methods.R")

spatial_probs = c(0.5,0.8)
# remove the first "Manual_DESpace" color
colours = c(all_colours[2:3], "white")

# points, borders:
shape_border = shape_border[1:2]
# points, fill:
shape_fill = shape_fill[1:2]

## For LIBD
dataset ="LIBD/"
gg_roc_LIBD = gg_fdr_LIBD = list()
#sample_names = c("151507", "151669", "151673")
patterns = c("bottom_pattern","circle_pattern",
             "Manual_clusters_pattern",
             "MixClusters_pattern","MixClusters_reverse_pattern")
for(i in seq_along(patterns)){
  gg = roc_fdr_plot_multi_samples(
    pattern_names = patterns[i],
    path = "~/Desktop/master_thesis/Paper/Data/Simulated/MultiSample/",
    spatial_probs = c(0.5,0.8), colours = colours,
    shape_border = shape_border, shape_fill = shape_fill,
    dataset = dataset)
  gg_roc_LIBD[[i]] = gg[[1]]
  gg_fdr_LIBD[[i]] = gg[[2]]
}

## For melanoma
dataset ="melanoma/"
gg_roc_melanoma = gg_fdr_melanoma = list()
sample_names = c("mel1_rep1", "mel2_rep1", "mel3_rep1", "mel4_rep1")
patterns = c("right_pattern","circle_pattern",
             "BayesSpace_clusters_pattern",
             "MixClusters_pattern","MixClusters_reverse_pattern")
for(i in seq_along(patterns)){
  gg = roc_fdr_plot_multi_samples(
    pattern_names = patterns[i],
    path = "~/Desktop/master_thesis/Paper/Data/Simulated/MultiSample/",
    spatial_probs = c(0.5,0.8), colours = colours,
    shape_border = shape_border, shape_fill = shape_fill,
    dataset = dataset)
  gg_roc_melanoma[[i]] = gg[[1]]
  gg_fdr_melanoma[[i]] = gg[[2]]
}

legend <- ggpubr::get_legend(gg_fdr_LIBD[[1]] +
                               guides(colour = guide_legend(ncol = 2, byrow = FALSE,
                                                            override.aes = list(shape = shape_fill,
                                                                                fill = colours[-length(colours)]) ) ) +
                               theme(legend.key.height=unit(3,"line"),
                                     legend.text=element_text(size=50),
                                     legend.key.width=unit(6,"line")))

ggpubr::as_ggplot(legend)

my_theme = theme(legend.position = "none", 
                 text = element_text(size = 16),
                 plot.title = element_text(hjust = 0, face = "bold", size=50),
                 plot.subtitle = element_text(hjust = 0.5, face = "bold", size=45))


AA1 = ggpubr::ggarrange( gg_fdr_LIBD[[1]] + xlab("") + labs(title = "LIBD", subtitle = "Bottom/Right") + my_theme,
                         gg_fdr_LIBD[[2]] + xlab("") + ylab("") + labs(subtitle = "Circular") + my_theme,
                         gg_fdr_LIBD[[3]] + xlab("") + ylab("") + labs(subtitle = "Annotations") + my_theme,
                         gg_fdr_LIBD[[4]] + xlab("") + ylab("") + labs(subtitle = "Mixture") + my_theme,
                         gg_fdr_LIBD[[5]] + xlab("") + ylab("") + labs(subtitle = "Inverted mixture") + my_theme,
                         ncol = 5, nrow = 1)
AA2 = ggpubr::ggarrange( gg_fdr_melanoma[[1]] + labs(title = "melanoma") + my_theme,
                         gg_fdr_melanoma[[2]] + ylab("")  + my_theme,
                         gg_fdr_melanoma[[3]] + ylab("") + my_theme,
                         gg_fdr_melanoma[[4]] + ylab("") + my_theme,
                         gg_fdr_melanoma[[5]] + ylab("") + my_theme,
                         ncol = 5, nrow = 1)

AA =  ggpubr::ggarrange( AA1, AA2, legend, heights=c(3,3,1),
                         ncol = 1, nrow = 3 )
ggsave(filename = paste0('MultipleSamples_',spatial_probs[1],'_',
                         spatial_probs[2],"_FDR.pdf"),
       plot = AA,
       device = "pdf",
       path = path_save,
       width = 45,
       height = 22,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
