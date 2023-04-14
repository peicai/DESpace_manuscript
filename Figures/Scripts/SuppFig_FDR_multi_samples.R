## Supp.Figure 11: FDR - multiple samples
rm(list = ls())
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Function.R")
path_save = "./Figures/Figures/Supplementary/"
path = "./Simulation/Output/"

# remove the first "Manual_DESpace" color and method
colours = c(all_colours[-1], "white")
methods_all = methods_all[-1]
## For LIBD
gg_roc_LIBD = gg_fdr_LIBD = list()
sample_names = c("151507", "151669", "151673")
patterns = c("bottom_pattern","circle_pattern",
             "Manual_clusters_pattern",
             "MixClusters_pattern","MixClusters_reverse_pattern")
dataset = "LIBD"
colours = c(all_colours[2:3], "white")

# points, borders:
shape_border = shape_border[1:2]
# points, fill:
shape_fill = shape_fill[1:2]
k <- 1
for(i in seq_along(patterns)){
  for(j in seq_along(sample_names)){
    gg = roc_fdr_plot_multi_samples(
      j,
      sample_names,
      pattern_names = patterns[i],
      path = "./Simulation/Output/MultiSample/",
      spatial_probs = c(0.5,0.8), colours = colours,
      shape_border = shape_border, shape_fill = shape_fill,
      dataset = dataset)
    gg_roc_LIBD[[k]] = gg[[1]]
    gg_fdr_LIBD[[k]] = gg[[2]]
    k <- k + 1
  }  
}


my_theme = theme(legend.position = "none", 
                 text = element_text(size = 9),
                 aspect.ratio = 0.8,
                 plot.title = element_text(hjust = 0, face = "bold", size=45),
                 plot.subtitle = element_text(hjust = 0.5, face = "bold", size=40),
                 axis.text.x = element_text(size = rel(4)),
                 axis.text.y=element_text(size=rel(4)),
                 axis.title.y = element_text(size=rel(5)),
                 axis.title.x = element_text(size=rel(5)))


legend <- ggpubr::get_legend(gg_fdr_LIBD[[1]] +
                               guides(colour = guide_legend(nrow = 1, byrow = FALSE,
                                                            override.aes = list(shape = shape_fill,
                                                                                fill = colours[-length(colours)]) ) ) +
                               theme(legend.key.height=unit(3,"line"),
                                     legend.text=element_text(size=40),
                                     legend.key.width=unit(6,"line")))

ggpubr::as_ggplot(legend)

AA1 = ggpubr::ggarrange(gg_fdr_LIBD[[1]] + xlab("") + labs(title = "Sample 151507", subtitle = "Bottom/Right") + my_theme,
                        gg_fdr_LIBD[[4]] + xlab("") + ylab("") + labs(subtitle = "Circular") + my_theme,
                        gg_fdr_LIBD[[7]] + xlab("") + ylab("") + labs(subtitle = "Annotations") + my_theme,
                        gg_fdr_LIBD[[10]] + xlab("") + ylab("") + labs(subtitle = "Mixture") + my_theme,
                        gg_fdr_LIBD[[13]] + xlab("") + ylab("") + labs(subtitle = "Inverted mixture") + my_theme,
                        
                        gg_fdr_LIBD[[2]] + xlab("") + labs(title = "Sample 151669", subtitle = "") + my_theme,
                        gg_fdr_LIBD[[5]] + xlab("") + ylab("") + labs(subtitle = "") + my_theme,
                        gg_fdr_LIBD[[8]] + xlab("") + ylab("") + labs(subtitle = "") + my_theme,
                        gg_fdr_LIBD[[11]] + xlab("") + ylab("") + labs(subtitle = "") + my_theme,
                        gg_fdr_LIBD[[14]] + xlab("") + ylab("") + labs(subtitle = "") + my_theme,
                        
                        gg_fdr_LIBD[[3]]+ labs(title = "Sample 151673", subtitle = "") + my_theme,
                        gg_fdr_LIBD[[6]]  + ylab("") + labs(subtitle = "") + my_theme,
                        gg_fdr_LIBD[[9]] + ylab("")+ labs(subtitle = "") + my_theme,
                        gg_fdr_LIBD[[12]] + ylab("")+ labs(subtitle = "") + my_theme,
                        gg_fdr_LIBD[[15]] + ylab("")+ labs(subtitle = "") + my_theme,
                        
                        ncol = 5, nrow = 3)
AA = ggpubr::ggarrange( AA1, legend, 
                        # font.label = list(color = "black", size = 40, face = "bold"),
                        #label.x = 0, labrl.y = 1,
                        #labels = c("LIBD", "melanoma", "mouse cerebellum"),
                        heights=c(9,1),
                        ncol = 1, nrow = 2 )

spatial_probs = c(0.5,0.8)
ggsave(filename = paste0('MultipleSamples_',spatial_probs[1],'_',
                         spatial_probs[2],"_FDR_LIBD.pdf"),
       plot = AA,
       device = "pdf",
       path = path_save,
       width = 35,
       height = 25,
       units = "in",
       dpi = 300,
       limitsize = TRUE)


## For melanoma
rm(list = ls())
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Function.R")
path_save = "./Figures/Figures/Supplementary/"
path = "./Simulation/Output/"
# remove the first "Manual_DESpace" color and method
colours = c(all_colours[2:3], "white")

# points, borders:
shape_border = shape_border[1:2]
# points, fill:
shape_fill = shape_fill[1:2]
dataset = "melanoma"
gg_roc_melanoma = gg_fdr_melanoma = list()
sample_names = c("mel1_rep1", "mel2_rep1", "mel3_rep1", "mel4_rep1")
patterns = c("right_pattern","circle_pattern",
             "BayesSpace_clusters_pattern",
             "MixClusters_pattern","MixClusters_reverse_pattern")
k=1
for(i in seq_along(patterns)){
  for(j in seq_along(sample_names)){
    gg = roc_fdr_plot_multi_samples(
      j,
      sample_names,
      pattern_names = patterns[i],
      path = "./Simulation/Output/MultiSample/",
      spatial_probs = c(0.5,0.8), colours = colours,
      shape_border = shape_border, shape_fill = shape_fill,
      dataset = dataset)
    gg_roc_melanoma[[k]] = gg[[1]]
    gg_fdr_melanoma[[k]] = gg[[2]]
    k <- k + 1
  }  
}
my_theme = theme(legend.position = "none", 
                 text = element_text(size = 9),
                 aspect.ratio = 0.8,
                 plot.title = element_text(hjust = 0, face = "bold", size=50),
                 plot.subtitle = element_text(hjust = 0.5, face = "bold", size=45),
                 axis.text.x = element_text(size = rel(4)),
                 axis.text.y=element_text(size=rel(4)),
                 axis.title.y = element_text(size=rel(6)),
                 axis.title.x = element_text(size=rel(6)))
legend <- ggpubr::get_legend(gg_fdr_melanoma[[1]] +
                               guides(colour = guide_legend(nrow = 1, byrow = FALSE,
                                                            override.aes = list(shape = shape_fill,
                                                                                fill = colours[-length(colours)]) ) ) +
                               theme(legend.key.height=unit(3,"line"),
                                     legend.text=element_text(size=50),
                                     legend.key.width=unit(6,"line")))

ggpubr::as_ggplot(legend)

AA1 = ggpubr::ggarrange(gg_fdr_melanoma[[1]] + xlab("") + labs(title = "Sample mel1_rep1", subtitle = "Bottom/Right") + my_theme,
                        gg_fdr_melanoma[[5]] + xlab("") + ylab("") + labs(subtitle = "Circular") + my_theme,
                        gg_fdr_melanoma[[9]] + xlab("") + ylab("") + labs(subtitle = "Annotations") + my_theme,
                        gg_fdr_melanoma[[13]] + xlab("") + ylab("") + labs(subtitle = "Mixture") + my_theme,
                        gg_fdr_melanoma[[17]] + xlab("") + ylab("") + labs(subtitle = "Inverted mixture") + my_theme,
                        ncol = 5, nrow = 1)

AA2 = ggpubr::ggarrange(gg_fdr_melanoma[[2]] + xlab("") + labs(title = "Sample mel2_rep1", subtitle = "") + my_theme,
                        gg_fdr_melanoma[[6]] + xlab("") + ylab("")+ labs(subtitle = "")  + my_theme,
                        gg_fdr_melanoma[[10]] + xlab("") + ylab("")+ labs(subtitle = "") + my_theme,
                        gg_fdr_melanoma[[14]] + xlab("") + ylab("") + labs(subtitle = "")+ my_theme,
                        gg_fdr_melanoma[[18]] + xlab("") + ylab("")+ labs(subtitle = "") + my_theme,
                        ncol = 5, nrow = 1)

AA3 = ggpubr::ggarrange(gg_fdr_melanoma[[3]] + xlab("") + labs(title = "Sample mel3_rep1", subtitle = "")  + my_theme,
                        gg_fdr_melanoma[[7]] + xlab("") + ylab("") + labs(subtitle = "")+ my_theme,
                        gg_fdr_melanoma[[11]] + xlab("") + ylab("")+ labs(subtitle = "") + my_theme,
                        gg_fdr_melanoma[[15]] + xlab("") + ylab("") + labs(subtitle = "")+ my_theme,
                        gg_fdr_melanoma[[19]] + xlab("") + ylab("") + labs(subtitle = "")+ my_theme,
                        ncol = 5, nrow = 1)

AA4 = ggpubr::ggarrange(gg_fdr_melanoma[[4]] + labs(title = "Sample mel4_rep1", subtitle = "")  + my_theme,
                        gg_fdr_melanoma[[8]] + ylab("") + labs(subtitle = "")+ my_theme,
                        gg_fdr_melanoma[[12]] + ylab("")+ labs(subtitle = "") + my_theme,
                        gg_fdr_melanoma[[16]] + ylab("")+ labs(subtitle = "") + my_theme,
                        gg_fdr_melanoma[[20]] + ylab("")+ labs(subtitle = "") + my_theme,
                        ncol = 5, nrow = 1)
AA = ggpubr::ggarrange( AA1, AA2, AA3, AA4, legend, 
                        # font.label = list(color = "black", size = 40, face = "bold"),
                        #label.x = 0, labrl.y = 1,
                        #labels = c("LIBD", "melanoma", "mouse cerebellum"),
                        heights=c(3,3,3,3,1),
                        widths = c(1,1,1,1,1),
                        ncol = 1, nrow = 5 )

spatial_probs = c(0.5,0.8)
ggsave(filename = paste0('MultipleSamples_',spatial_probs[1],'_',
                         spatial_probs[2],"_FDR_melanoma.pdf"),
       plot = AA,
       device = "pdf",
       path = path_save,
       width = 45,
       height = 35,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
