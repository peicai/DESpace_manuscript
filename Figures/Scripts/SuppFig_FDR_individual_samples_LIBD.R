rm(list = ls())
setwd("~/DESpace_manuscript")
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Functions.R")
library(stats)
path_save = "./Figures/Figures/Supplementary/"
path = "./Simulation/Output/"
# remove the first "Manual_DESpace" color and method
colours = c(all_colours[-1], "white")
# methods_all = methods_all[-1]
## For LIBD
gg_roc_LIBD = gg_fdr_LIBD = gg = list()
sample_names = c("151507", "151669", "151673")
patterns = c("bottom_patch","circle_patch",
             "Manual_clusters_patch",
             "mixture_patch","mixture_reverse_patch")
k <- 1
setwd("./DESpace_data")
new_data_path <- "./LIBD/"
for(i in seq_along(patterns)){
  for(j in seq_along(sample_names)){
    spatial_probs = c(0.5,0.9)
    merge_results(j,sample_names = sample_names, 
                                   pattern_names = patterns[i],
                                   path = path,new_data_path = new_data_path,
                                   spatial_probs = spatial_probs)
    gg <- roc_fdr_plot(j,
                             sample_names = sample_names, 
                             pattern_names = patterns[i],
                             path = path,colours = colours,new_data_path = new_data_path,
                             spatial_probs = spatial_probs,
                             methods_all = methods_all[-1],
                             methods_order = methods_order,
                             shape_border = shape_border,
                             shape_fill = shape_fill)
    gg_roc_LIBD[[k]] = gg[[1]]
    gg_fdr_LIBD[[k]] = gg[[2]]
    k <- k + 1
  }  
}

my_theme = theme(legend.position = "none", aspect.ratio = 1,
                 text = element_text(size = 9),
                 plot.title = element_text(hjust = 0, face = "bold", size=55),
                 plot.subtitle = element_text(hjust = 0.5, face = "bold", size=50))
legend <- ggpubr::get_legend(gg_fdr_LIBD[[1]] +
                              guides(colour = guide_legend(nrow = 3, byrow = FALSE,
                                                            override.aes = list(shape = shape_fill[-1],
                                                                                fill = colours[-c(length(colours))]) 
                                                           ) ) +
                               theme(legend.key.height=unit(3,"line"),
                                     legend.text=element_text(size=50),
                                     legend.key.width=unit(6,"line")))

ggpubr::as_ggplot(legend)

AA1 = ggpubr::ggarrange(gg_fdr_LIBD[[1]] + xlab("") + labs(title = "Sample 151507", subtitle = "Bottom/Right") + my_theme,
                        gg_fdr_LIBD[[2]] + xlab("") + ylab("") + labs(subtitle = "Circular") + my_theme,
                        gg_fdr_LIBD[[3]] + xlab("") + ylab("") + labs(subtitle = "Annotations") + my_theme,
                        gg_fdr_LIBD[[4]] + xlab("") + ylab("") + labs(subtitle = "Mixture") + my_theme,
                        gg_fdr_LIBD[[5]] + xlab("") + ylab("") + labs(subtitle = "Inverted mixture") + my_theme,
                        
                        gg_fdr_LIBD[[6]] + xlab("") + labs(title = "Sample 151669") + my_theme,
                        gg_fdr_LIBD[[7]] + xlab("") + ylab("") + my_theme,
                        gg_fdr_LIBD[[8]] + xlab("") + ylab("") + my_theme,
                        gg_fdr_LIBD[[9]] + xlab("") + ylab("") + my_theme,
                        gg_fdr_LIBD[[10]] + xlab("") + ylab("") + my_theme,
                        
                        gg_fdr_LIBD[[11]] + xlab("") + labs(title = "Sample 151673") + my_theme,
                        gg_fdr_LIBD[[12]] + xlab("") + ylab("") + my_theme,
                        gg_fdr_LIBD[[13]] + xlab("") + ylab("") + my_theme,
                        gg_fdr_LIBD[[14]] + xlab("") + ylab("") + my_theme,
                        gg_fdr_LIBD[[15]] + xlab("") + ylab("") + my_theme,
                        
                        ncol = 5, nrow = 3)
AA = ggpubr::ggarrange( AA1, legend, 
                        # font.label = list(color = "black", size = 40, face = "bold"),
                        #label.x = 0, labrl.y = 1,
                        #labels = c("LIBD", "melanoma", "mouse cerebellum"),
                        heights=c(9,1),
                        ncol = 1, nrow = 2 )

spatial_probs = c(0.5,0.9)
ggsave(filename = paste0('IndividualSamples_',spatial_probs[1],'_',
                         spatial_probs[2],"_FDR_LIBD.pdf"),
       plot = AA,
       device = "pdf",
       #path = path_save,
       width = 45,
       height = 35,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
