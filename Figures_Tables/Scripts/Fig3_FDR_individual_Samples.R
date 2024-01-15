rm(list = ls())
setwd("~/DESpace_manuscript")
path = "./Simulation/Output/"
path_save = "./Figures/Figures/Supplementary/"
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Functions.R")
library(stats)
# remove the first "Manual_DESpace" color and method
colours = c(all_colours[2:14], "white")
methods_all = methods_all[2:14]
#methods_order = methods_order
## For LIBD
gg_roc_LIBD = gg_fdr_LIBD = list()
sample_names = c("151507", "151669", "151673")
patterns = c("bottom_patch","circle_patch",
             "Manual_clusters_patch",
             "mixture_patch","mixture_reverse_patch")
k <- 1
setwd("./DESpace_data")
new_data_path <- "./LIBD/"
spatial_probs = c(0.5,0.9)
default = FALSE
rep = 1
for(i in seq_along(patterns)){
  gg = overall_roc_fdr_plot(
    sample_names = sample_names, 
    pattern_names = patterns[i],
    path = path,colours = colours,new_data_path = new_data_path,
    spatial_probs = c(0.5,0.9),default = FALSE,rep_i = 1,
    methods_all = methods_all,
    methods_order = methods_order,
    shape_border = shape_border,
    shape_fill = shape_fill)
  gg_roc_LIBD[[i]] = gg[[1]]
  gg_fdr_LIBD[[i]] = gg[[2]]
}


# for melanoma
gg_roc_melanoma = gg_fdr_melanoma = list()
sample_names = c("mel1_rep1", "mel2_rep1", "mel3_rep1", "mel4_rep1")
patterns = c("right_patch","circle_patch",
             "BayesSpace_clusters_patch",
             "mixture_patch","mixture_reverse_patch")
k <- 1
setwd("./DESpace_data")
new_data_path <- "./melanoma/"
for(i in seq_along(patterns)){
  gg = overall_roc_fdr_plot(
    sample_names = sample_names, 
    pattern_names = patterns[i],
    path = path,colours = colours,new_data_path = new_data_path,
    spatial_probs = c(0.5,0.9),default = FALSE,rep_i = 1,
    methods_all = methods_all,
    methods_order = methods_order,
    shape_border = shape_border,
    shape_fill = shape_fill)
  gg_roc_melanoma[[i]] = gg[[1]]
  gg_fdr_melanoma[[i]] = gg[[2]]
}

## slideseq v2
gg_roc_cerebellum = gg_fdr_cerebellum = list()
patterns = c("bottom_patch","circle_patch",
             "StLearn_clusters_patch",
             "mixture_patch","mixture_reverse_patch")
new_data_path <- "~/Desktop/SlideSeq2"
setwd("./DESpace_data")
source("./Figures/Scripts/Plots_Functions.R")

k <- 1
for(i in 1:5){
  spatial_probs = c(0.5,0.9)
  merge_results(j=NULL,sample_names = NULL, 
                pattern_names = patterns[i],
                path = path,new_data_path = new_data_path,
                spatial_probs = spatial_probs)
}
for(i in 1:3){
  gg <- overall_roc_fdr_plot_cerebellum(
    pattern_names = patterns[i],
    path = path,colours = colours,new_data_path = new_data_path,
    spatial_probs = c(0.5,0.9),default = FALSE,rep_i = 1,
    methods_all = methods_all,
    methods_order = methods_order,
    shape_border = shape_border,
    shape_fill = shape_fill)
  gg_roc_cerebellum[[i]] = gg[[1]]
  gg_fdr_cerebellum[[i]] = gg[[2]]
}
setwd("~/DESpace_manuscript")
source("./Figures/Scripts/Plots_Functions.R")
colours = c(all_colours[c(2:8,10:14)], "white")
methods_all = methods_all[-c(1,9)];
methods_order = methods_order[-8]
setwd("./DESpace_data")
for(i in 4:5){
  gg <- overall_roc_fdr_plot_cerebellum(
    pattern_names = patterns[i],
    path = path,colours = colours,new_data_path = new_data_path,
    spatial_probs = c(0.5,0.9),default = FALSE,rep_i = 1,
    methods_all = methods_all,
    methods_order = methods_order,
    shape_border = shape_border,
    shape_fill = shape_fill)
  gg_roc_cerebellum[[i]] = gg[[1]]
  gg_fdr_cerebellum[[i]] = gg[[2]]
}
source("./Figures/Scripts/Plots_Functions.R")
colours = c(all_colours[2:14], "white")
my_theme = theme(legend.position = "none", aspect.ratio = 1,
                 text = element_text(size = 9),
                 plot.title = element_text(hjust = 0, face = "bold", size=55),
                 plot.subtitle = element_text(hjust = 0.5, face = "bold", size=50))
legend <- ggpubr::get_legend(gg_fdr_LIBD[[1]] +
                               guides(colour = guide_legend(nrow = 3, byrow = FALSE,
                                                            override.aes = list(shape = shape_fill[2:14],
                                                                                fill = colours[-length(colours)]) 
                               ) ) +
                               theme(legend.key.height=unit(3,"line"),
                                     legend.text=element_text(size=50),
                                     legend.key.width=unit(6,"line")))

ggpubr::as_ggplot(legend)

AA1 = ggpubr::ggarrange( gg_fdr_LIBD[[1]] + xlab("") + labs(title = "LIBD", subtitle = "Bottom/Right") + my_theme,
                         gg_fdr_LIBD[[2]] + xlab("") + ylab("") + labs(subtitle = "Circular") + my_theme,
                         gg_fdr_LIBD[[3]] + xlab("") + ylab("") + labs(subtitle = "Annotations") + my_theme,
                         gg_fdr_LIBD[[4]] + xlab("") + ylab("") + labs(subtitle = "Mixture") + my_theme,
                         gg_fdr_LIBD[[5]] + xlab("") + ylab("") + labs(subtitle = "Inverted mixture") + my_theme,
                         ncol = 5, nrow = 1)
AA2 = ggpubr::ggarrange( gg_fdr_melanoma[[1]] + xlab("") + labs(title = "melanoma") + my_theme,
                         gg_fdr_melanoma[[2]] + xlab("") + ylab("") + my_theme,
                         gg_fdr_melanoma[[3]] + xlab("") + ylab("") + my_theme,
                         gg_fdr_melanoma[[4]] + xlab("") + ylab("") + my_theme,
                         gg_fdr_melanoma[[5]] + xlab("") + ylab("") + my_theme,
                         ncol = 5, nrow = 1)
AA3 = ggpubr::ggarrange(gg_fdr_cerebellum[[1]] + labs(title = "mouse cerebellum") + my_theme,
                        gg_fdr_cerebellum[[2]] + ylab("") + my_theme,
                        gg_fdr_cerebellum[[3]] + ylab("") + my_theme,
                        gg_fdr_cerebellum[[4]] + ylab("") + my_theme,
                        gg_fdr_cerebellum[[5]] + ylab("") + my_theme,
                        ncol = 5, nrow = 1)

AA =  ggpubr::ggarrange( AA1, AA2, AA3, legend, 
                         heights=c(3,3,3,1),
                         ncol = 1, nrow = 4 )

spatial_probs = c(0.5,0.9)
ggsave(filename = paste0('AggregatedIndividualSamples_',spatial_probs[1],'_',
                         spatial_probs[2],"_FDR.pdf"),
       plot = AA,
       device = "pdf",
       width = 45,
       height = 35,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
