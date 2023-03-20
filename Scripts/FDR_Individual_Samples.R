## Figure 3: FDR - individual samples
rm(list = ls())
setwd("~/Desktop/master_thesis/Paper/Data/Simulated")
path = "~/Desktop/master_thesis/Paper/Data/Simulated/"
path_save = "~/Desktop/master_thesis/Paper/Figures/Simulated/"
source("~/Desktop/master_thesis/Paper/Scripts/All_methods.R")
source("~/Desktop/master_thesis/Paper/Scripts/Plots_Function.R")
# remove the first "Manual_DESpace" color and method
colours = c(all_colours[-1], "white")
methods_all = methods_all[-1]
## For LIBD
gg_roc_LIBD = gg_fdr_LIBD = list()
sample_names = c("151507", "151669", "151673")
patterns = c("bottom_patch","circle_patch",
             "Manual_clusters_patch",
             "mixture_patch","mixture_reverse_patch")
for(i in seq_along(patterns)){
  gg = overall_roc_fdr_plot(
              sample_names = sample_names, 
              pattern_names = patterns[i],
              path = path,colours = colours,path_save = path_save,
              spatial_probs = c(0.5,0.9),default = FALSE,rep_i = 1,
              methods_all = methods_all,
              methods_order = methods_order,
              shape_border = shape_border,
              shape_fill = shape_fill,
              # roc_save_name='_ROC_final.pdf',
              # fdr_save_name = '_FDR_final.pdf',
              dataset = 'LIBD')
  gg_roc_LIBD[[i]] = gg[[1]]
  gg_fdr_LIBD[[i]] = gg[[2]]
  
  }

## For melanoma
gg_roc_melanoma = gg_fdr_melanoma = list()
sample_names = c("mel1_rep1", "mel2_rep1", "mel3_rep1", "mel4_rep1")
patterns = c("right_patch","circle_patch",
             "BayesSpace_clusters_patch",
             "mixture_patch","mixture_reverse_patch")
for(i in seq_along(patterns)){
  gg = overall_roc_fdr_plot(
    sample_names = sample_names, 
    pattern_names = patterns[i],
    path = path,colours = colours,path_save = path_save,
    spatial_probs = c(0.5,0.9),default = FALSE,rep_i = 1,
    methods_all = methods_all,
    methods_order = methods_order,
    shape_border = shape_border,
    shape_fill = shape_fill,
    # roc_save_name='_ROC_final.pdf',
    # fdr_save_name = '_FDR_final.pdf',
    dataset = 'melanoma')
  gg_roc_melanoma[[i]] = gg[[1]]
  gg_fdr_melanoma[[i]] = gg[[2]]
  
}

## For mouse_cerebellum
gg_roc_cerebellum = gg_fdr_cerebellum = list()
patterns = c("bottom_patch","circle_patch",
             "StLearn_clusters_patch",
             "mixture_patch","mixture_reverse_patch")
for(i in seq_along(patterns[1:3])){
  gg = roc_fdr_plot_cerebellum(
    sample_names = sample_names, 
    pattern_names = patterns[i],
    path = path,colours = colours,path_save = path_save,
    spatial_probs = c(0.5,0.9),default = FALSE,rep_i = 1,
    methods_all = methods_all,
    methods_order = methods_order,
    shape_border = shape_border,
    shape_fill = shape_fill)
  gg_roc_cerebellum[[i]] = gg[[1]]
  gg_fdr_cerebellum[[i]] = gg[[2]]
  
}
## remove SpaGCN as it's failed
for(i in c(4,5)){
  gg = roc_fdr_plot_cerebellum(
    sample_names = sample_names, 
    pattern_names = patterns[i],
    path = path,colours = colours[-8],path_save = path_save,
    spatial_probs = c(0.5,0.9),default = FALSE,rep_i = 1,
    methods_all = methods_all[-8],
    methods_order = methods_order[-8],
    shape_border = shape_border[-8],
    shape_fill = shape_fill[-8])
  
  gg_roc_cerebellum[[i]] = gg[[1]]
  gg_fdr_cerebellum[[i]] = gg[[2]]
  
}
### Common legend
legend <- ggpubr::get_legend(gg_fdr_LIBD[[1]] +
                               guides(colour = guide_legend(nrow = 2, byrow = FALSE,
                                                            override.aes = list(shape = shape_fill,
                                                                                fill = colours[-length(colours)]) ) ) +
                               theme(legend.key.height=unit(3,"line"),
                                     legend.text=element_text(size=50),
                                     legend.key.width=unit(6,"line")))

ggpubr::as_ggplot(legend)

my_theme = theme(legend.position = "none", 
                 text = element_text(size = 9),
                 plot.title = element_text(hjust = 0, face = "bold", size=55),
                 plot.subtitle = element_text(hjust = 0.5, face = "bold", size=50))


AA1 = ggpubr::ggarrange( gg_fdr_LIBD[[1]] + xlab("") + labs(title = "LIBD", subtitle = "Bottom/Right") + my_theme,
                            gg_fdr_LIBD[[2]] + xlab("") + ylab("") + labs(subtitle = "Circular") + my_theme,
                            gg_fdr_LIBD[[3]] + xlab("") + ylab("") + labs(subtitle = "Annotations") + my_theme,
                            gg_fdr_LIBD[[4]] + xlab("") + ylab("") + labs(subtitle = "Mixture") + my_theme,
                            gg_fdr_LIBD[[5]] + xlab("") + ylab("") + labs(subtitle = "Inverted mixture") + my_theme,
                         # hjust=-0.8,
                         # labels = "LIBD",font.label = list(color = "black", size = 40, face = "bold"),
                     ncol = 5, nrow = 1)
AA2 = ggpubr::ggarrange( gg_fdr_melanoma[[1]] + xlab("") + labs(title = "melanoma") + my_theme,
                             gg_fdr_melanoma[[2]] + xlab("") + ylab("") + my_theme,
                             gg_fdr_melanoma[[3]] + xlab("") + ylab("") + my_theme,
                             gg_fdr_melanoma[[4]] + xlab("") + ylab("") + my_theme,
                             gg_fdr_melanoma[[5]] + xlab("") + ylab("") + my_theme,
                         # hjust=-0.8,
                         # labels = "melanoma",font.label = list(color = "black", size = 40, face = "bold"),
                      ncol = 5, nrow = 1)
AA3 = ggpubr::ggarrange(gg_fdr_cerebellum[[1]] + labs(title = "mouse cerebellum") + my_theme,
                             gg_fdr_cerebellum[[2]] + ylab("") + my_theme,
                             gg_fdr_cerebellum[[3]] + ylab("") + my_theme,
                             gg_fdr_cerebellum[[4]] + ylab("") + my_theme,
                             gg_fdr_cerebellum[[5]] + ylab("") + my_theme,
                        # hjust=-0.8,
                        # labels = "mouse cerebellum",font.label = list(color = "black", size = 40, face = "bold"),
                      ncol = 5, nrow = 1)

AA =  ggpubr::ggarrange( AA1, AA2, AA3, legend, 
                        # font.label = list(color = "black", size = 40, face = "bold"),
                         #label.x = 0, labrl.y = 1,
                         #labels = c("LIBD", "melanoma", "mouse cerebellum"),
                     heights=c(3,3,3,1),
                     ncol = 1, nrow = 4 )
# BB = annotate_figure(
#   AA,left = text_grob("         mouse cerebellum                               melanoma                                                LIBD      ", 
#                       color = "red", face = "bold", size = 40,rot = 90))
spatial_probs = c(0.5,0.9)
ggsave(filename = paste0('AggregatedIndividualSamples_',spatial_probs[1],'_',
                         spatial_probs[2],"_FDR.pdf"),
       plot = AA,
       device = "pdf",
       path = path_save,
       width = 45,
       height = 35,
       units = "in",
       dpi = 300,
       limitsize = TRUE)


# AA = egg::ggarrange( plots = 
#                        list(gg_fdr_LIBD[[1]] + xlab("") + labs(title = "Bottom/Right") + my_theme,
#                             gg_fdr_LIBD[[2]] + xlab("") + ylab("") + labs(title = "Circular") + my_theme,
#                             gg_fdr_LIBD[[3]] + xlab("") + ylab("") + labs(title = "Annotations") + my_theme,
#                             gg_fdr_LIBD[[4]] + xlab("") + ylab("") + labs(title = "Mixture") + my_theme,
#                             gg_fdr_LIBD[[5]] + xlab("") + ylab("") + labs(title = "Inverted mixture") + my_theme,
#                             
#                             gg_fdr_melanoma[[1]] + xlab("") + my_theme,
#                             gg_fdr_melanoma[[2]] + xlab("") + ylab("") + my_theme,
#                             gg_fdr_melanoma[[3]] + xlab("") + ylab("") + my_theme,
#                             gg_fdr_melanoma[[4]] + xlab("") + ylab("") + my_theme,
#                             gg_fdr_melanoma[[5]] + xlab("") + ylab("") + my_theme,
#                             
#                             gg_fdr_cerebellum[[1]] + my_theme,
#                             gg_fdr_cerebellum[[2]]  + ylab("") + my_theme,
#                             gg_fdr_cerebellum[[3]]  + ylab("") + my_theme,
#                             gg_fdr_cerebellum[[4]]  + ylab("") + my_theme,
#                             gg_fdr_cerebellum[[5]]  + ylab("") + my_theme),
#                      bottom = legend,
#                      ncol = 5, nrow = 3)
# 
# ggsave(filename = paste0('AggregatedSamples_',spatial_probs[1],'_',
#                          spatial_probs[2],"_FDR.pdf"),
#        plot = AA,
#        device = "pdf",
#        path = path_save,
#        width = 35,
#        height = 30,
#        units = "in",
#        dpi = 300,
#        limitsize = TRUE)