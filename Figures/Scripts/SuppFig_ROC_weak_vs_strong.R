## Figure 7: ROC - individual samples weak vs. strong
rm(list = ls())
setwd("~/DESpace_manuscript")
path = "./Simulation/Output/"
path_save = "./Figures/Figures/Supplementary/"
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Functions.R")
# remove the first "Manual_DESpace" color and method
colours = c(all_colours[2:14], "white")
methods_all = methods_all[2:14]
## For LIBD
gg_roc_LIBD = gg_fdr_LIBD = list()
sample_names = c("151507", "151669", "151673")
patterns = c("bottom_patch","circle_patch",
             "Manual_clusters_patch",
             "mixture_patch","mixture_reverse_patch")

gg <- list()
setwd("./DESpace_data")
new_data_path <- "./LIBD/"
k=1

for(i in seq_along(patterns)){
  sample_names = c("151507", "151669", "151673")
  if(i %in% c(4,5)) spatial_probs = c(0.6,0.6)
  if(i %in% c(1,2,3)) spatial_probs = c(0.6,0.9)
  if(i == 1) sample_names <- sample_names[-3] # remove the sample 151673
  for(j in seq_along(sample_names)){
    merge_results(j,sample_names = sample_names, 
                  pattern_names = patterns[i],
                  path = path,new_data_path = new_data_path,
                  spatial_probs = spatial_probs)
  }
  gg = overall_roc_fdr_plot(
    sample_names = sample_names, 
    pattern_names = patterns[i],
    path = path,colours = colours,new_data_path = new_data_path,
    spatial_probs = spatial_probs,default = FALSE,rep_i = 1,
    methods_all = methods_all,
    methods_order = methods_order,
    shape_border = shape_border,
    shape_fill = shape_fill)
  
  gg_roc_LIBD[[i]] = gg[[1]]
  gg_fdr_LIBD[[i]] = gg[[2]]
}

## For 151673; bottom pattern: without 'SPARK'
# i <- 3; pattern_names <- c("bottom_patch")
# spatial_probs <- c(0.6, 0.9); default <- FALSE
# methods_all <- methods_all[-3] # remove 'SPARK'
# colours <- colours[-3]

## For melanoma
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Functions.R")
# remove the first "Manual_DESpace" color and method
colours = c(all_colours[2:14], "white")
methods_all = methods_all[2:14]

gg_roc_melanoma = gg_fdr_melanoma = list()
sample_names = c("mel1_rep1", "mel2_rep1", "mel3_rep1", "mel4_rep1")
patterns = c("right_patch","circle_patch",
             "BayesSpace_clusters_patch",
             "mixture_patch","mixture_reverse_patch")
setwd("./DESpace_data")
new_data_path <- "./melanoma/"
for(i in seq_along(patterns)){
  if(i %in% c(4,5)) spatial_probs = c(0.6,0.6)
  if(i %in% c(1,2,3)) spatial_probs = c(0.6,0.9)
  for(j in seq_along(sample_names)){
    merge_results(j,sample_names = sample_names, 
                  pattern_names = patterns[i],
                  path = path,new_data_path = new_data_path,
                  spatial_probs = spatial_probs)
  }
  gg = overall_roc_fdr_plot(
    sample_names = sample_names, 
    pattern_names = patterns[i],
    path = path,colours = colours,new_data_path = new_data_path,
    spatial_probs = spatial_probs,default = FALSE,rep_i = 1,
    methods_all = methods_all,
    methods_order = methods_order,
    shape_border = shape_border,
    shape_fill = shape_fill)
  gg_roc_melanoma[[i]] = gg[[1]]
  gg_fdr_melanoma[[i]] = gg[[2]]
}

my_theme = theme(legend.position = "none", 
                 text = element_text(size = 9),
                 axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(3)),
                 axis.text.y=element_text(size=rel(3)),
                 axis.title.y = element_text(size=rel(3)),
                 axis.title.x = element_text(size=rel(3)),
                 legend.title=element_text(size=rel(2)),
                 legend.text=element_text(size=rel(3)),
                 aspect.ratio = 0.7,
                 plot.title = element_text(hjust = 0.5, face = "bold", size=35),
                 plot.subtitle = element_text(hjust = 0, face = "bold", size=30))

legend <- ggpubr::get_legend(gg_roc_melanoma[[1]] +
                               guides(colour = guide_legend(nrow = 15, byrow = FALSE,
                                                            override.aes = list(shape = shape_fill[-1],size = 2.5,
                                                                                fill = colours[-length(colours)]) ))+#,
                               #linetype = guide_legend(override.aes = list(size = 2))) +
                               theme(legend.key.height=unit(4,"line"),
                                     legend.text=element_text(size=25),
                                     legend.key.width=unit(3,"line")))

ggpubr::as_ggplot(legend)

AA1 = ggpubr::ggarrange(gg_roc_LIBD[[1]]  + xlab("") + labs(title = "LIBD", subtitle = "Bottom/Right") + my_theme,
                        gg_roc_melanoma[[1]]   + xlab("") + ylab("") + labs(title = "melanoma", subtitle = "")+ my_theme,
                        ncol = 2, nrow = 1
) + 
  theme(plot.margin = margin(1, 0, -1, 0, "cm"))

AA2 = ggpubr::ggarrange(gg_roc_LIBD[[2]]   + xlab("")+ labs( subtitle = "Circular") + my_theme,
                        gg_roc_melanoma[[2]]  + xlab("") + ylab("") + labs(subtitle = "")  + my_theme,
                        ncol = 2, nrow = 1
) + 
  theme(plot.margin = margin(1, 0, -1, 0, "cm"))

AA3 = ggpubr::ggarrange(gg_roc_LIBD[[3]]  + xlab("") + labs( subtitle = "Annotations")  + my_theme,
                        gg_roc_melanoma[[3]]  + xlab("") + ylab("") + labs(subtitle = "") + my_theme,
                        ncol = 2, nrow = 1
) + 
  theme(plot.margin = margin(1, 0, -1, 0, "cm"))

AA4 = ggpubr::ggarrange(gg_roc_LIBD[[4]]  + xlab("") + labs( subtitle = "Mixture") + my_theme,
                        gg_roc_melanoma[[4]]  + xlab("") + ylab("") + labs(subtitle = "") + my_theme,
                        ncol = 2, nrow = 1
) + 
  theme(plot.margin = margin(1, 0, -1, 0, "cm"))
AA5 = ggpubr::ggarrange(gg_roc_LIBD[[5]]  + labs( subtitle = "Inverted mixture") + my_theme,
                        gg_roc_melanoma[[5]]  +  ylab("") + labs(subtitle = "") + my_theme,
                        ncol = 2, nrow = 1
) 


AA = ggpubr::ggarrange( ggarrange(AA1, AA2, AA3, AA4, AA5,
                                  heights=c(1,1,1,1,1),
                                  nrow = 5), legend, 
                        # font.label = list(color = "black", size = 40, face = "bold"),
                        #label.x = 0, labrl.y = 1,
                        #labels = c("Bottom/Right", "Circular", "Annotations", "Mixture", "Inverted mixture"),
                        widths = c(3,1.3),
                        nrow = 1 ) 

spatial_probs = c(0.6,0.9)
ggsave(filename = paste0('Aggregated_Weak_vs_strong_',spatial_probs[1],'_',
                         spatial_probs[2],"_ROC.pdf"),
       plot = AA,
       device = "pdf",
       width = 20,
       height = 30,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

