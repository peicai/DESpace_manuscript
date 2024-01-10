rm(list = ls())
<<<<<<< HEAD:Figures_Tables/Scripts/Fig_Barplot_Jarccard_index.R
=======
# library(S4Vectors)
# library(stringr)
# library(dplyr)
>>>>>>> origin/main:Figures/Scripts/Fig_Barplot_Jarccard_index.R
path = "./Real_data/results/"
path_save = "./Figures/Figures/Supplementary/"
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Functions.R")
############################ For LIBD ###################################
data = 'LIBD'
sample_names = c("151507", "151508", "151509", "151510",
                 "151669", "151670", "151671", "151672",
                 "151673", "151674", "151675", "151676")
list_result_LIBD <- list()
<<<<<<< HEAD:Figures_Tables/Scripts/Fig_Barplot_Jarccard_index.R
##################### scripts for generating "All_methods_genes_order.rda" ############################

=======
>>>>>>> origin/main:Figures/Scripts/Fig_Barplot_Jarccard_index.R
# list_result_LIBD <- process_data(
#    path = path,
#    data = 'LIBD', Manual = TRUE,
#    threshold = 2, # indicates that we do not set threshold
#    name_head = -10,name_tail = -5
# )
# names(list_result_LIBD) <- c("Manual_DESpace","BayesSpace_DESpace", "StLearn_DESpace",
#                               "Manual_scran","BayesSpace_scran", "StLearn_scran",
#                               "Manual_seurat","BayesSpace_seurat", "StLearn_seurat",
#                               "MERINGUE", "SPARK","nnSVG", "SPARK-X", "SpatialDE2","SpatialDE","SpaGCN",
#                               "All_results")

#save(list_result_LIBD, file = "./Real_data/results/LIBD/All_methods_genes_order.rda")
colors_method2 <- c(colors_method[1:10], "orchid",colors_method[11:12],"aquamarine",colors_method[13:14])
names(colors_method2)[c(11,14)] <- c("Manual_findMarkers", "Manual_FindAllMarkers")

load("./Real_data/results/LIBD/All_methods_genes_order.rda")
gg_Bar_LIBD <- ave_jaccard_btw_rep_LIBD(
  list_test_all=list_result_LIBD,
  path = path_save,
  sample_names=sample_names,
  top_genes = 1000,
  col_method = colors_method2,
  data = 'LIBD')


###################################### For melanoma ################
sample_names = c("mel1_rep1", "mel1_rep2",
                 "mel2_rep1", "mel2_rep2",
                 "mel3_rep1", "mel3_rep2",
                 "mel4_rep1", "mel4_rep2")
list_result_melanoma <- list()
<<<<<<< HEAD:Figures_Tables/Scripts/Fig_Barplot_Jarccard_index.R

##################### scripts for generating "All_methods_genes_order.rda" ############################
=======
>>>>>>> origin/main:Figures/Scripts/Fig_Barplot_Jarccard_index.R
#  list_result_melanoma <- process_data(
#    path = path,
#    data = 'melanoma', Manual = FALSE,
#    threshold = 2, # indicates that we do not set threshold
#    name_head = -13,name_tail = -5
#  )
#  names(list_result_melanoma) <- c("BayesSpace_DESpace", "StLearn_DESpace",
#                               "BayesSpace_scran", "StLearn_scran",
#                               "BayesSpace_seurat", "StLearn_seurat",
#                               "MERINGUE", "SPARK","nnSVG", "SPARK-X", "SpatialDE2","SpatialDE","SpaGCN",
#                               "All_results")
# #
#  save(list_result_melanoma, file = "./Real_data/results/melanoma/All_methods_genes_order.rda")

colors_method2 <- c(colors_method[1:10], "orchid",colors_method[11:12],"aquamarine",colors_method[13:14])
names(colors_method2)[c(11,14)] <- c("Manual_findMarkers", "Manual_FindAllMarkers")

load("./Real_data/results/melanoma/All_methods_genes_order.rda")

gg_Bar_melanoma <- ave_jaccard_btw_rep_melanoma(
  list_test_all=list_result_melanoma,
  path = path_save,
  sample_names=sample_names,
  top_genes = 200,
  col_method = colors_method2,
  data = 'melanoma')


d3_LIBD <- read.csv("./Real_data/results/LIBD/computation_time_Final.csv")

methods_all <- c("Manual_DESpace","BayesSpace_DESpace","StLearn_DESpace",
                 "MERINGUE", "nnSVG",
                 "SPARK","SPARK-X",
                 "SpatialDE", "SpatialDE2",
                 "SpaGCN",
                 "Manual_findMarkers","BayesSpace_findMarkers","StLearn_findMarkers",
                 "Manual_FindAllMarkers","BayesSpace_FindAllMarkers","StLearn_FindAllMarkers"
)
gg_legend <- ggplot(d3_LIBD,
                    #aes_string(y = "Time_Minutes", fill = "Method"),
                    aes(x = Method, y = Time_Minutes, fill = Method)) +  
  geom_bar( stat = "identity",position = position_dodge(0.8), width = 0.7) +  theme_bw() +   xlab("") +  ylab("Time (Minutes)") + 
  scale_y_sqrt( breaks = c(5,10,20,50,100,500,1000,2000) ) + 
  geom_text(aes(label=Time), vjust=-0.5,size = 3.5, fontface = "bold") +
  scale_fill_manual("Method", values =colors_method2, breaks = methods_all) +
  guides(colour = guide_legend(ncol = 5, byrow = FALSE,
                               override.aes = list(fill = colors_method2) )) +
  theme(axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1),   
        axis.text.y= element_text(size=rel(0.8),angle = 0, hjust = 1),  
        axis.title = element_text(size=rel(1)),      
        panel.grid.minor = element_blank(),         
        panel.grid.major.x = element_blank(),       
        panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),        
        aspect.ratio = 1,        
        legend.position = "bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=30),
        legend.key.width=unit(2, "cm"),
        legend.box="vertical",
        legend.key.height=unit(2.25,"line"),
        legend.margin=margin())

gg_Bar_LIBD2 <- gg_Bar_LIBD +
  theme(axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1),   
        axis.text.y= element_text(size=30,angle = 0, hjust = 1), 
  ) 
gg_Bar_melanoma2 <- gg_Bar_melanoma +
  theme(axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1),   
        axis.text.y= element_text(size=30,angle = 0, hjust = 1), 
  ) 
(AA = egg::ggarrange( plots = list(gg_Bar_LIBD2 +ylab("")+ labs(title = "LIBD") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                                   gg_Bar_melanoma2 +ylab("") +  labs(title = "melanoma") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size=35))
),
bottom = ggpubr::get_legend(gg_legend),
ncol = 2, nrow = 1))
ggsave(filename = paste0(path_save,"Barplot_Jaccrd.pdf"),
       plot = AA,
       device = "pdf",
       width = 25,
       height = 15,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
