## Figure 3: FDR - individual samples
rm(list = ls())
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Function.R")
path_save = "./Figures/Figures/Supplementary/"
raw_path = "./Simulation/Output/"
# remove the first "Manual_DESpace" color and method
colours = c(all_colours[-1], "white")
methods_all = methods_all[-1]
## For LIBD
gg_roc_LIBD = gg_fdr_LIBD = list()
sample_names = c("151507", "151669", "151673")
patterns = c("bottom_patch","circle_patch",
             "Manual_clusters_patch",
             "mixture_patch","mixture_reverse_patch")
k <- 1
gg <- list()
for(i in seq_along(patterns)){
  for(j in seq_along(sample_names)){
    if(i %in% c(4,5)) spatial_probs = c(0.6,0.6)
    if(i %in% c(1,2,3)) spatial_probs = c(0.6,0.9)
    if(i == 1 && j == 3){ ## 151673 with bottom_patch -> lose 'SPARK'
      gg <- NULL
    }else{
      gg = roc_fdr_plot(
        j,
        sample_names = sample_names, 
        pattern_names = patterns[i],
        path = path,colours = colours,path_save = path_save,
        spatial_probs = spatial_probs,default = FALSE,rep_i = 1,
        methods_all = methods_all,
        methods_order = methods_order,
        shape_border = shape_border,
        shape_fill = shape_fill,
        # roc_save_name='_ROC_final.pdf',
        # fdr_save_name = '_FDR_final.pdf',
        dataset = 'LIBD')
    }
    gg_roc_LIBD[[k]] = gg[[1]]
    gg_fdr_LIBD[[k]] = gg[[2]]
    print(paste0("i: ", i, "; j: ", j, "; k: ", k))
    k <- k + 1
  }  
}

## For 151673; bottom pattern: without 'SPARK'
i <- 3; pattern_names <- c("bottom_patch")
spatial_probs <- c(0.6, 0.9); default <- FALSE
methods_all <- methods_all[-3] # remove 'SPARK'
colours <- colours[-3]
dir = paste0(path,sample_names[i],'/',pattern_names,'/probs_',spatial_probs[1],'_',
             spatial_probs[2],default,'/')
#file = paste0(dir,rep_i,'_probs_',spatial_probs[1],'_',
# spatial_probs[2],'_results.txt')
file = paste0(dir,'results_all.csv')
result <- read.csv(file,header = TRUE,sep =  '\t')
#genes <- read.csv(paste0(dir,'probs_',spatial_probs[2],'_selected_genes.txt'),header = TRUE,sep =  '\t')
genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
all_genes <- as.data.frame((unique(result$genes)))

all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
genes_id <- apply(genes,1,function(x) substr(x, 10,18))
status <- ifelse(all_genes_id %in% genes_id, 1, 0)
truth <- cbind(all_genes, status)
# check
sum(truth$status) == dim(genes)[1]
if(sum(all_genes_id[length(all_genes_id)] == "_count")>=1){
  truth <- truth[1:dim(truth)[1]-1,]
}
result <- result %>% dplyr::distinct()
setDT(result)
data <- dcast(result, method ~ genes,value.var = 'adj.p.value', fun.aggregate=mean)
#data <- dcast(result, method ~ genes,value.var = 'adj.p.value')
df2 = truth[match(colnames(data)[-1], truth$`(unique(result$genes))`),]
df2 <- as.data.frame(df2[,2])
rownames(df2) <- colnames(data)[-1]
colnames(df2) <- 'status'
data2 <- dcast(result, method ~ genes,value.var = 'p.value', fun.aggregate=mean)
data <- data %>%
  mutate(method =  factor(method, levels = methods_order)) %>%
  arrange(method)
data2 <- data2 %>%
  mutate(method =  factor(method, levels = methods_order)) %>%
  arrange(method)

pval = data.frame(t(data2[, -1])
)
colnames(pval) <- methods_all
padj = data.frame(t(data[, -1])
)
colnames(padj) <- methods_all
DF_COBRA <- COBRAData(pval,
                      padj,
                      truth = data.frame(df2))
perf <- calculate_performance(DF_COBRA, binary_truth = "status",
                              thrs = c(0.01, 0.05, 0.1, 0.2))
cobra_plot <- prepare_data_for_plot(perf, colorscheme = colours, incloverall = FALSE,
                                    facetted = TRUE,conditionalfill = FALSE,
                                    keepmethods = c("BayesSpace_DESpace","StLearn_DESpace","SPARK","SPARK-X","SpatialDE","SpatialDE2","MERINGUE","SpaGCN","nnSVG"))
my.cols <- colours

(gg_roc_with_pval <- plot_roc(cobra_plot,linewidth=2) + scale_color_manual(values = my.cols,
                                                 name = "",
                                                 breaks=methods_all,
                                                 labels=methods_all) +
    #scale_fill_brewer(breaks=methods, type="qual", palette=3)+
    # scale_color_manual(values = cobra_plot@plotcolors, name = "", labels = methods_order,limits = force) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.position = "top") )
####
if(spatial_probs[1] == 0.6){gg_roc_with_pval <- gg_roc_with_pval + labs(y="TPR_strong", x = "TPR_weak")}

my_theme = theme(legend.position = "none", 
                 text = element_text(size = 9),
                 plot.title = element_text(hjust = 0, face = "bold", size=55),
                 plot.subtitle = element_text(hjust = 0.5, face = "bold", size=50))
legend <- ggpubr::get_legend(gg_roc_LIBD[[1]] +
                               guides(colour = guide_legend(nrow = 2, byrow = FALSE,
                                                            override.aes = list(shape = shape_fill,
                                                                                fill = colours[-length(colours)]) ) ) +
                               theme(legend.key.height=unit(3,"line"),
                                     legend.text=element_text(size=50),
                                     legend.key.width=unit(6,"line")))

ggpubr::as_ggplot(legend)

AA1 = ggpubr::ggarrange(gg_roc_LIBD[[1]] + xlab("") + labs(title = "Sample 151507", subtitle = "Bottom/Right") + my_theme,
                        gg_roc_LIBD[[2]] + xlab("") + ylab("") + labs(subtitle = "Circular") + my_theme,
                        gg_roc_LIBD[[1]] + xlab("") + ylab("") + labs(subtitle = "Annotations") + my_theme,
                        gg_roc_LIBD[[4]] + xlab("") + ylab("") + labs(subtitle = "Mixture") + my_theme,
                        gg_roc_LIBD[[5]] + xlab("") + ylab("") + labs(subtitle = "Inverted mixture") + my_theme,
                        
                        gg_roc_LIBD[[6]] + xlab("") + labs(title = "Sample 151507") + my_theme,
                        gg_roc_LIBD[[7]] + xlab("") + ylab("") + my_theme,
                        gg_roc_LIBD[[8]] + xlab("") + ylab("") + my_theme,
                        gg_roc_LIBD[[9]] + xlab("") + ylab("") + my_theme,
                        gg_roc_LIBD[[10]] + xlab("") + ylab("") + my_theme,
                        
                        gg_roc_LIBD[[11]] + xlab("") + labs(title = "Sample 151507") + my_theme,
                        gg_roc_LIBD[[12]] + xlab("") + ylab("") + my_theme,
                        gg_roc_LIBD[[13]] + xlab("") + ylab("") + my_theme,
                        gg_roc_LIBD[[14]] + xlab("") + ylab("") + my_theme,
                        gg_roc_LIBD[[15]] + xlab("") + ylab("") + my_theme,
                        
                        ncol = 5, nrow = 3)
AA = ggpubr::ggarrange( AA1, legend, 
                        # font.label = list(color = "black", size = 40, face = "bold"),
                        #label.x = 0, labrl.y = 1,
                        #labels = c("LIBD", "melanoma", "mouse cerebellum"),
                        heights=c(9,1),
                        ncol = 1, nrow = 2 )

spatial_probs = c(0.6,0.9)
ggsave(filename = paste0('Weak_vs_strong_',spatial_probs[1],'_',
                         spatial_probs[2],"_ROC_LIBD.pdf"),
       plot = AA,
       device = "pdf",
       path = path_save,
       width = 45,
       height = 35,
       units = "in",
       dpi = 300,
       limitsize = TRUE)


## For melanoma
rm(list = ls())
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Function.R")
path_save = "./Figures/Figures/Supplementary/"
raw_path = "./Simulation/Output/"
# remove the first "Manual_DESpace" color and method
colours = c(all_colours[-1], "white")
methods_all = methods_all[-1]

gg_roc_melanoma = gg_fdr_melanoma = list()
sample_names = c("mel1_rep1", "mel2_rep1", "mel3_rep1", "mel4_rep1")
patterns = c("right_patch","circle_patch",
             "BayesSpace_clusters_patch",
             "mixture_patch","mixture_reverse_patch")
k=1
for(i in seq_along(patterns)){
  for(j in seq_along(sample_names)){
    gg = roc_fdr_plot(
      j,
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
    gg_roc_melanoma[[k]] = gg[[1]]
    gg_fdr_melanoma[[k]] = gg[[2]]
    k <- k + 1
  }  
}
my_theme = theme(legend.position = "none", 
                 text = element_text(size = 9),
                 axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(4)),
                 axis.text.y=element_text(size=rel(4)),
                 axis.title.y = element_text(size=rel(4)),
                 axis.title.x = element_text(size=rel(4)),
                 legend.title=element_text(size=rel(2)),
                 legend.text=element_text(size=rel(4)),
                 plot.title = element_text(hjust = 0, face = "bold", size=35),
                 plot.subtitle = element_text(hjust = 0.5, face = "bold", size=30))
legend <- ggpubr::get_legend(gg_roc_melanoma[[1]] +
                               guides(colour = guide_legend(nrow = 2, byrow = FALSE,
                                                            override.aes = list(shape = shape_fill,size = 2.5,
                                                                                fill = colours[-length(colours)]) ))+#,
                                      #linetype = guide_legend(override.aes = list(size = 2))) +
                               theme(legend.key.height=unit(3,"line"),
                                     legend.text=element_text(size=30),
                                     legend.key.width=unit(6,"line")))

ggpubr::as_ggplot(legend)
AA1 = ggpubr::ggarrange(gg_roc_melanoma[[1]] + xlab("") + labs(title = "Sample mel1_rep1", subtitle = "Bottom/Right") + my_theme,
                        gg_roc_melanoma[[2]] + xlab("") + ylab("") + labs(subtitle = "Circular") + my_theme,
                        gg_roc_melanoma[[3]] + xlab("") + ylab("") + labs(subtitle = "Annotations") + my_theme,
                        gg_roc_melanoma[[4]] + xlab("") + ylab("") + labs(subtitle = "Mixture") + my_theme,
                        gg_roc_melanoma[[5]] + xlab("") + ylab("") + labs(subtitle = "Inverted mixture") + my_theme,
                        
                        gg_roc_melanoma[[6]] + xlab("") + labs(title = "Sample mel2_rep1") + my_theme,
                        gg_roc_melanoma[[7]] + xlab("") + ylab("") + my_theme,
                        gg_roc_melanoma[[8]] + xlab("") + ylab("") + my_theme,
                        gg_roc_melanoma[[9]] + xlab("") + ylab("") + my_theme,
                        gg_roc_melanoma[[10]] + xlab("") + ylab("") + my_theme,
                        
                        gg_roc_melanoma[[11]] + xlab("") + labs(title = "Sample mel3_rep1") + my_theme,
                        gg_roc_melanoma[[12]] + xlab("") + ylab("") + my_theme,
                        gg_roc_melanoma[[13]] + xlab("") + ylab("") + my_theme,
                        gg_roc_melanoma[[14]] + xlab("") + ylab("") + my_theme,
                        gg_roc_melanoma[[15]] + xlab("") + ylab("") + my_theme,
                        
                        gg_roc_melanoma[[16]] + xlab("") + labs(title = "Sample mel4_rep1") + my_theme,
                        gg_roc_melanoma[[17]] + xlab("") + ylab("") + my_theme,
                        gg_roc_melanoma[[18]] + xlab("") + ylab("") + my_theme,
                        gg_roc_melanoma[[19]] + xlab("") + ylab("") + my_theme,
                        gg_roc_melanoma[[20]] + xlab("") + ylab("") + my_theme,
                        
                        ncol = 5, nrow = 4)
AA = ggpubr::ggarrange( AA1, legend, 
                        # font.label = list(color = "black", size = 40, face = "bold"),
                        #label.x = 0, labrl.y = 1,
                        #labels = c("LIBD", "melanoma", "mouse cerebellum"),
                        heights=c(12,1),
                        ncol = 1, nrow = 2 )

spatial_probs = c(0.5,0.9)
ggsave(filename = paste0('Weak_vs_strong_',spatial_probs[1],'_',
                         spatial_probs[2],"_ROC_melanoma.pdf"),
       plot = AA,
       device = "pdf",
       path = path_save,
       width = 40,
       height = 30,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

