rm(list =ls())
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Functions.R")
library(qdap)
linetype_pattern <- c("solid", "dashed", "dotted", "dotdash", "longdash")
path_save = "./Figures/Figures/Supplementary/"
path = "./Simulation/Output/"
data = "LIBD"
sample_names <- c("151507", "151669","151673")
patterns = c("bottom_patch","circle_patch",
             "Manual_clusters_patch",
             "mixture_patch","mixture_reverse_patch")
all_results <- list()
spatial_probs = c(0.5,0.9);default = FALSE;rep_i = 1
k <- 1
setwd("./DESpace_data")
for(i in seq_along(patterns)){
  for(j in seq_along(sample_names)){
    dir = paste0(path,sample_names[j],'/',patterns[i],'/probs_',spatial_probs[1],'_',
                 spatial_probs[2],default,'/')
    file = paste0(dir,'results_all_revisions.csv')
    result <- read.csv(file,header = TRUE,sep =  ',')
    #genes <- read.csv(paste0(dir,'probs_',spatial_probs[2],'_selected_genes.txt'),header = TRUE,sep =  '\t')
    genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
    all_genes <- as.data.frame((unique(result$genes)))
    
    all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
    genes_id <- apply(genes,1,function(x) substr(x, 10,18))
    status <- ifelse(all_genes_id %in% genes_id, 0, 1)
    truth <- cbind(all_genes, status)
    # check
    sum(truth$status) == dim(genes)[1]
    if(sum(all_genes_id[length(all_genes_id)] == "_count")>=1){
      truth <- truth[1:dim(truth)[1]-1,]
    }
    result <- result %>% dplyr::distinct()
    setDT(result)
    result$NonSVGs <- truth$status[match(result$genes, truth$`(unique(result$genes))`)]
    # sum((result$NonSVGs == 1))/dim(result)[1] == 0.67 or 0.5 ## percentage of non-SVGs 
    result$method <- mgsub(methods_order,methods_all[-1],result$method)
    result <- result %>%
      filter(method != 'SpaGCN') %>%
      filter(NonSVGs == 1) %>%
      mutate(method =  factor(method, levels = methods_all)) %>%
      arrange(method)
    result$sample <- sample_names[j]
    result$pattern <- patterns[i]
    all_results[[k]] <- result
    k <- k + 1
  }
}
DF <- do.call(rbind, all_results)

(AA1 <- ggplot(data =  DF, aes(x = p.value, y = after_stat(ndensity), lty = interaction(sample, pattern),
                               col = method, fill = method)) +
    geom_density(adjust = 0.5, linewidth = 0.3, alpha = 0.1) +
    #facet_wrap(`Pattern` ~ `Sample`))
    facet_wrap(~method, ncol = 4) +
    scale_alpha_manual(values = c("TRUE" = 0.1, "FALSE" = 0.4)) +
    guides(col = FALSE,
           lty = guide_legend(ncol = 1, order = 2),
           fill = guide_legend(ncol = 2, order = 1,
                               override.aes = list(alpha = 1, col = NA))) +
    scale_x_continuous("p-value", breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
    scale_y_continuous("density", breaks = c(0, 1), expand = c(0, 0.1)) +
    .prettify("bw") + 
    my_theme +
    scale_colour_manual(values = colors_method) +
    scale_fill_manual(values = colors_method)) +
  scale_linetype_manual(values = rep(linetype_pattern, each = length(unique(DF$sample))))

###################################################
################## melanoma ###################
data = "melanoma"
all_results2 <- list()
sample_names = c("mel1_rep1", "mel2_rep1", "mel3_rep1", "mel4_rep1")
patterns = c("right_patch","circle_patch",
             "BayesSpace_clusters_patch",
             "mixture_patch","mixture_reverse_patch")
k=1
for(i in seq_along(patterns)){
  for(j in seq_along(sample_names)){
    dir = paste0(path,sample_names[j],'/',patterns[i],'/probs_',spatial_probs[1],'_',
                 spatial_probs[2],default,'/')
    file = paste0(dir,'results_all_revisions.csv')
    result <- read.csv(file,header = TRUE,sep =  ',')
    genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
    all_genes <- as.data.frame((unique(result$genes)))
    
    all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
    genes_id <- apply(genes,1,function(x) substr(x, 10,18))
    status <- ifelse(all_genes_id %in% genes_id, 0, 1)
    truth <- cbind(all_genes, status)
    # check
    sum(truth$status) == dim(genes)[1]
    if(sum(all_genes_id[length(all_genes_id)] == "_count")>=1){
      truth <- truth[1:dim(truth)[1]-1,]
    }
    result <- result %>% dplyr::distinct()
    setDT(result)
    result$NonSVGs <- truth$status[match(result$genes, truth$`(unique(result$genes))`)]
    # sum((result$NonSVGs == 1))/dim(result)[1] == 0.67 or 0.5 ## percentage of non-SVGs 
    result$method <- mgsub(methods_order,methods_all[-1],result$method)
    result <- result %>%
      filter(method != 'SpaGCN') %>%
      filter(NonSVGs == 1) %>%
      mutate(method =  factor(method, levels = methods_all)) %>%
      arrange(method)
    result$sample <- sample_names[j]
    result$pattern <- patterns[i]
    all_results2[[k]] <- result
    k <- k + 1
  }
}

DF2 <- do.call(rbind, all_results2)
(AA2 <- ggplot(data =  DF2, aes(x = p.value, y = after_stat(ndensity), lty = interaction(sample, pattern),
                                col = method, fill = method)) +
    geom_density(adjust = 0.5, linewidth = 0.3, alpha = 0.1) +
    #facet_wrap(`Pattern` ~ `Sample`))
    facet_wrap(~method, ncol = 4) +
    scale_alpha_manual(values = c("TRUE" = 0.1, "FALSE" = 0.4)) +
    guides(col = FALSE,
           lty = guide_legend(ncol = 1, order = 2),
           fill = guide_legend(ncol = 2, order = 1,
                               override.aes = list(alpha = 1, col = NA))) +
    scale_x_continuous("p-value", breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
    scale_y_continuous("density", breaks = c(0, 1), expand = c(0, 0.1)) +
    .prettify("bw") + 
    my_theme +
    scale_colour_manual(values = colors_method) +
    scale_fill_manual(values = colors_method)) +
  scale_linetype_manual(values = rep(linetype_pattern, each = length(unique(DF2$sample))))


###################################################
################## mouse cerebellum ###################
data = "Slideseq2"
all_results3 <- list()
patterns = c("bottom_patch","circle_patch",
             "StLearn_clusters_patch",
             "mixture_patch","mixture_reverse_patch")
k=1
for(i in seq_along(patterns)){
  dir = paste0(path,'SlideSeq2/',patterns[i],'/probs_',spatial_probs[1],'_',
               spatial_probs[2],default,'/')
  #file = paste0(dir,rep_i,'_probs_',spatial_probs[1],'_',
  # spatial_probs[2],'_results.txt')
  file = paste0(dir,'results_all_revisions.csv')
  result <- read.csv(file,header = TRUE,sep =  ',')
  #genes <- read.csv(paste0(dir,'probs_',spatial_probs[2],'_selected_genes.txt'),header = TRUE,sep =  '\t')
  genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
  all_genes <- as.data.frame((unique(result$genes)))
  
  all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
  genes_id <- apply(genes,1,function(x) substr(x, 10,18))
  status <- ifelse(all_genes_id %in% genes_id, 0, 1)
  truth <- cbind(all_genes, status)
  # check
  sum(truth$status) == dim(genes)[1]
  if(sum(all_genes_id[length(all_genes_id)] == "_count")>=1){
    truth <- truth[1:dim(truth)[1]-1,]
  }
  result <- result %>% dplyr::distinct()
  setDT(result)
  result$NonSVGs <- truth$status[match(result$genes, truth$`(unique(result$genes))`)]
  # sum((result$NonSVGs == 1))/dim(result)[1] == 0.67 or 0.5 ## percentage of non-SVGs 
  result$method <- mgsub(methods_order,methods_all[-1],result$method)
  result <- result %>%
    filter(method != 'SpaGCN') %>%
    filter(NonSVGs == 1) %>%
    mutate(method =  factor(method, levels = methods_all)) %>%
    arrange(method)
  result$pattern <- patterns[i]
  all_results3[[k]] <- result
  k <- k + 1
  
}

DF3 <- do.call(rbind, all_results3)
(AA3 <- ggplot(data =  DF3, aes(x = p.value, y = after_stat(ndensity), lty = pattern,
                                col = method, fill = method)) +
    geom_density(adjust = 0.5, linewidth = 0.3, alpha = 0.1) +
    #facet_wrap(`Pattern` ~ `Sample`))
    facet_wrap(~method, ncol = 4) +
    scale_alpha_manual(values = c("TRUE" = 0.1, "FALSE" = 0.4)) +
    guides(col = FALSE,
           lty = guide_legend(ncol = 1, order = 2),
           fill = guide_legend(ncol = 2, order = 1,
                               override.aes = list(alpha = 1, col = NA))) +
    scale_x_continuous("p-value", breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
    scale_y_continuous("density", breaks = c(0, 1), expand = c(0, 0.1)) +
    .prettify("bw") + 
    my_theme +
    scale_colour_manual(values = colors_method) +
    scale_fill_manual(values = colors_method)) +
  scale_linetype_manual(values = linetype_pattern)

##################################
############ legend ############
DF3$Pattern <- mgsub(patterns,c("Bottom/Right", "Circular", "Annotations", 
                                "Mixture", "Inverted mixture"),DF3$pattern)
(legend_all <- ggplot(data =  DF3, aes(x = p.value, y = after_stat(ndensity), lty =Pattern,
                                       col = method, fill = method)) +
    geom_density(adjust = 0.5, linewidth = 0.3, alpha = 0.1) +
    #facet_wrap(`Pattern` ~ `Sample`))
    facet_wrap(~method, ncol = 4) +
    scale_alpha_manual(values = c("TRUE" = 0.1, "FALSE" = 0.4)) +
    guides(col = FALSE,
           lty = guide_legend(ncol = 1, order = 2),
           fill = guide_legend(ncol = 2, order = 1,title = "Method",
                               override.aes = list(alpha = 1, col = NA))) +
    scale_x_continuous("p-value", breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
    scale_y_continuous("density", breaks = c(0, 1), expand = c(0, 0.1)) +
    .prettify("bw") + 
    my_theme + #theme(legend.direction = "vertical", legend.box = "vertical") +
    scale_colour_manual(values = colors_method) +
    scale_fill_manual(values = colors_method)) +
  scale_linetype_manual(values = linetype_pattern)


######################################################
############################# output #####################
AA = grid.arrange(ggpubr::ggarrange(AA1 + labs(title = "LIBD") + theme(plot.title = element_text(hjust = 0, face = "bold", size=10),
                                                                       legend.position = "none",aspect.ratio=0.3)),
                  ggpubr::ggarrange(AA2 +  labs(title = "melanoma") + theme(plot.title = element_text(hjust = 0, face = "bold", size=10),
                                                                            legend.position = "none",aspect.ratio=0.3)),
                  ggpubr::ggarrange(AA3 +  labs(title = "mouse cerebellum") + theme(plot.title = element_text(hjust = 0, face = "bold", size=10),
                                                                                    legend.position = "none",aspect.ratio=0.3)),
                  ggpubr::get_legend(legend_all),
                  heights = c(1,1,1,0.6),
                  ncol = 1, nrow = 4 )

ggsave(filename = "IndividualSample_null_density.pdf",
       plot = AA,
       device = "pdf",
       width = 10,
       height = 12,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
