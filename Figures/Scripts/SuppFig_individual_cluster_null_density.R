rm(list =ls())
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Function.R")
path_save = "./Figures/Figures/Supplementary/"
raw_path = "./Simulation/Output/"

library(DESpace)
data = "LIBD"
sample_names <- c("151507", "151669", "151673")
spatial_probs = c(0.5,0.9)
default = FALSE
#pattern = "mixture_patch"
pattern_list <- c('mixture_patch','mixture_reverse_patch')

Vectors <- expand.grid(sample_names,
                       pattern_list)
rank_match <- function(row){
  x = tab2$in_cluster[row]
  # Non-SV genes
  if(is.na(x)){
    rank = tab2[row,(13+num_cluster+1):(13+2*num_cluster)]
    return(names(rank[which(rank==1)]))
  }else{
    # SV genes -> x = real "in" cluster number
    return(tab2[[x]][row]) # -> the rank of the "x" cluster
  }
  
} 

`%notin%` <- Negate(`%in%`)

DF_NONSV_pval <- .GetPVal(Vectors = Vectors,
                     raw_path = raw_path,
                     data = data,
                     sample = sample,
                     pattern = pattern,
                     cluster_name = "stLearn",
                     OPTION = "without",
                     spatial_probs = c(0.5, 0.9))

DF <- do.call("rbind", DF_NONSV_pval)
DF$Method <- rep("StLearn_DESpace", dim(DF)[1])
require(dplyr)
DF <- 
  mutate(DF, Pattern = case_when(
  Pattern == 'mixture_patch' ~ factor("Mixture"), 
  Pattern == 'mixture_reverse_patch' ~ factor("Inverted mixture"), 
  TRUE ~ Pattern 
))

(AA1 <- ggplot(data =  DF, aes(x = Pvalue, y = ..ndensity.., lty = Cluster,
                             col = Method, fill = Method)) +
    geom_density(adjust = 0.5, size = 0.3, alpha = 0.1) +
    #facet_wrap(`Pattern` ~ `Sample`))
    facet_grid(Pattern ~ Sample) +
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
  scale_fill_manual(values = colors_method))


DF_NONSV_pval <- .GetPVal(Vectors = Vectors,
                          raw_path = raw_path,
                          data = data,
                          sample = sample,
                          pattern = pattern,
                          cluster_name = "BayesSpace",
                          OPTION = "without",
                          spatial_probs = c(0.5, 0.9))

DF <- do.call("rbind", DF_NONSV_pval)
DF$Method <- rep("BayesSpace_DESpace", dim(DF)[1])
require(dplyr)
DF <- 
  mutate(DF, Pattern = case_when(
    Pattern == 'mixture_patch' ~ factor("Mixture"), 
    Pattern == 'mixture_reverse_patch' ~ factor("Inverted mixture"), 
    TRUE ~ Pattern 
  ))

(AA2 <- ggplot(data =  DF, aes(x = Pvalue, y = ..ndensity.., lty = Cluster,
                               col = Method, fill = Method)) +
    geom_density(adjust = 0.5, size = 0.3, alpha = 0.1) +
    #facet_wrap(`Pattern` ~ `Sample`))
    facet_grid(Pattern ~ Sample) +
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
    scale_fill_manual(values = colors_method))

############################################# With recomputed dispersion ###########################################
DF_NONSV_pval <- .GetPVal(Vectors = Vectors,
                          raw_path = raw_path,
                          data = data,
                          sample = sample,
                          pattern = pattern,
                          cluster_name = "stLearn",
                          OPTION = "with",
                          spatial_probs = c(0.5, 0.9))

DF <- do.call("rbind", DF_NONSV_pval)
DF$Method <- rep("StLearn_DESpace", dim(DF)[1])
require(dplyr)
DF <- 
  mutate(DF, Pattern = case_when(
    Pattern == 'mixture_patch' ~ factor("Mixture"), 
    Pattern == 'mixture_reverse_patch' ~ factor("Inverted mixture"), 
    TRUE ~ Pattern 
  ))

(AA3 <- ggplot(data =  DF, aes(x = Pvalue, y = ..ndensity.., lty = Cluster,
                               col = Method, fill = Method)) +
    geom_density(adjust = 0.5, size = 0.3, alpha = 0.1) +
    #facet_wrap(`Pattern` ~ `Sample`))
    facet_grid(Pattern ~ Sample) +
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
    scale_fill_manual(values = colors_method))

DF_NONSV_pval <- .GetPVal(Vectors = Vectors,
                          raw_path = raw_path,
                          data = data,
                          sample = sample,
                          pattern = pattern,
                          cluster_name = "BayesSpace",
                          OPTION = "with",
                          spatial_probs = c(0.5, 0.9))

DF <- do.call("rbind", DF_NONSV_pval)
DF$Method <- rep("BayesSpace_DESpace", dim(DF)[1])
require(dplyr)
DF <- 
  mutate(DF, Pattern = case_when(
    Pattern == 'mixture_patch' ~ factor("Mixture"), 
    Pattern == 'mixture_reverse_patch' ~ factor("Inverted mixture"), 
    TRUE ~ Pattern 
  ))

(AA4 <- ggplot(data =  DF, aes(x = Pvalue, y = ..ndensity.., lty = Cluster,
                               col = Method, fill = Method)) +
    geom_density(adjust = 0.5, size = 0.3, alpha = 0.1) +
    #facet_wrap(`Pattern` ~ `Sample`))
    facet_grid(Pattern ~ Sample) +
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
    scale_fill_manual(values = colors_method))

################# Legend ######################
DF$Method[1:15126] <- c("StLearn_DESpace")
(legend_all <- ggplot(data =  DF, aes(x = Pvalue, y = ..ndensity.., 
                               col = Method, fill = Method)) +
    geom_density(adjust = 0.5, size = 0.3, alpha = 0.1) +
    #facet_wrap(`Pattern` ~ `Sample`))
    facet_grid(Pattern ~ Sample) +
    scale_alpha_manual(values = c("TRUE" = 0.1, "FALSE" = 0.4)) +
    guides(col = FALSE,
           fill = guide_legend(nrow = 1, order = 1,title = "",
                               override.aes = list(alpha = 1, col = NA))) +
    scale_x_continuous("p-value", breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
    scale_y_continuous("density", breaks = c(0, 1), expand = c(0, 0.1)) +
    .prettify("bw") + 
    my_theme +
    scale_colour_manual(values = colors_method) +
    scale_fill_manual(values = colors_method))

################# Save plots ##################
(AA = egg::ggarrange( plots = list(AA2 +ylab("")+ labs(title = "BayesSpace_DESpace") + theme(plot.title = element_text(hjust = 0, face = "bold", size=10),
                                                                                             legend.position = "none"),
                                   AA1 +  labs(title = "StLearn_DESpace") + theme(plot.title = element_text(hjust = 0, face = "bold", size=10),
                                                                                     legend.position = "none")
),
bottom = ggpubr::get_legend(legend_all),
ncol = 1, nrow = 2))
ggsave(filename = "IndividualCluster_null_density_LIBD_without.pdf",
       plot = AA,
       device = "pdf",
       path = path_save,
       width = 9,
       height = 6,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

(AA = egg::ggarrange( plots = list(AA4 + ylab("")+ labs(title = "BayesSpace_DESpace") + theme(plot.title = element_text(hjust = 0, face = "bold", size=10),
                                                                                             legend.position = "none"),
                                   AA3 +  labs(title = "StLearn_DESpace") + theme(plot.title = element_text(hjust = 0, face = "bold", size=10),
                                                                                  legend.position = "none")
),
bottom = ggpubr::get_legend(legend_all),
ncol = 1, nrow = 2))
ggsave(filename = "IndividualCluster_null_density_LIBD_with.pdf",
       plot = AA,
       device = "pdf",
       path = path_save,
       width = 9,
       height = 6,
       units = "in",
       dpi = 300,
       limitsize = TRUE)




.GetPVal <- function(Vectors = Vectors,
                     raw_path = raw_path,
                     data = data,
                     sample = sample,
                     pattern = pattern,
                     cluster_name = cluster_name,
                     OPTION = OPTION,
                     spatial_probs = c(0.5, 0.9)
                     ){
  DF_NONSV_pval <- list()
  for(v in 1:nrow(Vectors)){
    sample <- Vectors[v,1]
    pattern <- Vectors[v,2]
    #i <- ceiling(v/2) # the ith sample
    path = paste0(raw_path,"individual_SV_cluster/Version_use_BayesSpace_stLearn/", data,
                  "/",sample, "_", pattern, "_",cluster_name,"_results_",OPTION,".rda")
    ## Results
    load(path)
    if(OPTION == "without") results = results2
    results <- purrr::map(results, ~ rename_with(., ~ ifelse(
      grepl("gene_name", .), "gene_id", .
    )))
    
    ll <- lapply(seq_len(length(results)), function(i){
      data.frame(gene_id = results[[i]]$gene_id, 
                 value  = get("PValue", results[[i]]))})
    names(ll) <- names(results)
    
    ll2 <- lapply(ll, function(LL) lapply(LL, `[`, order(LL$gene_id)) )
    genesID <- ll2[[1]]$gene_id
    ll2 <-lapply(ll2, "[",  "value")
    
    com_results <- as.data.frame(do.call(cbind, lapply(ll2, as.data.frame)))
    rownames(com_results) <- genesID
    colnames(com_results) <- names(results)
    
    ## SV genes
    dir = paste0(raw_path,sample,'/',pattern,'/probs_',spatial_probs[1],'_',
                 spatial_probs[2],default,'/')
    genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
    colnames(genes) <- "gene_id"
    
    DF_NONSV <- com_results %>% filter(rownames(com_results) %notin% genes$gene_id)
    DF_NONSV_pval[[v]] <- DF_NONSV
    DF_NONSV_pval[[v]] <- reshape2::melt(DF_NONSV_pval[[v]]) # as rownames, i.e., gene names, are not indentical -> would be discarded
    colnames(DF_NONSV_pval[[v]]) <- c("Cluster", "Pvalue")
    #DF_NONSV_pval[[v]]$Layer <- gsub("_FDR.*$","", DF_NONSV_pval[[v]]$Layer)
    DF_NONSV_pval[[v]]$Cluster <- gsub("_PValue.*$","", DF_NONSV_pval[[v]]$Cluster)
    DF_NONSV_pval[[v]]$Sample <- paste("Sample", sample)
    DF_NONSV_pval[[v]]$Pattern <- pattern
  }
  return(DF_NONSV_pval)
}
