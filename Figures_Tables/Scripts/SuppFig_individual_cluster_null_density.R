rm(list =ls())
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Functions.R")
path_save = "./Figures/Figures/Supplementary/"
colors_method <- c("plum2","mediumpurple",
                   "turquoise","turquoise4"
)
names(colors_method) <- c("BayesSpace_findMarkers", "StLearn_findMarkers",
                          "BayesSpace_FindAllMarkers", "StLearn_FindAllMarkers")
#library(DESpace)
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
    rank = tab2[row,(11+2*num_cluster+1):dim(tab2)[2]]
    return(names(rank[which(rank==1)]))
  }else{
    # SV genes -> x = real "in" cluster number
    return(tab2[[x]][row]) # -> the rank of the "x" cluster
  }
  
}
.top_results_scran <- function(results, select = "p.value"){
  cluster_results <- lapply(results, as.data.frame)
  cluster_results <- lapply(cluster_results, function(x) x[order(rownames(x)),])
  ll <- lapply(1:length(cluster_results), function(i){
    data.frame(gene_id = rownames(cluster_results[[i]]), 
               value  = get(select, cluster_results[[i]]))})
  names(ll) <- names(cluster_results)
  
  ll2 <- lapply(ll, function(LL) lapply(LL, `[`, order(LL$gene_id)) )
  genesID <- ll2[[1]]$gene_id
  ll2 <-lapply(ll2, "[",  "value")
  
  com_results <- as.data.frame(do.call(cbind, lapply(ll2, as.data.frame)))
  rownames(com_results) <- genesID
  colnames(com_results) <- names(cluster_results)
  # checked lapply(com_results, is.numeric) are TRUE
  sel_layer_min <- apply(com_results,1,which.min)
  rank_results<-as.data.frame(colnames(com_results)[sel_layer_min])
  colnames(rank_results) <- "Layer"
  rank_results <- cbind(rank_results, 
                        lapply(1:nrow(com_results), function(x) cluster_results[[sel_layer_min[x]]][x,]) %>% 
                          bind_rows()) 
  colnames(rank_results)[3] <- paste0("Layer_", colnames(rank_results)[3])
  colnames(com_results) <- paste0(names(cluster_results), "_", select)
  com_results['gene_id'] <- rownames(com_results)
  rank_results['gene_id'] <- rownames(rank_results)
  rownames(rank_results) <- rank_results$gene_id
  com_results <- merge(rank_results, com_results, by = "gene_id") 
  rownames(com_results) <- com_results$gene_id
  # Sort results based on FDR
  com_results <- as.data.frame(com_results) %>% arrange(Layer_FDR)
  return(com_results)
}
.GetPVal_scran <- function(Vectors = Vectors,
               data = data,
               sample = sample,
               pattern = pattern,
               cluster_name = "BayesSpace",
               methodMarker = "scran",
               spatial_probs = c(0.5, 0.9)){
  DF_NONSV_pval <- list()
  for(v in 1:nrow(Vectors)){
    sample <- Vectors[v,1]
    pattern <- Vectors[v,2]
    #i <- ceiling(v/2) # the ith sample
    path = paste0("./individual_SV_cluster/", data,
                  "/",sample, "_", pattern, "_", cluster_name, "_results_",methodMarker, ".rda")
    ## Results
    load(path)
    select = "p.value"
    FDR_results1 = .top_results_scran(results = results, select = select)
    (num_cluster = (dim(FDR_results1)[2]-5)/2)
    ## SV genes
    dir = paste0("./simulation/results/",sample,'/',pattern,'/probs_',spatial_probs[1],'_',
                 spatial_probs[2],default,'/')
    genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
    colnames(genes) <- "gene_id"
    
    DF_NONSV <- FDR_results1 %>% filter(gene_id %notin% genes$gene_id)
    DF_NONSV_pval[[v]] <- DF_NONSV[,(5+num_cluster+1):(dim(DF_NONSV)[2])]
    DF_NONSV_pval[[v]] <- reshape2::melt(DF_NONSV_pval[[v]]) # as rownames, i.e., gene names, are not indentical -> would be discarded
    colnames(DF_NONSV_pval[[v]]) <- c("Cluster", "Pvalue")
    DF_NONSV_pval[[v]]$Cluster <- gsub("_PValue.*$","", DF_NONSV_pval[[v]]$Cluster)
    DF_NONSV_pval[[v]]$Sample <- paste("Sample", sample)
    DF_NONSV_pval[[v]]$Pattern <- pattern
  }
  return(DF_NONSV_pval)
}
`%notin%` <- Negate(`%in%`)

DF_NONSV_pval <- .GetPVal_scran(Vectors = Vectors,
                          data = data,
                          sample = sample,
                          pattern = pattern,
                          cluster_name = "BayesSpace",
                          methodMarker = "scran",
                          spatial_probs = c(0.5, 0.9))

DF <- do.call("rbind", DF_NONSV_pval)
DF$Method <- rep("BayesSpace_findMarkers", dim(DF)[1])
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


DF_NONSV_pval <- .GetPVal_scran(Vectors = Vectors,
                                data = data,
                                sample = sample,
                                pattern = pattern,
                                cluster_name = "StLearn",
                                methodMarker = "scran",
                                spatial_probs = c(0.5, 0.9))

DF <- do.call("rbind", DF_NONSV_pval)
DF$Method <- rep("StLearn_findMarkers", dim(DF)[1])
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

################# Legend ######################
DF$Method[1:20] <- c("BayesSpace_findMarkers")
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
(AA = egg::ggarrange( plots = list(AA1 +ylab("")+ labs(title = "BayesSpace_findMarkers") + theme(plot.title = element_text(hjust = 0, face = "bold", size=10),
                                                                                           legend.position = "none"),
                                   AA2 +  labs(title = "StLearn_findMarkers") + theme(plot.title = element_text(hjust = 0, face = "bold", size=10),
                                                                                legend.position = "none")
),
bottom = ggpubr::get_legend(legend_all),
ncol = 1, nrow = 2))
ggsave(filename = "IndividualCluster_null_density_LIBD_scran.pdf",
       plot = AA,
       device = "pdf",
       width = 9,
       height = 6,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

############################################# With Seurat ###########################################
library(stats)
rank_match <- function(row){
  x = tab2$in_cluster[row]
  # Non-SV genes
  if(is.na(x)){
    rank = tab2[row,(11+2*num_cluster+1):dim(tab2)[2]]
    return(names(rank[which(rank==1)]))
  }else{
    # SV genes -> x = real "in" cluster number
    return(tab2[[x]][row]) # -> the rank of the "x" cluster
  }
  
}
# select = "p_val_adj"
.top_results_seurat <- function(results, select = "p_val"){
  cluster_results <- reshape(results, idvar = "gene", timevar = "cluster", direction = "wide")
  colnames(cluster_results)[1] <- "gene_id"
  rownames(cluster_results) <- cluster_results$genesID
  col_ind <- grepl( paste0(select,"."), colnames(cluster_results), fixed = TRUE)
  colnames(cluster_results)[col_ind] <- gsub(paste0(select,"."), "", colnames(cluster_results)[col_ind])
  com_results <- cluster_results[,col_ind]
  sel_layer_min <- apply(com_results,1,which.min)
  com_results['gene_id'] <- cluster_results$gene_id
  rank_results<-as.data.frame(colnames(com_results)[sel_layer_min])
  colnames(rank_results) <- "Layer"
  rank_results <- cbind(rank_results, 
                        lapply(1:nrow(com_results), function(x) {
                          top_cluster <- rank_results[x,]
                          col_ind <- grepl( paste0(".", top_cluster), colnames(cluster_results), fixed = TRUE)
                          res <- cbind(cluster_results[x, top_cluster], cluster_results[x,col_ind])
                          colnames(res)[1] <- paste0(select,".",top_cluster)
                          res <- res[,c(1,2,5)] # only keep p_val, avg_log2FC and p_val_adj
                          colnames(res) <- gsub(paste0(".", top_cluster), "", colnames(res))
                          res}) %>% 
                          bind_rows())[,1:4] # only keep rank, p_val, avg_log2FC and p_val_adj
  colnames(rank_results)[2:4] <- paste0("Layer_", colnames(rank_results)[2:4])
  rank_results['gene_id'] <- com_results$gene_id
  rownames(rank_results) <- rank_results$gene_id
  com_results <- merge(rank_results, com_results, by = "gene_id") 
  rownames(com_results) <- com_results$gene_id
  # Sort results based on FDR
  com_results <- as.data.frame(com_results) %>% arrange(Layer_p_val_adj)
  return(com_results)
}
.GetPVal_seurat <- function(Vectors = Vectors,
                           data = data,
                           sample = sample,
                           pattern = pattern,
                           cluster_name = "BayesSpace",
                           methodMarker = "seurat",
                           spatial_probs = c(0.5, 0.9)){
  DF_NONSV_pval <- list()
  for(v in 1:nrow(Vectors)){
    sample <- Vectors[v,1]
    pattern <- Vectors[v,2]
    #i <- ceiling(v/2) # the ith sample
    path = paste0("./individual_SV_cluster/", data,
                  "/",sample, "_", pattern, "_", cluster_name, "_results_",methodMarker, ".rda")
    ## Results
    load(path)
    select = "p_val"
    FDR_results1 = .top_results_seurat(results = results, select = select)
    (num_cluster = dim(FDR_results1)[2]-5)
    ## SV genes
    dir = paste0("./simulation/results/",sample,'/',pattern,'/probs_',spatial_probs[1],'_',
                 spatial_probs[2],default,'/')
    genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
    colnames(genes) <- "gene_id"
    
    DF_NONSV <- FDR_results1 %>% filter(gene_id %notin% genes$gene_id)
    DF_NONSV_pval[[v]] <- DF_NONSV[,(5+1):(dim(DF_NONSV)[2])]
    DF_NONSV_pval[[v]] <- reshape2::melt(DF_NONSV_pval[[v]]) # as rownames, i.e., gene names, are not indentical -> would be discarded
    colnames(DF_NONSV_pval[[v]]) <- c("Cluster", "Pvalue")
    DF_NONSV_pval[[v]]$Cluster <- gsub("_PValue.*$","", DF_NONSV_pval[[v]]$Cluster)
    DF_NONSV_pval[[v]]$Sample <- paste("Sample", sample)
    DF_NONSV_pval[[v]]$Pattern <- pattern
  }
  return(DF_NONSV_pval)
}

DF_NONSV_pval <- .GetPVal_seurat(Vectors = Vectors,
                          data = data,
                          sample = sample,
                          pattern = pattern,
                          cluster_name = "StLearn",
                          methodMarker = "seurat",
                          spatial_probs = c(0.5, 0.9))

DF <- do.call("rbind", DF_NONSV_pval)
DF$Method <- rep("StLearn_FindAllMarkers", dim(DF)[1])
require(dplyr)
DF <- 
  mutate(DF, Pattern = case_when(
    Pattern == 'mixture_patch' ~ factor("Mixture"), 
    Pattern == 'mixture_reverse_patch' ~ factor("Inverted mixture"), 
    TRUE ~ Pattern 
  ))
library(ggplot2)
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

DF_NONSV_pval <- .GetPVal_seurat(Vectors = Vectors,
                          data = data,
                          sample = sample,
                          pattern = pattern,
                          cluster_name = "BayesSpace",
                          methodMarker = "seurat",
                          spatial_probs = c(0.5, 0.9))

DF <- do.call("rbind", DF_NONSV_pval)
DF$Method <- rep("BayesSpace_FindAllMarkers", dim(DF)[1])
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
DF$Method[1:20] <- c("StLearn_FindAllMarkers")
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

(AA = egg::ggarrange( plots = list(AA4 + ylab("")+ labs(title = "BayesSpace_FindAllMarkers") + theme(plot.title = element_text(hjust = 0, face = "bold", size=10),
                                                                                              legend.position = "none"),
                                   AA3 +  labs(title = "StLearn_FindAllMarkers") + theme(plot.title = element_text(hjust = 0, face = "bold", size=10),
                                                                                  legend.position = "none")
),
bottom = ggpubr::get_legend(legend_all),
ncol = 1, nrow = 2))
ggsave(filename = "IndividualCluster_null_density_LIBD_Seurat.pdf",
       plot = AA,
       device = "pdf",
       width = 9,
       height = 6,
       units = "in",
       dpi = 300,
       limitsize = TRUE)




