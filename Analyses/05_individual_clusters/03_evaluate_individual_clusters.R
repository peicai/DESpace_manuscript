rm(list =ls())
source("./05_individual_clusters/01_Plots_function.R")
source("./05_individual_clusters/01_Individual_clusters_source_code.R")

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
wrapped_func <- function(v, OPTION, data, sample_names,
                         spatial_probs, default, pattern_list, DF_pval,
                         DF_rank, cluster_method){
  #for(v in 1:nrow(Vectors)){
  sample <- Vectors[v,1]
  pattern <- Vectors[v,2]
  path = paste0("./Data/Simulation/Ouput/individual_SV_cluster/", data,
                "/",sample, "_", pattern, "_", cluster_method, "_results_",OPTION,".rda")
  ## Results
  load(path)
  if(OPTION == "without") results = results2
  # library(dplyr)
  # library(purrr)
  
  results <- map(results, ~ rename_with(., ~ ifelse(
    grepl("gene_name", .), "gene_id", .
  )))
  
  FDR_results1 = .top_results(cluster_results = results, select = "FDR")
  (num_cluster = dim(FDR_results1)[2]-7)
  all_genes <- as.data.frame((unique(FDR_results1$gene_id)))
  ## Get index positions of all minima -> column name "Layer"
  #which(FDR_results1[,8:(7+num_cluster)] == min(FDR_results1[,8:(7+num_cluster)]))
  #FDR_results1$in_cluster <- apply(FDR_results1[,8:(7+num_cluster)], 1, function(x) which(x == min(x)))
  
  ## SV genes
  dir = paste0("./Data/Simulation/Output/",sample,'/',pattern,'/probs_',spatial_probs[1],'_',
               spatial_probs[2],default,'/')
  genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
  colnames(genes) <- "gene_id"
  
  ## Select SV genes only:
  tab = as.data.frame( merge(FDR_results1, genes, by="gene_id"))
  # ind <- ifelse(pattern == "mixture_reverse_patch",7,8)
  # if(pattern == "mixture_patch" && sample == "151669") {ind <- 7}
  ind <- 7
  # as.numeric(as.factor(tab$Layer)) + 1 == FDR_results1$in_cluster
  #IN = as.numeric(tab[ cbind(1:nrow(tab), as.numeric(as.factor(tab$Layer)) + 1 + 7 ) ]) # as.numeric("layer5") = 4, then plus 1 to get 5; start from the 7th column
  IN = as.numeric(tab[ cbind(1:nrow(tab), as.numeric(as.factor(tab$Layer))  + ind ) ]) # as.numeric("layer5") = 4, then plus 1 to get 5; start from the 7th column
  
  head(IN); length(IN)
  ## remove the in cluster column
  # x[7+num_cluster+1] -> 'in_cluster'
  # 1:7 -> gene_name  Layer Layer_LR Layer_logCPM Layer_logFC  Layer_PValue  Layer_FDR
  # +7 -> remove first 7th columns
  OUT = as.numeric(c(apply(tab, 1, function(x){
    id = as.numeric(as.factor(x[2]))
    x[-c(1:7,id+7)]
  })))
  DF_pval[[v]] = data.frame(p_values = c(IN, OUT),
                            top_SV = c(rep(TRUE, length(IN)),
                                       rep(FALSE, length(OUT)) ),
                            pattern = pattern,
                            sample = sample)
  
  # make wider
  # ranking:
  # For each gene
  tab2 = data.frame(gene_id = FDR_results1$gene_id,
                    pattern = pattern,
                    sample = sample,
                    cluster_method)
  tab2$SV <- ifelse(tab2$gene_id %in% genes$gene_id, TRUE, FALSE)
  
  ## Ground truth (rownames): pattern2 -> layer3; pattern3 -> layer4 and so on; pattern6 -> WM
  pat <- "^pattern_.*?([0-9]).*"
  genes$GoundTruth <-  unlist(lapply(rownames(genes),function(x) gsub(pat, "\\1", x)))
  genes$GoundTruth <- paste0("Layer",as.numeric(genes$GoundTruth)+1)
  genes$GoundTruth[genes$GoundTruth == "Layer7"] <- "WM"
  tab2<- merge(tab2, genes, all=TRUE, by = "gene_id")
  in_cluster <- tab2$GoundTruth
  ## in_cluster: Match the real "in" cluster to results
  ## Based on plots, manually match the cluster number
  if(sample == "151507"){
    if(pattern == "mixture_patch"){
      if(cluster_method == "BayesSpace"){
        in_cluster[in_cluster == "Layer3"] <- 4
        in_cluster[in_cluster == "Layer4"] <- 7
        in_cluster[in_cluster == "Layer5"] <- 2
        in_cluster[in_cluster == "Layer6"] <- 6
        in_cluster[in_cluster == "WM"] <- 3
      } # end if BayesSpce
      if(cluster_method == "stLearn"){
        in_cluster[in_cluster == "Layer3"] <- 1
        in_cluster[in_cluster == "Layer4"] <- 4
        in_cluster[in_cluster == "Layer5"] <- 2
        in_cluster[in_cluster == "Layer6"] <- 0
        in_cluster[in_cluster == "WM"] <- 3
      } # end if stLearn
    } # end if mixture pattern
    if(pattern == "mixture_reverse_patch"){
      if(cluster_method == "BayesSpace"){
        in_cluster[in_cluster == "Layer3"] <- 4
        in_cluster[in_cluster == "Layer4"] <- 6
        in_cluster[in_cluster == "Layer5"] <- 2
        in_cluster[in_cluster == "Layer6"] <- 5
        in_cluster[in_cluster == "WM"] <- 7
      } # end if BayesSpce
      if(cluster_method == "stLearn"){
        in_cluster[in_cluster == "Layer3"] <- 2#1
        in_cluster[in_cluster == "Layer4"] <- 1#4
        in_cluster[in_cluster == "Layer5"] <- 4#2
        in_cluster[in_cluster == "Layer6"] <- 0
        in_cluster[in_cluster == "WM"] <- 3
      } # end if stLearn
    } # end if mixture reverse pattern
  } # end if 151507
  if(sample == "151669"){
    if(pattern == "mixture_patch"){
      if(cluster_method == "BayesSpace"){
        in_cluster[in_cluster == "Layer3"] <- 2
        in_cluster[in_cluster == "Layer4"] <- 5
        in_cluster[in_cluster == "Layer5"] <- 1
        in_cluster[in_cluster == "Layer6"] <- 6
        in_cluster[in_cluster == "WM"] <- 4
      } # end if BayesSpce
      if(cluster_method == "stLearn"){
        in_cluster[in_cluster == "Layer3"] <- 4
        in_cluster[in_cluster == "Layer4"] <- 2
        in_cluster[in_cluster == "Layer5"] <- 1
        in_cluster[in_cluster == "Layer6"] <- 3
        in_cluster[in_cluster == "WM"] <- 0
      } # end if stLearn
    } # end if mixture pattern
    if(pattern == "mixture_reverse_patch"){
      if(cluster_method == "BayesSpace"){
        in_cluster[in_cluster == "Layer3"] <- 2
        in_cluster[in_cluster == "Layer4"] <- 5
        in_cluster[in_cluster == "Layer5"] <- 1
        in_cluster[in_cluster == "Layer6"] <- 6
        in_cluster[in_cluster == "WM"] <- 4
      } # end if BayesSpce
      if(cluster_method == "stLearn"){
        in_cluster[in_cluster == "Layer3"] <- 4
        in_cluster[in_cluster == "Layer4"] <- 0
        in_cluster[in_cluster == "Layer5"] <- 1
        in_cluster[in_cluster == "Layer6"] <- 2
        in_cluster[in_cluster == "WM"] <- 3
      } # end if stLearn
    } # end if mixture reverse pattern
  } # end if 151669
  if(sample == "151673"){
    if(pattern == "mixture_patch"){
      if(cluster_method == "BayesSpace"){
        in_cluster[in_cluster == "Layer3"] <- 4
        in_cluster[in_cluster == "Layer4"] <- 6
        in_cluster[in_cluster == "Layer5"] <- 5
        in_cluster[in_cluster == "Layer6"] <- 1
        in_cluster[in_cluster == "WM"] <- 2
      } # end if BayesSpce
      if(cluster_method == "stLearn"){
        in_cluster[in_cluster == "Layer3"] <- 3
        in_cluster[in_cluster == "Layer4"] <- 1
        in_cluster[in_cluster == "Layer5"] <- 4
        in_cluster[in_cluster == "Layer6"] <- 2
        in_cluster[in_cluster == "WM"] <- 0
      } # end if stLearn
    } # end if mixture pattern
    if(pattern == "mixture_reverse_patch"){
      if(cluster_method == "BayesSpace"){
        in_cluster[in_cluster == "Layer3"] <- 5
        in_cluster[in_cluster == "Layer4"] <- 4
        in_cluster[in_cluster == "Layer5"] <- 6
        in_cluster[in_cluster == "Layer6"] <- 1
        in_cluster[in_cluster == "WM"] <- 2
      } # end if BayesSpce
      if(cluster_method == "stLearn"){
        in_cluster[in_cluster == "Layer3"] <- 2
        in_cluster[in_cluster == "Layer4"] <- 0
        in_cluster[in_cluster == "Layer5"] <- 3
        in_cluster[in_cluster == "Layer6"] <- 1
        in_cluster[in_cluster == "WM"] <- 4
      } # end if stLearn
    } # end if mixture reverse pattern
  } # end if 151669
  tab2$in_cluster <- in_cluster
  tab2 <- merge(tab2, FDR_results1, by = "gene_id")
  # for each gene, rank based on p values; 1 -> smallest p values
  ord = as.data.frame(t(apply(tab2[,-c(1:13,13+num_cluster+1)], 1, function(x){rank(x)})))
  pat <- "^([0-9]).*"
  colnames(ord) <- gsub(pat, "\\1", colnames(ord))
  tab2 <- cbind(tab2, ord)
  rank_in <- apply(as.data.frame(1:dim(tab2)[1]), 1, rank_match)
  tab2$rank <- rank_in
  
  ## Only for SV genes
  tab2 <- tab2 %>% filter(SV == TRUE)
  IN = as.numeric(tab2$rank)
  head(IN); length(IN)
  OUT = as.numeric(c(apply(tab2, 1, function(x){
    in_clus = x[7] # tab2$in_cluster
    id = which(names(x) == as.character(in_clus))
    x[-c(1:(13+num_cluster), id, dim(tab2)[2])]}))) # only keep ranking columns except the top ranked cluster (remove first 13+num_cluster columns and the last column ("rank"))
  
  DF_rank[[v]] = data.frame(rank = c(IN, OUT),
                            top_SV = c(rep(TRUE, length(IN)),
                                       rep(FALSE, length(OUT)) ),
                            pattern = pattern,
                            sample = sample,
                            cluster_method)
  return(list(DF_pval, DF_rank))
}

# end function rank_match

########################################################################################################
#################################### without re-compute dispersion (stLearn) ##########################################
DF_pval <- list()
DF_rank <- list()
OPTION = "with" ## with dispersion results
cluster_method = "stLearn"
for(v in 1:nrow(Vectors)){
  wrapped_func(v, OPTION, data, sample_names, spatial_probs, default, pattern_list, DF_pval, DF_rank, cluster_method)
}
DF_pval_all <- do.call(rbind, DF_pval)
DF_rank_all <- do.call(rbind, DF_rank)

save(DF_pval_all, file = paste0('./Data/Simulation/Ouput/individual_cluster/LIBD_',cluster_method, '_pval_results_',OPTION,'.rda'))
save(DF_rank_all, file = paste0('./Data/Simulation/Ouput/individual_cluster/LIBD_',cluster_method, '_rank_results_',OPTION,'.rda'))

DF_rank_SV <- DF_rank_all %>% filter(top_SV)
mean(DF_rank_SV$rank == 1)

(table <- aggregate(DF_rank_SV$rank == 1, list(DF_rank_SV$sample, DF_rank_SV$pattern), FUN=mean))
colnames(table) <- c("Sample","Pattern","Accuracy")
write.table(table, sep = ",", file = paste0('./Data/Simulation/Ouput/individual_cluster/LIBD_',cluster_method, '_rank_table_score_',OPTION,'.csv'))

DF_NONSV_pval <- list()
cluster_name = "stLearn"
for(v in 1:nrow(Vectors)){
  sample <- Vectors[v,1]
  pattern <- Vectors[v,2]
  #i <- ceiling(v/2) # the ith sample
  path = paste0("./Data/Simulation/Ouput/individual_SV_cluster/", data,
                "/",sample, "_", pattern, "_",cluster_name,"_results_",OPTION,".rda")
  ## Results
  load(path)
  if(OPTION == "without") results = results2
  results <- map(results, ~ rename_with(., ~ ifelse(
    grepl("gene_name", .), "gene_id", .
  )))
  FDR_results1 = .top_results(cluster_results = results, select = "p_val")
  (num_cluster = dim(FDR_results1)[2]-7)
  ## SV genes
  dir = paste0("./Data/Simulation/Output/",sample,'/',pattern,'/probs_',spatial_probs[1],'_',
               spatial_probs[2],default,'/')
  genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
  colnames(genes) <- "gene_id"
  
  DF_NONSV <- FDR_results1 %>% filter(gene_id %notin% genes$gene_id)
  DF_NONSV_pval[[v]] <- DF_NONSV[,(7+1):(7+num_cluster)]
  DF_NONSV_pval[[v]] <- reshape2::melt(DF_NONSV_pval[[v]]) # as rownames, i.e., gene names, are not indentical -> would be discarded
  colnames(DF_NONSV_pval[[v]]) <- c("Layer", "Pvalue")
  #DF_NONSV_pval[[v]]$Layer <- gsub("_FDR.*$","", DF_NONSV_pval[[v]]$Layer)
  DF_NONSV_pval[[v]]$Layer <- gsub("_PValue.*$","", DF_NONSV_pval[[v]]$Layer)
  DF_NONSV_pval[[v]]$Sample <- paste("Sample", sample)
  DF_NONSV_pval[[v]]$Pattern <- pattern
}

library(ggplot2)
DF <- do.call("rbind", DF_NONSV_pval)
(p <- ggplot(data =  DF, aes(x = Pvalue, y = ..ndensity.., lty = Layer)) +
    geom_density(adjust = 0.5, size = 0.3, alpha = 0.1) +
    #facet_wrap(`Pattern` ~ `Sample`))
    facet_grid(Pattern ~ Sample))