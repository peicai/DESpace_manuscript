rm(list =ls())
# source("./05_individual_clusters/01_Plots_function.R")
# source("./05_individual_clusters/01_Individual_clusters_source_code.R")
library(dplyr)
data = "LIBD"
sample_names <- c("151507", "151669", "151673")
spatial_probs = c(0.5,0.9)
default = FALSE
#pattern = "mixture_patch"
pattern_list <- c('mixture_patch','mixture_reverse_patch')

Vectors <- expand.grid(sample_names,
                       pattern_list)
`%notin%` <- Negate(`%in%`)
DF_rank <- list()
methodMarker = "Seurat" 
cluster_method = "BayesSpace"
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
for(v in 1:nrow(Vectors)){
  sample <- Vectors[v,1]
  pattern <- Vectors[v,2]
  path = paste0("./individual_SV_cluster/", data,
                "/",sample, "_", pattern, "_", cluster_method, "_results_",methodMarker, ".rda")
  ## Results
  load(path)
  select = "p_val"
  com_results <- .top_results_seurat(results = results, select = select)
  num_cluster <- dim(com_results)[2]-5
  sub_results <- com_results[,-c(1:5)]
  sub_results$gene_id <- rownames(sub_results)
  ## SV genes
  dir = paste0("./simulation/results/",sample,'/',pattern,'/probs_',spatial_probs[1],'_',
               spatial_probs[2],default,'/')
  genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
  colnames(genes) <- "gene_id"
  
  ## Select SV genes only:
  tab = as.data.frame( merge(sub_results, genes, by="gene_id"))
  pat <- "^pattern_.*?([0-9]).*"
  genes$GoundTruth <-  unlist(lapply(rownames(genes),function(x) gsub(pat, "\\1", x)))
  genes$GoundTruth <- paste0("Layer",as.numeric(genes$GoundTruth)+1)
  genes$GoundTruth[genes$GoundTruth == "Layer7"] <- "WM"
  in_cluster <- genes$GoundTruth
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
  } # end if 151673
  genes$in_cluster <- in_cluster
  # rank based on p values
  ord = as.data.frame(t(apply(tab[,-1], 1, rank)))
  ord$gene_id <- tab$gene_id
  tab2 <- merge(ord, genes, by = "gene_id")
  in_rank <- apply(data.frame(1:dim(tab2)[1]), 1, function(row){
    in_clu <- tab2$in_cluster[row]
    (in_cluster_rank <- tab2[row,in_clu])
  })
  DF_rank[[v]] = data.frame(rank = in_rank,
                            pattern = pattern,
                            sample = sample,
                            cluster_method)
}
DF_rank_all <- do.call(rbind, DF_rank)
(table <- aggregate(DF_rank_all$rank == 1, list(DF_rank_all$sample, DF_rank_all$pattern), FUN=mean))
colnames(table) <- c("Sample","Pattern","Accuracy")
write.table(table, sep = ",", file = paste0('./individual_SV_cluster/LIBD_',cluster_method, '_rank_table_score_',methodMarker,'.csv'))

(table <- aggregate(DF_rank_all$rank %in% c(1, 1.5), list(DF_rank_all$sample, DF_rank_all$pattern), FUN=mean))
colnames(table) <- c("Sample","Pattern","Accuracy")

################################# stLearn ###############
###########

DF_rank <- list()
methodMarker = "Seurat" 
cluster_method = "StLearn"
for(v in 1:nrow(Vectors)){
  sample <- Vectors[v,1]
  pattern <- Vectors[v,2]
  path = paste0("./individual_SV_cluster/", data,
                "/",sample, "_", pattern, "_", cluster_method, "_results_",methodMarker, ".rda")
  ## Results
  load(path)
  select = "p_val"
  com_results <- .top_results_seurat(results = results, select = select)
  num_cluster <- dim(com_results)[2]-5
  sub_results <- com_results[,-c(1:5)]
  sub_results$gene_id <- rownames(sub_results)
  ## SV genes
  dir = paste0("./simulation/results/",sample,'/',pattern,'/probs_',spatial_probs[1],'_',
               spatial_probs[2],default,'/')
  genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
  colnames(genes) <- "gene_id"
  
  ## Select SV genes only:
  tab = as.data.frame( merge(sub_results, genes, by="gene_id"))
  pat <- "^pattern_.*?([0-9]).*"
  genes$GoundTruth <-  unlist(lapply(rownames(genes),function(x) gsub(pat, "\\1", x)))
  genes$GoundTruth <- paste0("Layer",as.numeric(genes$GoundTruth)+1)
  genes$GoundTruth[genes$GoundTruth == "Layer7"] <- "WM"
  in_cluster <- genes$GoundTruth
  if(sample == "151507"){
    if(pattern == "mixture_patch"){
      if(cluster_method == "BayesSpace"){
        in_cluster[in_cluster == "Layer3"] <- 4
        in_cluster[in_cluster == "Layer4"] <- 7
        in_cluster[in_cluster == "Layer5"] <- 2
        in_cluster[in_cluster == "Layer6"] <- 6
        in_cluster[in_cluster == "WM"] <- 3
      } # end if BayesSpce
      if(cluster_method == "StLearn"){
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
      if(cluster_method == "StLearn"){
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
      if(cluster_method == "StLearn"){
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
      if(cluster_method == "StLearn"){
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
      if(cluster_method == "StLearn"){
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
      if(cluster_method == "StLearn"){
        in_cluster[in_cluster == "Layer3"] <- 2
        in_cluster[in_cluster == "Layer4"] <- 0
        in_cluster[in_cluster == "Layer5"] <- 3
        in_cluster[in_cluster == "Layer6"] <- 1
        in_cluster[in_cluster == "WM"] <- 4
      } # end if stLearn
    } # end if mixture reverse pattern
  } # end if 151673
  genes$in_cluster <- in_cluster
  # rank based on p values
  ord = as.data.frame(t(apply(tab[,-1], 1, rank)))
  ord$gene_id <- tab$gene_id
  tab2 <- merge(ord, genes, by = "gene_id")
  in_rank <- apply(data.frame(1:dim(tab2)[1]), 1, function(row){
    in_clu <- tab2$in_cluster[row]
    (in_cluster_rank <- tab2[row,in_clu])
  })
  DF_rank[[v]] = data.frame(rank = in_rank,
                            pattern = pattern,
                            sample = sample,
                            cluster_method)
}
DF_rank_all <- do.call(rbind, DF_rank)
(table <- aggregate(DF_rank_all$rank == 1, list(DF_rank_all$sample, DF_rank_all$pattern), FUN=mean))
colnames(table) <- c("Sample","Pattern","Accuracy")
write.table(table, sep = ",", file = paste0('./individual_SV_cluster/LIBD_',cluster_method, '_rank_table_score_',methodMarker,'.csv'))

(table <- aggregate(DF_rank_all$rank %in% c(1, 1.5), list(DF_rank_all$sample, DF_rank_all$pattern), FUN=mean))
colnames(table) <- c("Sample","Pattern","Accuracy")

