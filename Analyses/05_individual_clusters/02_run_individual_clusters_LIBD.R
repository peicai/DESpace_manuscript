rm(list = ls())
source("./Analyses/03_simulations/02_simulation_function_data_specific.R")
source("./Analyses/05_individual_clusters/01_Individual_clusters_source_code.R")

.individual_test_wrap <- function(sce, layer_names, sample_names, pattern_names, times_df){
  cluster_method <- NULL
  if(layer_names == "spatial.cluster") cluster_method <- "BayesSpace"
  if(layer_names == "stLearn_pca_kmeans") cluster_method <- "stLearn"
  
  layer = colData(sce)[[layer_names]]
  colData(sce)[[layer_names]] = as.factor(colData(sce)[[layer_names]])
  sce = sce[, !is.na(layer)]
  one_layer = list()
  
  y <- DGEList(counts=assays(sce)$counts, 
               genes=rownames(assays(sce)$counts))
  
  y$samples$lib.size <- colSums(y$counts)
  y <- calcNormFactors(y)
  times = c()
  ################################ Option1: with estimated dispersion
  set.seed(123)
  registerDoParallel(cores=4)
  layer = colData(sce)[[layer_names]]
  design = data.frame(condition = layer)
  design$condition <- droplevels(design$condition)
  design_model <- model.matrix(~design$condition)
  rownames(design_model) <- colnames(y)
  # print(design_model)
  # 
  y1 <- estimateDisp(y, design_model, robust=TRUE, 
                     BPPARAM = BPPARAM)
  a = system.time({cluster_results <- individual_test(sce, edgeR_y = y1, layer_col = layer_names)})
  times[1] = a[3]
  names(times[1]) = "With_dispersion"
  results = cluster_results
  save(results, file = paste0("./Data/Simulation/Ouput/individual_SV_cluster/LIBD/",sample_names,"_",pattern_names,"_",cluster_method,"_results_with.rda"))
  
  ################################ Option2: without estimated dispersion
  
  set.seed(123)
  registerDoParallel(cores=4)
  # test 
  b = system.time({ cluster_results2 <- individual_test(sce, edgeR_y = NULL, layer_col = layer_names)})
  times[2] = b[3]
  names(times[2]) = "Without_disperison" 
  
  
  results2 = cluster_results2
  save(results2, file = paste0("./Data/Simulation/Output/individual_SV_cluster/LIBD/",sample_names,"_",pattern_names,"_",cluster_method,"_results_without.rda")) 
  
  
  times[3] = sample_names
  names(times[3]) = "sample"
  # 
  times[4] = pattern_names
  names(times[4]) = "pattern"
  # 
  times[5] = cluster_method
  names(times[5]) = "cluster_source"
  # 
  times_df <- rbind(times_df,times)
  return(times_df)
}

sample_names_all  = c("151507","151669","151673")
pattern_names_all = c("mixture_patch", "mixture_reverse_patch")
cluster_levels <-list("Layer3", "Layer4", "Layer5", "Layer6", "WM")

Vectors <- expand.grid(vector1 = c(1:3),
                       vector2 = c(1:2))
times_df1 = NULL
times_df2 = NULL
for(v in 1:nrow(Vectors)){
  i <- Vectors[v,1]
  j <- Vectors[v,2]
  sample_names = sample_names_all[i]
  pattern_names = pattern_names_all[j]
  
  ## load data
  path = paste0("./Data/Simulation/Output/",sample_names,"/",pattern_names,"/probs_0.5_0.9FALSE/")
  genes = read.table(paste0(path,"probs_0.9_selected_genes.txt"))
  #load(paste0(path, "probs_0.5_0.9_final_object.rda"))
  # sce object contains clusters information from BayesSpace (obtained from the main simulation)
  load(paste0(path, "sce_edgeR.rda"))
  layer_names = "spatial.cluster"
  print("BayesSpace")
  times_df1 <- .individual_test_wrap(sce, layer_names, sample_names, pattern_names, times_df1)
  rm(sce)
  load(paste0(path, "stLearn_sce_edgeR.rda"))
  layer_names = "stLearn_pca_kmeans"
  print("StLearn")
  times_df2 <- .individual_test_wrap(sce, layer_names, sample_names, pattern_names, times_df2)
}
times_df <- rbind(times_df1, times_df2)
write.csv(times_df, quote = FALSE, row.names = TRUE,
          file = paste0("LIBD_single_cluster_times.csv"))