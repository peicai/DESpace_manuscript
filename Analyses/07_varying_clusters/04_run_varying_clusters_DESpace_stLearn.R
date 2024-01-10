rm(list = ls())
################## stLearn (after getting stLearn cluster results) #########

save_dir <- "./Simulation/Output/151507/"
spatial_probs <- c(0.5,0.9); default = FALSE
for(i in c(2,4,6,8,10,12)){
  pattern_name <- paste0("split_mixture_patch_", i)
  load(paste0(save_dir, pattern_name, "/probs_0.5_0.9FALSE/sce_edgeR.rda"))
  set.seed(123)
  prob = method = adj.p.value = time = NULL
  simobject = sce
  subdir = paste0(save_dir,'/',pattern_name,'/')
  resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
  
  all_results = list()
  
  rep_list = list()
  st_result_path <- paste0(resultdir, "stLearn_results.csv")
  
  stLearn_results <- read.csv(st_result_path, sep = ',', header = TRUE)
  colnames(stLearn_results) <- c("cell_ID","imagecol","imagerow","stLearn_pca_kmeans","stLearn_louvain")
  
  final_object <- addMetadata(
    final_object,
    stLearn_results,
    by_column = TRUE,
    column_cell_ID = "cell_ID"
  )
  
  
  spatial_gene_results <- SV_edgeR_counts(final_object,
                                          layer_names = "stLearn_louvain",
                                          original_layer_names = "pattern",
                                          dir = resultdir,           
                                          pattern_name = pattern_name,
                                          cluster_method='stLearn',
                                          platform = "Visium",
                                          return_object = 'data.table')
  spatial_gene_results[, method := 'stLearn']
  save(spatial_gene_results, file = paste0(resultdir,"result_SV_edgeR_stLearn_louvain.rda"))
  print(i)
}
