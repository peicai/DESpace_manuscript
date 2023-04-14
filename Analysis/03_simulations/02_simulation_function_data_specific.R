source("./Analysis/03_simulations/01_simulation_basic_function.R")

Melanoma_simulation <- function(i,List_sce_one,
                                coordinate_name = c("row","col"),
                                working_directory = paste0("~/simulation/results/",sample_names[i]),
                                pattern_type = 'BayesSpace_clusters_patch',
                                my_dir = paste0("~/Desktop/simulation/results/",sample_names[i]),
                                spatial_probs = c(0.6,0.9),
                                SV_percent = 0.5,
                                spat_methods = c('SV_edgeR_NoNorm',
                                                 'SV_edgeR_CPMs','SV_edgeR_CPMs_NoNorm'
                                ),
                                spat_methods_params = list(NA,NA, NA),
                                spat_methods_names = c('SV_edgeR_NoNorm',
                                                       'SV_edgeR_CPMs','SV_edgeR_CPMs_NoNorm'
                                ),
                                layer_names = "spatial.cluster",
                                original_layer_names = "spatial.cluster",
                                default = T,
                                reps =1,
                                cluster_index=NULL, # for multiple patterns, list
                                platform = "ST",
                                ...
){
  set.seed(123)
  sce_one = List_sce_one[[i]]
  exprs <- assays(sce_one)$counts
  sce_one <- addPerCellQCMetrics(sce_one)
  sce_one <- logNormCounts(sce_one)
  CD <- colData(sce_one)
  cell_ID = rownames(colData(sce_one))
  colData(sce_one) <- cbind(CD, cell_ID)
  sample_genes = rownames(exprs)
  
  object_subset = sce_one[sample_genes,]
  pattern_type <- pattern_type # 'right_patch', 'circle_patch', 'BayesSpace_clusters_patch'
  if(i==1){combined_cluster <- c(1,4)}
  if(i==2){combined_cluster <- c(2,4)}
  if(i==3){combined_cluster <- c(1,4)}
  if(i==4){combined_cluster <- c(1,3)}
  if(i==5){combined_cluster <- c(1,4)}
  if(i==6){combined_cluster <- c(1,3)}
  if(i==7){combined_cluster <- c(1,2)}
  if(i==8){combined_cluster <- c(1,2)}
  print(pattern_type)
  covariates = cbind.data.frame(x=as.numeric(as.character(colData(object_subset)[[coordinate_name[1]]])),
                                y=as.numeric(as.character(colData(object_subset)[[coordinate_name[2]]])),
                                cell_ID = rownames(colData(object_subset))
  )
  colnames(covariates) <- c("sdimx","sdimy","cell_ID")
  if(pattern_type == 'right_patch'){
    # pattern 1: bottom right stripe
    
    pattern = covariates %>% filter(sdimx > 0 & sdimy < mean(covariates$sdimy))
    pattern_ids = pattern$cell_ID
    pattern_name = 'right_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    #cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    #cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    
    cluster_index = NULL
  } else if (pattern_type == 'circle_patch'){
    # pattern 2: circle pattern
    mx <- round(mean(covariates$sdimx))
    my <- round(mean(covariates$sdimy))
    dist_x <- (covariates$sdimx - mx)^2 
    dist_y <-  (covariates$sdimy - my)^2
    pattern = covariates[dist_x + dist_y < 25,]
    pattern_ids = pattern$cell_ID
    pattern_name = 'circle_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    #cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    #cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    cluster_index = NULL
  } else if (pattern_type == 'BayesSpace_clusters_patch'){
    # pattern 3: BayesSpace clusters
    colData(sce_one)$spatial.cluster
    #spatPlot2D(gobject = giotto_object_subset, save_plot = T, #cell_color_code = colData(sce_one)$spatial.cluster,
    #             point_size = 0.9, cell_color = as.factor(colData(sce_one)$spatial.cluster), default_save_name  = "BayesSpace_clusters")
    pattern = covariates[as.numeric(colData(sce_one)$spatial.cluster) %in% combined_cluster,]
    pattern_ids = pattern$cell_ID
    pattern_name = 'BayesSpace_clusters_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    #cell_meta = as.data.table(colData(object_subset))
    # cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    #cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    cluster_index = NULL
  } else if(pattern_type == "MixClusters_patch"){
    pattern_ids <- list()
    pattern_name = 'mixture_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #  cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    pattern_ids <- lapply(cluster_index,
                          function(x) covariates[as.numeric(colData(sce_one)$spatial.cluster) %in% x,'cell_ID'] )
    print(cluster_index)
    print(length(pattern_ids))
    #for(ii in 1:length(pattern_ids)){
    #  x = pattern_ids[[ii]]
    #  print(x)
    #  cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
    
    #}
  } else if(pattern_type == "MixClusters_reverse_patch"){
    pattern_ids <- list()
    pattern_name = 'mixture_reverse_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #  cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    pattern_ids <- lapply(cluster_index,
                          function(x) covariates[as.numeric(colData(sce_one)$spatial.cluster) %notin% x,'cell_ID'] )
    print(cluster_index)
    print(length(pattern_ids))
    # for(ii in 1:length(pattern_ids)){
    #  x = pattern_ids[[ii]]
    #print(x)
    #  cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
    #}
  } 
  print(dim(counts(object_subset)))
  subdir = paste0(my_dir,'/',pattern_name,'/')
  resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
  if(!file.exists(resultdir)) dir.create(path = resultdir, recursive = TRUE)
  save_name = 'probs'
  if(file.exists(paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))) {
    message("Simulation dataset has been created. 
            Skip data simulation part.
            Load the simulation dataset...")
    load(paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))
    if(exists('final_object')) 
      final_object <-  final_object
    else 
      final_object <- combined_object
  }
  
  selected_genes = sample_genes
  my_dir = my_dir
  
  save_dir = my_dir
  runPatternSimulation(sce_object = object_subset,
                       rawobject = object_subset,
                       coordinate_name = coordinate_name,
                       pattern_name = pattern_name,
                       pattern_cell_ids = pattern_ids,
                       gene_names = selected_genes,
                       spatial_probs = spatial_probs,
                       SV_percent = SV_percent,
                       spat_methods = spat_methods,
                       spat_methods_params = spat_methods_params,
                       spat_methods_names = spat_methods_names,
                       layer_names = layer_names,
                       original_layer_names = original_layer_names,
                       reps = reps,
                       save_plot = T,
                       save_data = TRUE,
                       max_col = 2,
                       save_dir = my_dir,
                       default = default,
                       cluster_index = cluster_index,
                       platform = platform)
  
  
  
}

LIBD_simulation <- function(i,List_sce_one,
                            coordinate_name = c("imagerow","imagecol"),
                            working_directory = paste0("~/simulation/results/",sample_names[i]),
                            pattern_type = 'BayesSpace_clusters_patch',
                            my_dir = paste0("~/Desktop/simulation/results/",sample_names[i]),
                            spatial_probs = c(0.6,0.9),
                            SV_percent = 0.5,
                            spat_methods = c('SV_edgeR_NoNorm',
                                             'SV_edgeR_CPMs','SV_edgeR_CPMs_NoNorm'
                            ),
                            spat_methods_params = list(NA,NA, NA),
                            spat_methods_names = c('SV_edgeR_NoNorm',
                                                   'SV_edgeR_CPMs','SV_edgeR_CPMs_NoNorm'
                            ),
                            layer_names = "spatial.cluster",
                            original_layer_names = "spatial.cluster",
                            default = T,
                            reps =1,
                            cluster_index=NULL, # for multiple patterns, list
                            platform = "Visium",
                            ...
){
  set.seed(123)
  sce_one = List_sce_one[[i]]
  exprs <- assays(sce_one)$counts
  sce_one <- addPerCellQCMetrics(sce_one)
  sce_one <- logNormCounts(sce_one)
  CD <- colData(sce_one)
  cell_ID = rownames(colData(sce_one))
  colData(sce_one) <- cbind(CD, cell_ID)
  sample_genes = rownames(exprs)
  
  object_subset = sce_one[sample_genes,]
  pattern_type <- pattern_type # 'right_patch', 'circle_patch', 'BayesSpace_clusters_patch'
  # Manual
  if(i == 1){combined_cluster <- c(1,4,5)}
  if(i == 5){combined_cluster <- c(4,7)}
  if(i == 9){combined_cluster <- c(2,3,4)}
  
  # BayesSpace
  if(i == 1){combined_cluster2 <- c(1,5,7)}
  if(i == 5){combined_cluster2 <- c(2,3,4,5)}
  if(i == 9){combined_cluster2 <- c(3,4,6,7)}
  print(pattern_type)
  covariates = cbind.data.frame(x=as.numeric(as.character(colData(object_subset)[[coordinate_name[1]]])),
                                y=as.numeric(as.character(colData(object_subset)[[coordinate_name[2]]])),
                                cell_ID = rownames(colData(object_subset))
  )
  colnames(covariates) <- c("sdimx","sdimy","cell_ID")
  if(pattern_type == 'bottom_patch'){
    # pattern 1: bottom right stripe
    
    pattern = covariates %>% filter(sdimx > 0 & sdimy < mean(covariates$sdimy))
    pattern_ids = pattern$cell_ID
    pattern_name = 'bottom_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    #cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    #cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    
    cluster_index = NULL
  } else if (pattern_type == 'circle_patch'){
    # pattern 2: circle pattern
    mx <- round(mean(covariates$sdimx))
    my <- round(mean(covariates$sdimy))
    dist_x <- (covariates$sdimx - mx)^2 
    dist_y <-  (covariates$sdimy - my)^2
    pattern = covariates[dist_x + dist_y < 10000,]
    pattern_ids = pattern$cell_ID
    pattern_name = 'circle_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    #cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    #cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    cluster_index = NULL
  } else if (pattern_type == 'Manual_clusters_patch'){
    # pattern 3: BayesSpace clusters
    colData(sce_one)$spatial.cluster
    #spatPlot2D(gobject = giotto_object_subset, save_plot = T, #cell_color_code = colData(sce_one)$spatial.cluster,
    #             point_size = 0.9, cell_color = as.factor(colData(sce_one)$spatial.cluster), default_save_name  = "BayesSpace_clusters")
    pattern = covariates[as.numeric(colData(sce_one)$layer_guess_reordered) %in% combined_cluster,]
    pattern_ids = pattern$cell_ID
    pattern_name = 'Manual_clusters_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    #cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    #cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    cluster_index = NULL
  } else if(pattern_type == "MixClusters_patch"){
    pattern_ids <- list()
    pattern_name = 'mixture_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    pattern_ids <- lapply(cluster_index,
                          function(x) covariates[as.numeric(colData(sce_one)$layer_guess_reordered) %in% x,'cell_ID'] )
    print(cluster_index)
    print(length(pattern_ids))
    # for(ii in 1:length(pattern_ids)){
    #  x = pattern_ids[[ii]]
    #   print(x)
    #   cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
    
  } else if(pattern_type == "MixClusters_reverse_patch"){
    pattern_ids <- list()
    pattern_name = 'mixture_reverse_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    pattern_ids <- lapply(cluster_index,
                          function(x) covariates[as.numeric(colData(sce_one)$layer_guess_reordered) %notin% x,'cell_ID'] )
    print(cluster_index)
    print(length(pattern_ids))
    #for(ii in 1:length(pattern_ids)){
    #  x = pattern_ids[[ii]]
    #print(x)
    #  cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
    #}
  } 
  
  print(dim(counts(object_subset)))
  subdir = paste0(my_dir,'/',pattern_name,'/')
  resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
  if(!file.exists(resultdir)) dir.create(path = resultdir, recursive = TRUE)
  save_name = 'probs'
  if(file.exists(paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))) {
    message("Simulation dataset has been created. \n
            Skip data simulation part.\n
            Load the simulation dataset...\n")
    load(paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))
    if(exists('final_object')) 
      final_object <-  final_object
    else 
      final_object <- combined_object
    
    selected_genes <- rownames(counts(final_object))
    selected_cells <- colnames(counts(final_object))
    object_subset = object_subset[,selected_cells]
    print(dim(counts(object_subset)))
    
  }else {
    set.seed(123)
    selected_genes = sample(sample_genes, 5000)
  }
  
  mini_object = object_subset[selected_genes,]
  
  
  
  my_dir = my_dir
  
  save_dir = my_dir
  runPatternSimulation(sce_object = mini_object,
                       rawobject = mini_object,
                       coordinate_name = coordinate_name,
                       pattern_name = pattern_name,
                       pattern_cell_ids = pattern_ids,
                       gene_names = selected_genes,
                       spatial_probs = spatial_probs,
                       SV_percent = SV_percent,
                       spat_methods = spat_methods,
                       spat_methods_params = spat_methods_params,
                       spat_methods_names = spat_methods_names,
                       layer_names = layer_names,
                       original_layer_names = original_layer_names,
                       reps = reps,
                       save_plot = T,
                       save_data = TRUE,
                       max_col = 2,
                       save_dir = my_dir,
                       default = default,
                       cluster_index = cluster_index,
                       platform = platform)
  
  
  
}

Slideseq_simulation <- function(sce_one,
                                coordinate_name = c("row","col"),
                                working_directory = paste0("~/simulation/results/",sample_names[i]),
                                pattern_type = 'BayesSpace_clusters_patch',
                                my_dir = paste0("~/Desktop/simulation/results/",sample_names[i]),
                                spatial_probs = c(0.6,0.9),
                                SV_percent = 0.5,
                                spat_methods = c('SV_edgeR_NoNorm',
                                                 'SV_edgeR_CPMs','SV_edgeR_CPMs_NoNorm'
                                ),
                                spat_methods_params = list(NA,NA, NA),
                                spat_methods_names = c('SV_edgeR_NoNorm',
                                                       'SV_edgeR_CPMs','SV_edgeR_CPMs_NoNorm'
                                ),
                                layer_names = "spatial.cluster",
                                original_layer_names = "spatial.cluster",
                                default = T,
                                reps =1,
                                cluster_index=NULL, # for multiple patterns, list
                                platform = "Visium",
                                ...
){
  set.seed(123)
  exprs <- assays(sce_one)$counts
  sce_one <- addPerCellQCMetrics(sce_one)
  sce_one <- logNormCounts(sce_one)
  CD <- colData(sce_one)
  cell_ID = rownames(colData(sce_one))
  colData(sce_one) <- cbind(CD, cell_ID)
  st_result_path <- paste0('./Simulation/Input/mouse_cerebellum/100_stLearn_results.csv')
  #st_result_path <- paste0('/cluster/scratch/peicai/SlideSeq2/stLearn_q4.csv')
  stLearn_results <- read.csv(st_result_path, sep = ',', header = TRUE)
  colnames(stLearn_results) <- c("cell_ID","imagecol","imagerow","stLearn_pca_kmeans")
  
  CD <- colData(sce_one)
  colData(sce_one) <- cbind(CD, stLearn_results$stLearn_pca_kmeans)
  
  sample_genes = rownames(exprs)
  
  object_subset = sce_one[sample_genes,]
  pattern_type <- pattern_type # 'right_patch', 'circle_patch', 'BayesSpace_clusters_patch'
  combined_cluster <- c(3)
  print(pattern_type)
  covariates = cbind.data.frame(x=as.numeric(as.character(colData(object_subset)[[coordinate_name[1]]])),
                                y=as.numeric(as.character(colData(object_subset)[[coordinate_name[2]]])),
                                cell_ID = rownames(colData(object_subset))
  )
  colnames(covariates) <- c("sdimx","sdimy","cell_ID")
  if(pattern_type == 'bottom_patch'){
    # pattern 1: bottom right stripe
    #thres = summary(covariates$sdimy)[6]
    pattern = covariates %>% filter(sdimy > 0 & sdimx > 4470)#thres) # 15.6%
    pattern_ids = pattern$cell_ID
    pattern_name = 'bottom_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    
    cluster_index = NULL
    ggplot(cell_meta, aes(x = col, y = row)) +
      geom_point(aes(color=bottom_patch),size=0.1,shape=".")+coord_flip()+ scale_x_reverse()
    
  } else if (pattern_type == 'circle_patch'){
    # pattern 2: circle pattern
    mx <- round(mean(covariates$sdimx)+1300)
    my <- round(mean(covariates$sdimy)-700)
    dist_x <- (covariates$sdimx - mx)^2 
    dist_y <-  (covariates$sdimy - my)^2
    pattern = covariates[dist_x + dist_y < 450000,]
    pattern_ids = pattern$cell_ID
    pattern_name = 'circle_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    cluster_index = NULL
    ggplot(cell_meta, aes(x = col, y = row)) +
      geom_point(aes(color=circle_patch),size=0.1,shape=".")+coord_flip()+ scale_x_reverse()
    
  } else if (pattern_type == 'StLearn_clusters_patch'){
    # pattern 3: BayesSpace clusters
    pattern = covariates[as.numeric(colData(sce_one)$`stLearn_results$stLearn_pca_kmeans`) %in% combined_cluster,]
    pattern_ids = pattern$cell_ID
    pattern_name = 'StLearn_clusters_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'out', 'in')]
    cluster_index = NULL
    ggplot(cell_meta, aes(x = col, y = row)) +
      geom_point(aes(color=StLearn_clusters_patch),size=0.1,shape=".")+coord_flip()+ scale_x_reverse()
    cell_number = nrow(cell_meta[get(pattern_name) == 'in'])
    print(cell_number)
    
  } else if(pattern_type == "MixClusters_patch"){
    pattern_ids <- list()
    pattern_name = 'mixture_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    pattern_ids <- lapply(cluster_index,
                          function(x) covariates[as.numeric(colData(sce_one)[['stLearn_results$stLearn_pca_kmeans']]) %in% x,'cell_ID'] )
    
    for(ii in 1:length(pattern_ids)){
      x = pattern_ids[[ii]]
      print(x)
      cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
      
    }} else if(pattern_type == "MixClusters_reverse_patch"){
      pattern_ids <- list()
      pattern_name = 'mixture_reverse_patch'
      pattern_colors = c('in' = 'green', 'out' = 'red')
      cell_meta = as.data.table(colData(object_subset))
      #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
      pattern_ids <- lapply(cluster_index,
                            function(x) covariates[as.numeric(colData(sce_one)[['stLearn_results$stLearn_pca_kmeans']]) %notin% x,'cell_ID'] )
      print(cluster_index)
      print(length(pattern_ids))
      for(ii in 1:length(pattern_ids)){
        x = pattern_ids[[ii]]
        #print(x)
        cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
      }} 
  cell_meta <- DataFrame(cell_meta)
  rownames(cell_meta) <- rownames(colData(object_subset))
  colData(object_subset) <-  cell_meta
  print(dim(counts(object_subset)))
  subdir = paste0(my_dir,'/',pattern_name,'/')
  resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
  if(!file.exists(resultdir)) dir.create(path = resultdir, recursive = TRUE)
  save_name = 'probs'
  if(file.exists(paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))) {
    message("Simulation dataset has been created. 
            Skip data simulation part.
            Load the simulation dataset...")
    load(paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))
    if(exists('final_object')) 
      final_object <-  final_object
    else 
      final_object <- combined_object
    
    selected_genes <- rownames(counts(final_object))
    selected_cells <- colnames(counts(final_object))
    object_subset = object_subset[,selected_cells]
    print(dim(counts(object_subset)))
    
  }else {
    set.seed(123)
    selected_genes = sample(sample_genes, 5000)
  }
  
  mini_object = object_subset[selected_genes,]
  
  
  
  my_dir = my_dir
  
  save_dir = my_dir
  runPatternSimulation(sce_object = mini_object,
                       rawobject = mini_object,
                       coordinate_name = coordinate_name,
                       pattern_name = pattern_name,
                       pattern_cell_ids = pattern_ids,
                       gene_names = selected_genes,
                       spatial_probs = spatial_probs,
                       SV_percent = SV_percent,
                       spat_methods = spat_methods,
                       spat_methods_params = spat_methods_params,
                       spat_methods_names = spat_methods_names,
                       layer_names = layer_names,
                       original_layer_names = original_layer_names,
                       reps = reps,
                       save_plot = T,
                       save_data = TRUE,
                       max_col = 2,
                       save_dir = my_dir,
                       default = default,
                       cluster_index = cluster_index,
                       platform = platform)
  
  
  
}

Slideseq_sub_simulation <- function(sce_one,
                                    coordinate_name = c("row","col"),
                                    working_directory = paste0("~/simulation/results/",sample_names[i]),
                                    pattern_type = 'BayesSpace_clusters_patch',
                                    my_dir = paste0("~/Desktop/simulation/results/",sample_names[i]),
                                    spatial_probs = c(0.6,0.9),
                                    SV_percent = 0.5,
                                    spat_methods = c('SV_edgeR_NoNorm',
                                                     'SV_edgeR_CPMs','SV_edgeR_CPMs_NoNorm'
                                    ),
                                    spat_methods_params = list(NA,NA, NA),
                                    spat_methods_names = c('SV_edgeR_NoNorm',
                                                           'SV_edgeR_CPMs','SV_edgeR_CPMs_NoNorm'
                                    ),
                                    layer_names = "spatial.cluster",
                                    original_layer_names = "spatial.cluster",
                                    default = T,
                                    reps =1,
                                    cluster_index=NULL, # for multiple patterns, list
                                    platform = "Visium",
                                    ...
){
  set.seed(123)
  exprs <- assays(sce_one)$counts
  sce_one <- addPerCellQCMetrics(sce_one)
  sce_one <- logNormCounts(sce_one)
  CD <- colData(sce_one)
  cell_ID = rownames(colData(sce_one))
  colData(sce_one) <- cbind(CD, cell_ID)
  #st_result_path <- paste0('~/Desktop/Real_data/data/Slideseq2/stLearn_q4.csv')
  st_result_path <- paste0('./Simulation/Input/mouse_cerebellum/100_stLearn_results.csv')
  #st_result_path <- paste0('/cluster/scratch/peicai/SlideSeq2/stLearn_q4.csv')
  stLearn_results <- read.csv(st_result_path, sep = ',', header = TRUE)
  colnames(stLearn_results) <- c("cell_ID","imagecol","imagerow","stLearn_pca_kmeans")
  
  CD <- colData(sce_one)
  colData(sce_one) <- cbind(CD, stLearn_results$stLearn_pca_kmeans)
  
  sample_genes = rownames(exprs)
  
  object_subset = sce_one[sample_genes,]
  pattern_type <- pattern_type # 'right_patch', 'circle_patch', 'BayesSpace_clusters_patch'
  combined_cluster <- c(0,2)
  print(pattern_type)
  covariates = cbind.data.frame(x=as.numeric(as.character(colData(object_subset)[[coordinate_name[1]]])),
                                y=as.numeric(as.character(colData(object_subset)[[coordinate_name[2]]])),
                                cell_ID = rownames(colData(object_subset))
  )
  colnames(covariates) <- c("sdimx","sdimy","cell_ID")
  if(pattern_type == 'bottom_patch'){
    # pattern 1: bottom right stripe
    
    pattern = covariates %>% filter(sdimx > 0 & sdimy < mean(covariates$sdimy))
    pattern_ids = pattern$cell_ID
    pattern_name = 'bottom_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    
    cluster_index = NULL
  } else if (pattern_type == 'circle_patch'){
    # pattern 2: circle pattern
    mx <- round(mean(covariates$sdimx))
    my <- round(mean(covariates$sdimy)-500)
    dist_x <- (covariates$sdimx - mx)^2 
    dist_y <-  (covariates$sdimy - my)^2
    pattern = covariates[dist_x + dist_y < 850000,]
    pattern_ids = pattern$cell_ID
    pattern_name = 'circle_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    cluster_index = NULL
    #ggplot(cell_meta, aes(x = col, y = row)) +
    #  geom_point(aes(color=circle_patch),size=0.1,shape=".")+coord_flip()+ scale_x_reverse()
    
  } else if (pattern_type == 'StLearn_clusters_patch'){
    # pattern 3: BayesSpace clusters
    pattern = covariates[as.numeric(colData(sce_one)$stLearn_pca_kmeans) %in% combined_cluster,]
    pattern_ids = pattern$cell_ID
    pattern_name = 'StLearn_clusters_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    cluster_index = NULL
    #ggplot(cell_meta, aes(x = col, y = row)) +
    #    geom_point(aes(color=StLearn_clusters_patch),size=0.1,shape=".")+coord_flip()+ scale_x_reverse()
    
  } else if(pattern_type == "MixClusters_patch"){
    pattern_ids <- list()
    pattern_name = 'mixture_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    pattern_ids <- lapply(cluster_index,
                          function(x) covariates[as.numeric(colData(sce_one)[['stLearn_results$stLearn_pca_kmeans']]) %in% x,'cell_ID'] )
    
    for(ii in 1:length(pattern_ids)){
      x = pattern_ids[[ii]]
      print(x)
      cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
      
    }} else if(pattern_type == "MixClusters_reverse_patch"){
      pattern_ids <- list()
      pattern_name = 'mixture_reverse_patch'
      pattern_colors = c('in' = 'green', 'out' = 'red')
      cell_meta = as.data.table(colData(object_subset))
      #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
      pattern_ids <- lapply(cluster_index,
                            function(x) covariates[as.numeric(colData(sce_one)[['stLearn_results$stLearn_pca_kmeans']]) %notin% x,'cell_ID'] )
      print(cluster_index)
      print(length(pattern_ids))
      for(ii in 1:length(pattern_ids)){
        x = pattern_ids[[ii]]
        #print(x)
        cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
      }} 
  
  print(dim(counts(object_subset)))
  subdir = paste0(my_dir,'/',pattern_name,'/')
  resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
  if(!file.exists(resultdir)) dir.create(path = resultdir, recursive = TRUE)
  save_name = 'probs'
  if(file.exists(paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))) {
    message("Simulation dataset has been created. 
            Skip data simulation part.
            Load the simulation dataset...")
    load(paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))
    if(exists('final_object')) 
      final_object <-  final_object
    else 
      final_object <- combined_object
    
    selected_genes <- rownames(counts(final_object))
    selected_cells <- colnames(counts(final_object))
    object_subset = object_subset[,selected_cells]
    print(dim(counts(object_subset)))
    
  }else {
    set.seed(123)
    selected_genes = sample(sample_genes, 5000)
  }
  
  mini_object = object_subset[selected_genes,]
  
  
  
  my_dir = my_dir
  
  save_dir = my_dir
  runStLearnPatternSimulation(sce_object = mini_object,
                              rawobject = mini_object,
                              coordinate_name = coordinate_name,
                              #st_result_path = paste0('~/Desktop/stLearn/Slideseq/simulation/',pattern_name,'/probs_',
                              #                        spatial_probs[1],'_',spatial_probs[2],default,"_stLearn_results.csv"),
                              st_result_path = paste0('./Simulation/Input/mouse_cerebellum/',pattern_name,'/probs_',
                                                      spatial_probs[1],'_',spatial_probs[2],default,"/stLearn_results.csv"),
                              pattern_name = pattern_name,
                              pattern_cell_ids = pattern_ids,
                              gene_names = selected_genes,
                              spatial_probs = spatial_probs,
                              SV_percent = SV_percent,
                              spat_methods = spat_methods,
                              spat_methods_params = spat_methods_params,
                              spat_methods_names = spat_methods_names,
                              layer_names = layer_names,
                              original_layer_names = original_layer_names,
                              reps = reps,
                              save_plot = T,
                              save_data = TRUE,
                              max_col = 2,
                              save_dir = my_dir,
                              default = default,
                              cluster_method = 'stLearn',
                              cluster_index = cluster_index,
                              platform = platform)
  
  
  
}


LIBD_sub_simulation <- function(i,
                                List_sce_one,
                                coordinate_name = c("row","col"),
                                st_result_path = paste0("./Simulation/Input/LIBD/",
                                                        pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],default,"/1k/",sample_names[i],
                                                        "_stLearn_results.csv"),
                                
                                working_directory = paste0("~/simulation/results/",sample_names[i]),
                                pattern_type = 'BayesSpace_clusters_patch',
                                my_dir = paste0("~/Desktop/simulation/results/",sample_names[i]),
                                spatial_probs = c(0.6,0.9),
                                SV_percent = 0.5,
                                spat_methods = c('SV_edgeR_NoNorm',
                                                 'SV_edgeR_CPMs','SV_edgeR_CPMs_NoNorm'
                                ),
                                spat_methods_params = list(NA,NA, NA),
                                spat_methods_names = c('SV_edgeR_NoNorm',
                                                       'SV_edgeR_CPMs','SV_edgeR_CPMs_NoNorm'
                                ),
                                layer_names = "spatial.cluster",
                                original_layer_names = "spatial.cluster",
                                default = T,
                                reps =1,
                                cluster_index=NULL, # for multiple patterns, list
                                platform = "ST",
                                ...
){
  set.seed(123)
  sce_one = List_sce_one[[i]]
  exprs <- assays(sce_one)$counts
  sce_one <- addPerCellQCMetrics(sce_one)
  sce_one <- logNormCounts(sce_one)
  CD <- colData(sce_one)
  cell_ID = rownames(colData(sce_one))
  colData(sce_one) <- cbind(CD, cell_ID)
  sample_genes = rownames(exprs)
  
  object_subset = sce_one[sample_genes,]
  pattern_type <- pattern_type # 'right_patch', 'circle_patch', 'BayesSpace_clusters_patch'
  # Manual
  if(i == 1){combined_cluster <- c(1,4,5)}
  if(i == 5){combined_cluster <- c(4,7)}
  if(i == 9){combined_cluster <- c(2,3,4)}
  
  # BayesSpace
  if(i == 1){combined_cluster2 <- c(1,5,7)}
  if(i == 5){combined_cluster2 <- c(2,3,4,5)}
  if(i == 9){combined_cluster2 <- c(3,4,6,7)}
  print(pattern_type)
  covariates = cbind.data.frame(x=as.numeric(as.character(colData(object_subset)[[coordinate_name[1]]])),
                                y=as.numeric(as.character(colData(object_subset)[[coordinate_name[2]]])),
                                cell_ID = rownames(colData(object_subset))
  )
  colnames(covariates) <- c("sdimx","sdimy","cell_ID")
  if(pattern_type == 'bottom_patch'){
    # pattern 1: bottom right stripe
    
    pattern = covariates %>% filter(sdimx > 0 & sdimy < mean(covariates$sdimy))
    pattern_ids = pattern$cell_ID
    pattern_name = 'bottom_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    
    cluster_index = NULL
  } else if (pattern_type == 'circle_patch'){
    # pattern 2: circle pattern
    mx <- round(mean(covariates$sdimx))
    my <- round(mean(covariates$sdimy))
    dist_x <- (covariates$sdimx - mx)^2 
    dist_y <-  (covariates$sdimy - my)^2
    pattern = covariates[dist_x + dist_y < 10000,]
    pattern_ids = pattern$cell_ID
    pattern_name = 'circle_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    cluster_index = NULL
  } else if (pattern_type == 'Manual_clusters_patch'){
    # pattern 3: BayesSpace clusters
    colData(sce_one)$spatial.cluster
    #spatPlot2D(gobject = giotto_object_subset, save_plot = T, #cell_color_code = colData(sce_one)$spatial.cluster,
    #             point_size = 0.9, cell_color = as.factor(colData(sce_one)$spatial.cluster), default_save_name  = "BayesSpace_clusters")
    pattern = covariates[as.numeric(colData(sce_one)$layer_guess_reordered) %in% combined_cluster,]
    pattern_ids = pattern$cell_ID
    pattern_name = 'Manual_clusters_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    cluster_index = NULL
  } else if(pattern_type == "MixClusters_patch"){
    pattern_ids <- list()
    pattern_name = 'mixture_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    pattern_ids <- lapply(cluster_index,
                          function(x) covariates[as.numeric(colData(sce_one)$layer_guess_reordered) %in% x,'cell_ID'] )
    print(cluster_index)
    print(length(pattern_ids))
    for(ii in 1:length(pattern_ids)){
      x = pattern_ids[[ii]]
      print(x)
      cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
      
    }} else if(pattern_type == "MixClusters_reverse_patch"){
      pattern_ids <- list()
      pattern_name = 'mixture_reverse_patch'
      pattern_colors = c('in' = 'green', 'out' = 'red')
      cell_meta = as.data.table(colData(object_subset))
      #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
      pattern_ids <- lapply(cluster_index,
                            function(x) covariates[as.numeric(colData(sce_one)$layer_guess_reordered) %notin% x,'cell_ID'] )
      print(cluster_index)
      print(length(pattern_ids))
      for(ii in 1:length(pattern_ids)){
        x = pattern_ids[[ii]]
        #print(x)
        cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
      }} 
  
  print(dim(counts(object_subset)))
  subdir = paste0(my_dir,'/',pattern_name,'/')
  resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
  if(!file.exists(resultdir)) dir.create(path = resultdir, recursive = TRUE)
  save_name = 'probs'
  if(file.exists(paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))) {
    message("Simulation dataset has been created. 
            Skip data simulation part.
            Load the simulation dataset...")
    load(paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))
    if(exists('final_object')) 
      final_object <-  final_object
    else 
      final_object <- combined_object
    
    selected_genes <- rownames(counts(final_object))
    selected_cells <- colnames(counts(final_object))
    object_subset = object_subset[,selected_cells]
    print(dim(counts(object_subset)))
    
  }else {
    set.seed(123)
    selected_genes = sample(sample_genes, 5000)
  }
  
  mini_object = object_subset[selected_genes,]
  
  save_dir = my_dir
  runStLearnPatternSimulation(sce_object = mini_object,
                              rawobject = mini_object,
                              coordinate_name = coordinate_name,
                              st_result_path = st_result_path,
                              pattern_name = pattern_name,
                              pattern_cell_ids = pattern_ids,
                              gene_names = selected_genes,
                              spatial_probs = spatial_probs,
                              SV_percent = SV_percent,
                              spat_methods = spat_methods,
                              spat_methods_params = spat_methods_params,
                              spat_methods_names = spat_methods_names,
                              layer_names = layer_names,
                              original_layer_names = original_layer_names,
                              reps = reps,
                              save_plot = T,
                              save_data = TRUE,
                              max_col = 2,
                              save_dir = my_dir,
                              default = default,
                              cluster_method = 'stLearn',
                              cluster_index = cluster_index,
                              platform = platform)
}


melanoma_sub_simulation <- function(i,
                                    List_sce_one,
                                    coordinate_name = c("row","col"),
                                    st_result_path = paste0("./Simulation/Input/melanoma/",
                                                            pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],default,"/",sample_names[i],
                                                            "_stLearn_results.csv"),
                                    
                                    working_directory = paste0("~/simulation/results/",sample_names[i]),
                                    pattern_type = 'BayesSpace_clusters_patch',
                                    my_dir = paste0("~/Desktop/simulation/results/",sample_names[i]),
                                    spatial_probs = c(0.6,0.9),
                                    SV_percent = 0.5,
                                    spat_methods = c('SV_edgeR_NoNorm',
                                                     'SV_edgeR_CPMs','SV_edgeR_CPMs_NoNorm'
                                    ),
                                    spat_methods_params = list(NA,NA, NA),
                                    spat_methods_names = c('SV_edgeR_NoNorm',
                                                           'SV_edgeR_CPMs','SV_edgeR_CPMs_NoNorm'
                                    ),
                                    layer_names = "spatial.cluster",
                                    original_layer_names = "spatial.cluster",
                                    default = T,
                                    reps =1,
                                    cluster_index=NULL, # for multiple patterns, list
                                    platform = "ST",
                                    ...
){
  set.seed(123)
  sce_one = List_sce_one[[i]]
  exprs <- assays(sce_one)$counts
  sce_one <- addPerCellQCMetrics(sce_one)
  sce_one <- logNormCounts(sce_one)
  CD <- colData(sce_one)
  cell_ID = rownames(colData(sce_one))
  colData(sce_one) <- cbind(CD, cell_ID)
  sample_genes = rownames(exprs)
  
  object_subset = sce_one[sample_genes,]
  pattern_type <- pattern_type # 'right_patch', 'circle_patch', 'BayesSpace_clusters_patch'
  if(i==1){combined_cluster <- c(1,4)}
  if(i==2){combined_cluster <- c(2,4)}
  if(i==3){combined_cluster <- c(1,4)}
  if(i==4){combined_cluster <- c(1,3)}
  if(i==5){combined_cluster <- c(1,4)}
  if(i==6){combined_cluster <- c(1,3)}
  if(i==7){combined_cluster <- c(1,2)}
  if(i==8){combined_cluster <- c(1,2)}
  print(pattern_type)
  covariates = cbind.data.frame(x=as.numeric(as.character(colData(object_subset)[[coordinate_name[1]]])),
                                y=as.numeric(as.character(colData(object_subset)[[coordinate_name[2]]])),
                                cell_ID = rownames(colData(object_subset))
  )
  colnames(covariates) <- c("sdimx","sdimy","cell_ID")
  if(pattern_type == 'right_patch'){
    # pattern 1: bottom right stripe
    
    pattern = covariates %>% filter(sdimx > 0 & sdimy < mean(covariates$sdimy))
    pattern_ids = pattern$cell_ID
    pattern_name = 'right_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    #cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    #cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    
    cluster_index = NULL
  } else if (pattern_type == 'circle_patch'){
    # pattern 2: circle pattern
    mx <- round(mean(covariates$sdimx))
    my <- round(mean(covariates$sdimy))
    dist_x <- (covariates$sdimx - mx)^2 
    dist_y <-  (covariates$sdimy - my)^2
    pattern = covariates[dist_x + dist_y < 25,]
    pattern_ids = pattern$cell_ID
    pattern_name = 'circle_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    #cell_meta = as.data.table(colData(object_subset))
    #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    #cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    cluster_index = NULL
  } else if (pattern_type == 'BayesSpace_clusters_patch'){
    # pattern 3: BayesSpace clusters
    colData(sce_one)$spatial.cluster
    #spatPlot2D(gobject = giotto_object_subset, save_plot = T, #cell_color_code = colData(sce_one)$spatial.cluster,
    #             point_size = 0.9, cell_color = as.factor(colData(sce_one)$spatial.cluster), default_save_name  = "BayesSpace_clusters")
    pattern = covariates[as.numeric(colData(sce_one)$spatial.cluster) %in% combined_cluster,]
    pattern_ids = pattern$cell_ID
    pattern_name = 'BayesSpace_clusters_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    #cell_meta = as.data.table(colData(object_subset))
    # cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    #cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
    cluster_index = NULL
  } else if(pattern_type == "MixClusters_patch"){
    pattern_ids <- list()
    pattern_name = 'mixture_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #  cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    pattern_ids <- lapply(cluster_index,
                          function(x) covariates[as.numeric(colData(sce_one)$spatial.cluster) %in% x,'cell_ID'] )
    print(cluster_index)
    print(length(pattern_ids))
    #for(ii in 1:length(pattern_ids)){
    #  x = pattern_ids[[ii]]
    #  print(x)
    #  cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
    
    #}
  } else if(pattern_type == "MixClusters_reverse_patch"){
    pattern_ids <- list()
    pattern_name = 'mixture_reverse_patch'
    pattern_colors = c('in' = 'green', 'out' = 'red')
    cell_meta = as.data.table(colData(object_subset))
    #  cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
    pattern_ids <- lapply(cluster_index,
                          function(x) covariates[as.numeric(colData(sce_one)$spatial.cluster) %notin% x,'cell_ID'] )
    print(cluster_index)
    print(length(pattern_ids))
    # for(ii in 1:length(pattern_ids)){
    #  x = pattern_ids[[ii]]
    #print(x)
    #  cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
    #}
  } 
  print(dim(counts(object_subset)))
  subdir = paste0(my_dir,'/',pattern_name,'/')
  resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
  if(!file.exists(resultdir)) dir.create(path = resultdir, recursive = TRUE)
  save_name = 'probs'
  if(file.exists(paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))) {
    message("Simulation dataset has been created. 
            Skip data simulation part.
            Load the simulation dataset...")
    load(paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))
    if(exists('final_object')) 
      final_object <-  final_object
    else 
      final_object <- combined_object
  }
  
  selected_genes = sample_genes
  my_dir = my_dir
  
  save_dir = my_dir
  runStLearnPatternSimulation(sce_object = object_subset,
                              rawobject = object_subset,
                              coordinate_name = coordinate_name,
                              st_result_path = st_result_path,
                              pattern_name = pattern_name,
                              pattern_cell_ids = pattern_ids,
                              gene_names = selected_genes,
                              spatial_probs = spatial_probs,
                              SV_percent = SV_percent,
                              spat_methods = spat_methods,
                              spat_methods_params = spat_methods_params,
                              spat_methods_names = spat_methods_names,
                              layer_names = layer_names,
                              original_layer_names = original_layer_names,
                              reps = reps,
                              save_plot = T,
                              save_data = TRUE,
                              max_col = 2,
                              save_dir = my_dir,
                              default = default,
                              cluster_method = 'stLearn',
                              cluster_index = cluster_index,
                              platform = platform)
}