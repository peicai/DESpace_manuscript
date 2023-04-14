rm(list = ls())
source("./Analysis/03_simulations/02_simulation_function_data_specific.R")
source("./Analysis/04_multiple_samples/01_multiple_samples_simulation_function.R")
## download FULL data:
melanoma_files <- list.files(path = "./DESpace_data/Data/melanoma",
                             pattern = "^mel *.*rda$", full.names = TRUE)
sample_names <- substr(gsub(".rda.*$", "", gsub("^.*/", "", melanoma_files)),1,9)
List_sce_one = list()
# load filtered melanoma data
set.seed(169612)
setwd("./DESpace_data/Data/melanoma")
for(ii in 1:length(sample_names)){
  load(melanoma_files[ii])
  sample <- sample_names[ii]
  layer = colData(sce_one)$spatial.cluster
  List_sce_one[[ii]] = sce_one
}


#pattern_name = 'bottom_pattern'#'MixClusters_reverse_patch'#'bottom_pattern'#'MixClusters_reverse_patch'#'MixClusters_patch'#'mixture_patch'#"Manual_clusters_patch"#"circle_pattern"


Multi_sample_pip <- function(pattern_name, List_sce_one){
  if(pattern_name %in% c('right_pattern','circle_pattern','BayesSpace_clusters_pattern')) {
    SV_percent = 0.33
    spatial_probs = c(0.5,0.8)
    cluster_index <- NULL
    num_clust <- 2
  }else{
    spatial_probs = c(0.5,0.8,0.8,0.8)
    SV_percent = c(0.5,1/6,1/6,1/6)
    cluster_index = list(1,c(2,3),4)
    num_clust <- 4
  }
  default = F
  save_dir = my_dir = "./DESpace_data/Simulation/Output/MultiSample/melanoma/"
  pattern_type = pattern_name
  default = FALSE
  save_name = 'probs'
  cluster_index = NULL
  if(pattern_name %notin% c('right_pattern','circle_pattern','BayesSpace_clusters_pattern')) cluster_index <- list(1,c(2,3),4)
  print(cluster_index)
  sample_names <- c("mel1_rep1", "mel2_rep1", "mel3_rep1", "mel4_rep1")
  layer_names = "spatial.cluster"
  coordinate_name = c("row","col")
  pattern_cell_ids_all <- list()
  rep_i = 1
  num_clust <- ifelse(pattern_name %in% c('right_pattern','circle_pattern','BayesSpace_clusters_pattern'), 2, 4)
  
  ###SvGene
  ###########################################################################
  sce1 = List_sce_one[[1]]
  sce2 = List_sce_one[[3]]
  sce3 = List_sce_one[[5]]
  sce4 = List_sce_one[[7]]
  ## common genes
  a = rownames(counts(sce1))
  b = rownames(counts(sce2))
  c = rownames(counts(sce3))
  d = rownames(counts(sce4))
  CommonGene <- Reduce(intersect, list(a,b,c,d))

  set.seed(123)
  selected_genes = CommonGene
  ## SV genes
  if(length(spatial_probs) == 2){
    SvGene <- sample(as.vector(selected_genes), round(length(selected_genes)*SV_percent))
    folder_path = paste0(save_dir, '/', pattern_name)
    if(!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
    subdir = paste0(save_dir,'/',pattern_name,'/')
    resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
    if(!file.exists(resultdir)) dir.create(resultdir, recursive = TRUE)
    write.table(x = as.matrix(SvGene),
                file = paste0(resultdir, save_name,'_',spatial_probs[2], '_selected_genes.txt'),
                sep = '\t')
  }else {
    SvGene <- split(selected_genes,sample(1:length(SV_percent),length(selected_genes),replace=TRUE,prob=SV_percent))
    folder_path = paste0(save_dir, '/', pattern_name)
    if(!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
    subdir = paste0(save_dir,'/',pattern_name,'/')
    resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
    if(!file.exists(resultdir)) dir.create(resultdir, recursive = TRUE)
    
    names(SvGene) <- paste0(rep(c("pattern_"),length(spatial_probs)), 
                            rep(1:length(spatial_probs)))
    save(SvGene, file = paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_all_genes.rda'))
    write.table(x = as.matrix(unlist(SvGene[-1])),
                file = paste0(resultdir, save_name,'_',spatial_probs[2], '_selected_genes.txt'),
                sep = '\t')
    
  }
  
  ## subset sce
  keep_col <- c("row","col","sizeFactor","sum","detected","total","cluster.init","spatial.cluster")
  sce1 <- sce1[selected_genes,] 
  colData(sce1) <- DataFrame(as.data.frame(colData(sce1)) %>% 
                               dplyr::select(all_of(keep_col)))
  colData(sce1)$"cell_ID" <- rownames(colData(sce1))
  colData(sce1)$"sample_name" <- "mel1_rep1"
  colData(sce1)$row <- as.numeric(as.character(colData(sce1)$row))
  colData(sce1)$col <- as.numeric(as.character(colData(sce1)$col))
  
  sce2 <- sce2[selected_genes,]
  colData(sce2) <- DataFrame(as.data.frame(colData(sce2)) %>% 
                               dplyr::select(all_of(keep_col)))
  colData(sce2)$"cell_ID" <- rownames(colData(sce2))
  colData(sce2)$"sample_name" <- "mel2_rep1"
  colData(sce2)$row <- as.numeric(as.character(colData(sce2)$row))
  colData(sce2)$col <- as.numeric(as.character(colData(sce2)$col))
  
  sce3 <- sce3[selected_genes,]
  colData(sce3) <- DataFrame(as.data.frame(colData(sce3)) %>% 
                               dplyr::select(all_of(keep_col)))
  colData(sce3)$"cell_ID" <- rownames(colData(sce3))
  colData(sce3)$"sample_name" <- "mel3_rep1"
  colData(sce3)$row <- as.numeric(as.character(colData(sce3)$row))
  colData(sce3)$col <- as.numeric(as.character(colData(sce3)$col))
  
  
  sce4 <- sce4[selected_genes,]
  colData(sce4) <- DataFrame(as.data.frame(colData(sce4)) %>% 
                               dplyr::select(all_of(keep_col)))
  colData(sce4)$"cell_ID" <- rownames(colData(sce4))
  colData(sce4)$"sample_name" <- "mel4_rep1"
  colData(sce4)$row <- as.numeric(as.character(colData(sce4)$row))
  colData(sce4)$col <- as.numeric(as.character(colData(sce4)$col))
  
  
  gene_names = selected_genes
  ## input
  sce_list <- list()
  sce_list[[1]] = sce1
  sce_list[[2]] = sce2
  sce_list[[3]] = sce3
  sce_list[[4]] = sce4
  ## Combine into 1 SCE and preprocess
  sce.combined = SingleCellExperiment::cbind(sce1, sce2, sce3, sce4, deparse.level = 1)
  
  
  PatternCellId <- function(j,object_subset,
                            coordinate_name,pattern_type,
                            combined_cluster,
                            cluster_index = NULL){
    covariates = cbind.data.frame(x=as.numeric(as.character(colData(object_subset)[[coordinate_name[1]]])),
                                  y=as.numeric(as.character(colData(object_subset)[[coordinate_name[2]]])),
                                  cell_ID = rownames(colData(object_subset))
    )
    colnames(covariates) <- c("sdimx","sdimy","cell_ID")
    if(pattern_type == 'right_pattern'){
      # pattern 1: bottom right stripe
      
      pattern = covariates %>% filter(sdimx > 0 & sdimy < mean(covariates$sdimy))
      pattern_ids = pattern$cell_ID
      pattern_name = 'right_pattern'
      pattern_colors = c('in' = 'green', 'out' = 'red')
      cell_meta = as.data.table(colData(object_subset))
      #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
      cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
      
      cluster_index = NULL
    } else if (pattern_type == 'circle_pattern'){
      # pattern 2: circle pattern
      mx <- round(mean(covariates$sdimx))
      my <- round(mean(covariates$sdimy))
      dist_x <- (covariates$sdimx - mx)^2 
      dist_y <-  (covariates$sdimy - my)^2
      pattern = covariates[dist_x + dist_y < 15,]
      pattern_ids = pattern$cell_ID
      pattern_name = 'circle_pattern'
      pattern_colors = c('in' = 'green', 'out' = 'red')
      cell_meta = as.data.table(colData(object_subset))
      #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
      cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
      cluster_index = NULL
      
    } else if (pattern_type == 'BayesSpace_clusters_pattern'){
      # pattern 3: BayesSpace clusters
      pattern = covariates[as.numeric(colData(object_subset)$spatial.cluster) %in% combined_cluster,]
      pattern_ids = pattern$cell_ID
      pattern_name = 'BayesSpace_clusters_pattern'
      pattern_colors = c('in' = 'green', 'out' = 'red')
      cell_meta = as.data.table(colData(object_subset))
      #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
      cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
      cluster_index = NULL
    } else if(pattern_type == "MixClusters_pattern"){
      pattern_ids <- list()
      pattern_name = 'mixture_pattern'
      pattern_colors = c('in' = 'green', 'out' = 'red')
      cell_meta = as.data.table(colData(object_subset))
      covariates = as.data.table(covariates)
      #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
      pattern_ids <- lapply(cluster_index,
                            function(x) covariates[as.numeric(colData(object_subset)$spatial.cluster) %in% x,cell_ID] )
      print(cluster_index)
      print(length(pattern_ids))
      for(ii in 1:length(pattern_ids)){
        x = pattern_ids[[ii]]
        print(x)
        cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
        
      }} else if(pattern_type == "MixClusters_reverse_pattern"){
        pattern_ids <- list()
        pattern_name = 'mixture_reverse_pattern'
        pattern_colors = c('in' = 'green', 'out' = 'red')
        cell_meta = as.data.table(colData(object_subset))
        covariates = as.data.table(covariates)
        #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
        pattern_ids <- lapply(cluster_index,
                              function(x) covariates[as.numeric(colData(object_subset)$spatial.cluster) %notin% x,cell_ID] )
        print(cluster_index)
        print(length(pattern_ids))
        for(ii in 1:length(pattern_ids)){
          x = pattern_ids[[ii]]
          #print(x)
          cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
        }} 
    return(list(pattern_ids, cell_meta))
  }
  
  print("Split pattern cell ids...")
  for(j in 1:4){
    if(j==1){combined_cluster <- c(1,4)}
    if(j==2){combined_cluster <- c(1,4)}
    if(j==3){combined_cluster <- c(1,4)}
    if(j==4){combined_cluster <- c(1,2)}
    ## BayesSpace-pattern: choose 1 cluster as in only (maybe melanoma or a smaller one if melanoma doesn’t work);
    #if(j==1){combined_cluster <- c(1)}
    #if(j==3){combined_cluster <- c(2)}
    #if(j==5){combined_cluster <- c(4)}
    #if(j==7){combined_cluster <- c(1)}
    
    pattern_cell_ids_all[[j]] <- PatternCellId(j,object_subset = sce_list[[j]],
                                               coordinate_name = coordinate_name,
                                               pattern_type = pattern_type,
                                               combined_cluster = combined_cluster,
                                               cluster_index = cluster_index)[[1]]
    cell_metadata[[j]] <- PatternCellId(j,object_subset = sce_list[[j]],
                                        coordinate_name = coordinate_name,
                                        pattern_type = pattern_type,
                                        combined_cluster = combined_cluster,
                                        cluster_index = cluster_index)[[2]]
  }
  # ###############################################################################
  print("Start simulation")
  
  library(doParallel)
  detectCores()
  registerDoParallel(3)
  
  runPatternSimulation(sce_list = sce_list, rawobject = NULL,
                       coordinate_name = coordinate_name,
                       SvGene = SvGene,
                       sample_names = sample_names,
                       pattern_name =  pattern_name,
                       save_dir = save_dir,
                       pattern_colors = c('in' = 'green', 'out' = 'red'),
                       pattern_cell_ids_all = pattern_cell_ids_all,
                       gene_names = gene_names,
                       spatial_probs = spatial_probs,
                       reps = 1,
                       SV_percent = SV_percent,
                       layer_names=layer_names,
                       verbose = F,
                       original_layer_names = "spatial.cluster",
                       cluster_index=cluster_index,
                       platform = "ST",
                       default = F)
  
}

pattern_list <- c('right_patch','circle_patch','BayesSpace_clusters_patch','MixClusters_patch','MixClusters_reverse_patch')


library(doParallel)
detectCores()
registerDoParallel(3)


foreach(p = c(1:5)) %dopar% {
  pattern_name <- pattern_list[p] 
  Multi_sample_pip(pattern_name = pattern_name, List_sce_one = List_sce_one)
}

