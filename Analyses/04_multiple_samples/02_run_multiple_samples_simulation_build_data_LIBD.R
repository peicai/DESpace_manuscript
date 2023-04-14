rm(list = ls())
source("./Analysis/03_simulations/02_simulation_function_data_specific.R")
source("./Analysis/04_multiple_samples/01_multiple_samples_simulation_function.R")
## download FULL data:
LIBD_files <- list.files(path = "./DESpace_data/Data/LIBD",
                         pattern = "^15 *.*rda$", full.names = TRUE)

sample_names <- substr(gsub(".rda.*$", "", gsub("^.*/", "", LIBD_files)),1,6)
List_sce_one = list()


# load filtered melanoma data
set.seed(169612)
#setwd("~/Desktop/master_thesis/Data/LIBD_BayesSpace")
for(ii in 1:length(sample_names)){
  load(LIBD_files[ii])
  sample <- sample_names[ii]
  layer = colData(sce_one)$layer_guess_reordered
  # remove NAs in layers:
  #sce_one = sce_one[, !is.na(layer)]
  List_sce_one[[ii]] = sce_one
}

#pattern_name = 'bottom_pattern'#'MixClusters_reverse_patch'#'bottom_pattern'#'MixClusters_reverse_patch'#'MixClusters_patch'#'mixture_patch'#"Manual_clusters_patch"#"circle_pattern"


Multi_sample_pip <- function(pattern_name, List_sce_one){
  if(pattern_name %in% c('bottom_pattern','circle_pattern','Manual_clusters_pattern')) {
    SV_percent = 0.33
    spatial_probs = c(0.5,0.8)
    cluster_index <- NULL
    num_clust <- 2
  }else{
    spatial_probs = c(0.5,rep(0.8,5))
    SV_percent = c(0.5,rep(0.1,5))
    cluster_index <- list(3,4,5,6,7)
    num_clust <- 5
  }
  
  save_dir = my_dir = "./DESpace_data/Simulation/Output/MultiSample/LIBD"
  pattern_type = pattern_name
  SubGeneNum = 5000
  save_name = 'probs'
  print(cluster_index)
  sample_names <- c("151507","151669","151673")
  coordinate_name = c("imagerow","imagecol")
  layer_names = "spatial.cluster"
  pattern_cell_ids_all <- list()
  rep_i = 1
  ###SvGene
  ###########################################################################
  #sample_names = sample_names[c(1,5,9)]
  sce1 = List_sce_one[[1]]
  sce2 = List_sce_one[[5]]
  sce3 = List_sce_one[[9]]
  ## common genes
  a = rownames(counts(sce1))
  b = rownames(counts(sce2))
  c = rownames(counts(sce3))
  CommonGene <- Reduce(intersect, list(a,b,c))
  
  ## Randomly select SubGeneNum genes; default = 1000 for the LIBD data
  set.seed(123)
  selected_genes = sample(CommonGene, SubGeneNum)
  
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
  keep_col <- c("barcode","sample_name","tissue","row","col","imagerow","imagecol",
                "layer_guess_reordered","sum","detected","total","sizeFactor")
  sce1 <- sce1[selected_genes,] 
  colData(sce1) <- DataFrame(as.data.frame(colData(sce1)) %>% 
                               dplyr::select(all_of(keep_col)))
  colnames(colData(sce1))[1] <- "cell_ID"
  sce2 <- sce2[selected_genes,]
  colData(sce2) <- DataFrame(as.data.frame(colData(sce2)) %>% 
                               dplyr::select(all_of(keep_col)))
  colnames(colData(sce2))[1] <- "cell_ID"
  sce3 <- sce3[selected_genes,]
  colData(sce3) <- DataFrame(as.data.frame(colData(sce3)) %>% 
                               dplyr::select(all_of(keep_col)))
  colnames(colData(sce3))[1] <- "cell_ID"
  
  gene_names = selected_genes
  ## input
  sce_list <- list()
  sce_list[[1]] = sce1
  sce_list[[2]] = sce2
  sce_list[[3]] = sce3
  ## Combine into 1 SCE and preprocess
  sce.combined = SingleCellExperiment::cbind(sce1, sce2, sce3, deparse.level = 1)
  
  PatternCellId <- function(j,object_subset,
                            coordinate_name,pattern_type,
                            combined_cluster,
                            cluster_index = NULL){
    covariates = cbind.data.frame(x=as.numeric(as.character(colData(object_subset)[[coordinate_name[1]]])),
                                  y=as.numeric(as.character(colData(object_subset)[[coordinate_name[2]]])),
                                  cell_ID = rownames(colData(object_subset))
    )
    colnames(covariates) <- c("sdimx","sdimy","cell_ID")
    if(pattern_type == 'bottom_pattern'){
      # pattern 1: bottom right stripe
      
      pattern = covariates %>% filter(sdimx > 0 & sdimy < mean(covariates$sdimy))
      pattern_ids = pattern$cell_ID
      pattern_name = 'bottom_pattern'
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
      pattern = covariates[dist_x + dist_y < 10000,]
      pattern_ids = pattern$cell_ID
      pattern_name = 'circle_pattern'
      pattern_colors = c('in' = 'green', 'out' = 'red')
      cell_meta = as.data.table(colData(object_subset))
      #cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
      cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
      cluster_index = NULL
    } else if (pattern_type == 'Manual_clusters_pattern'){
      # pattern 3: BayesSpace clusters
      pattern = covariates[as.numeric(colData(object_subset)$layer_guess_reordered) %in% combined_cluster,]
      pattern_ids = pattern$cell_ID
      pattern_name = 'Manual_clusters_pattern'
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
                            function(x) covariates[as.numeric(colData(object_subset)$layer_guess_reordered) %in% x,cell_ID] )
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
                              function(x) covariates[as.numeric(colData(object_subset)$layer_guess_reordered) %notin% x,cell_ID] )
        print(cluster_index)
        print(length(pattern_ids))
        for(ii in 1:length(pattern_ids)){
          x = pattern_ids[[ii]]
          #print(x)
          cell_meta[, (pattern_name) := ifelse(cell_ID %in% x, 'in', 'out')]
        }} 
    return(pattern_ids)
  }
  
  print("Split pattern cell ids...")
  for(j in 1:3){
    if(j == 1){combined_cluster <- c(1,4,5)}
    if(j == 2){combined_cluster <- c(4,7)}
    if(j == 3){combined_cluster <- c(2,3,4)}
    pattern_cell_ids_all[[j]] <- PatternCellId(j,object_subset = sce_list[[j]],
                                               coordinate_name = coordinate_name,
                                               pattern_type = pattern_type,
                                               combined_cluster = combined_cluster,
                                               cluster_index = cluster_index)
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
                       platform = "Visium",
                       default = F)

}

pattern_list <- c('bottom_pattern','circle_pattern','Manual_clusters_pattern',
                  'MixClusters_pattern','MixClusters_reverse_pattern')

library(doParallel)
detectCores()
registerDoParallel(3)


foreach(p = c(1:5)) %dopar% {
  pattern_name <- pattern_list[p] 
  Multi_sample_pip(pattern_name = pattern_name, List_sce_one = List_sce_one)
}


