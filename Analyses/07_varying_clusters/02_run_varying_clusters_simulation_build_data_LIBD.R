rm(list = ls())
source("./Analysis/07_varying_clusters/01_varying_clusters_simulation_function.R")
library(spdep)
LIBD_files <- list.files(path = "./DESpace_data/Data/LIBD",
                         pattern = "^15 *.*rda$", full.names = TRUE)
sample_names <- substr(gsub(".rda.*$", "", gsub("^.*/", "", LIBD_files)),1,6)
List_sce_one = list()

# load filtered melanoma data
set.seed(169612)
for(ii in 1:length(sample_names)){
  load(LIBD_files[ii])
  sample <- sample_names[ii]
  layer = colData(sce_one)$layer_guess_reordered
  List_sce_one[[ii]] = sce_one
} 
i=1
pattern_type = 'split_mixture_patch'
my_dir = paste0("./Simulation/Output/",sample_names[i])
spat_methods = c('SV_edgeR_counts')
spat_methods_params = list(NA)
spat_methods_names = c('SV_edgeR_counts')
layer_names = "spatial.cluster"
original_layer_names = "layer_guess_reordered"
default = F
reps =1
coordinate_name = c("imagerow","imagecol")
platform = "Visium"

set.seed(123)
sce_one = List_sce_one[[i]]
sce_one = sce_one[, !is.na(sce_one$layer_guess_reordered)]
exprs <- assays(sce_one)$counts
sce_one <- addPerCellQCMetrics(sce_one)
sce_one <- logNormCounts(sce_one)
CD <- colData(sce_one)
cell_ID = rownames(colData(sce_one))
colData(sce_one) <- cbind(CD, cell_ID)
sample_genes = rownames(exprs)

object_subset = sce_one[sample_genes,]
covariates = cbind.data.frame(x=as.numeric(as.character(colData(object_subset)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(object_subset)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(object_subset))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
CD <- as.data.frame(CD)

## split the original 7 clusters (manual annotations) into 2, 4, 6, 8, 10 and 12 clusters

CD$pattern_ids <- with(CD, sign(imagerow + 1.1*imagecol - 635))
CD$pattern_ids2 <- with(CD, sign(imagerow + 1.1*imagecol - 480))
CD$pattern_ids3 <- with(CD, sign(imagerow + 1.1*imagecol - 700))

for(n_clus in c(2,4,6,8,10,12)){
  CD$layer_23 <- with(CD, sign(imagerow - 0.7*imagecol + 120))
  CD$layer_new <- as.character(CD$layer_guess_reordered)
  CD$layer_new <- ifelse(CD$layer_guess_reordered %in% c("Layer2", "Layer3") &
                           CD$layer_23 == 1,
                         paste0(CD$layer_guess_reordered, "_2nd"),
                         as.character(CD$layer_guess_reordered))
  if(n_clus == 12){
    CD$layer_new[CD$layer_new %in% c("Layer3_2nd", "Layer2_2nd")] <- "Layer32_2nd"
    CD$layer_new[CD$layer_new %in% c("Layer3", "Layer2")] <- "Layer32"
    CD$layer_new[CD$layer_new %in% c("Layer5", "Layer4")] <- "Layer45"
    CD$layer_new[CD$layer_new %in% c("WM", "Layer6")] <- "WM6"
    CD$pattern <- CD$layer_new
    CD$pattern <- ifelse(CD$pattern %in% c("Layer32_2nd", "Layer45",
                                           "Layer1") &
                           CD$pattern_ids3 == 1,
                         as.character(paste0(CD$pattern, "1")),
                         CD$pattern)
    CD$pattern <- ifelse(CD$pattern %in% c("Layer32_2nd", "Layer45",
                                           "Layer1") &
                           CD$pattern_ids3 == -1&
                           CD$pattern_ids2 == 1,
                         as.character(paste0(CD$pattern, "0")),
                         CD$pattern)
    CD$pattern <- ifelse(CD$pattern %in% c("Layer32_2nd", "Layer45",
                                           "Layer1") &
                           CD$pattern_ids2 == 1,
                         as.character(paste0(CD$pattern, "-1")),
                         CD$pattern)
    CD$pattern <- ifelse(CD$pattern %in% c("WM6") &
                           CD$pattern_ids3 == 1,
                         "WM6-1",
                         CD$pattern)
  }else if(n_clus == 10){
    CD$layer_new[CD$layer_new %in% c("Layer3_2nd", "Layer2_2nd")] <- "Layer32_2nd"
    CD$layer_new[CD$layer_new %in% c("Layer3", "Layer2")] <- "Layer32"
    CD$layer_new[CD$layer_new %in% c("Layer5", "Layer4")] <- "Layer45"
    CD$layer_new[CD$layer_new %in% c("WM", "Layer6")] <- "WM6"
    CD$pattern <- paste0(CD$layer_new, CD$pattern_ids)
    CD$pattern[CD$pattern %in% c("WM61", "WM6-1")] <- "WM6"
    CD$pattern[CD$pattern %in% c("Layer32-1", "Layer321")] <- "Layer32"
    CD$pattern[CD$pattern %in% c("Layer32_2nd-1", "Layer32_2nd1")] <- "Layer32_2nd"
    CD$pattern <- ifelse(CD$pattern %in% c("Layer32_2nd") &
                           CD$pattern_ids3 == 1,
                         "Layer32_2nd1",
                         CD$pattern)
    CD$pattern <- ifelse(CD$pattern %in% c("Layer32_2nd") &
                           CD$pattern_ids3 == -1&
                           CD$pattern_ids2 == 1,
                         "Layer32_2nd-1",
                         CD$pattern)
    CD$pattern <- ifelse(CD$pattern %in% c("Layer32_2nd") &
                           CD$pattern_ids2 == 1,
                         "Layer32_2nd1",
                         CD$pattern)
    CD$pattern <- ifelse(CD$pattern %in% c("WM6") &
                           CD$pattern_ids3 == 1,
                         "WM6-1",
                         CD$pattern)
  }else if(n_clus == 8){
    CD$layer_new[CD$layer_new %in% c("Layer3_2nd", "Layer2_2nd")] <- "Layer32_2nd"
    CD$layer_new[CD$layer_new %in% c("Layer3", "Layer2")] <- "Layer32"
    CD$layer_new[CD$layer_new %in% c("Layer5", "Layer4")] <- "Layer45"
    CD$layer_new[CD$layer_new %in% c("WM", "Layer6")] <- "WM6"
    CD$pattern <- paste0(CD$layer_new, CD$pattern_ids)
    
  }else if(n_clus == 6){
    CD$pattern <- paste0(CD$layer_new, CD$pattern_ids)
    CD$pattern[CD$pattern %in% c("WM1","Layer61", "Layer51","Layer41")] <- "A"
    CD$pattern[CD$pattern %in% c("WM-1","Layer6-1", "Layer5-1","Layer4-1")] <- "B"
    CD$pattern[CD$pattern %in% c("Layer3_2nd-1", "Layer2_2nd-1")] <- "C"
    CD$pattern[CD$pattern %in% c("Layer3_2nd1", "Layer2_2nd1")] <- "D"
    CD$pattern[CD$pattern %in% c("Layer31", "Layer21","Layer11")] <- "E"
    CD$pattern[CD$pattern %in% c("Layer3-1", "Layer2-1","Layer1-1")] <- "F"
  }else if(n_clus == 4){
    CD$pattern <- paste0(CD$layer_guess_reordered, CD$pattern_ids)
    CD$pattern[CD$pattern %in% c("WM-1", "Layer6-1", "Layer5-1","Layer4-1")] <- "A"
    CD$pattern[CD$pattern %in% c("WM1", "Layer61", "Layer51","Layer41")] <- "B"
    CD$pattern[CD$pattern %in% c("Layer1-1", "Layer2-1","Layer3-1")] <- "C"
    CD$pattern[CD$pattern %in% c("Layer11", "Layer21","Layer31")] <- "D"
  }else if(n_clus == 2){
    CD$pattern <- as.character(CD$pattern_ids)
  }


cluster_index <- unique(CD$pattern)
pattern_ids <- list()
pattern_name = paste0('split_mixture_patch_', length(cluster_index))
pattern_colors = c('in' = 'green', 'out' = 'red')
pattern_ids <- lapply(cluster_index,
                      function(x) covariates[CD$pattern %in% x,'cell_ID'] )
spatial_probs = c(0.5,rep(0.9,length(pattern_ids)))
SV_percent = c(0.5,rep(round(0.5/length(pattern_ids),4),length(pattern_ids)))
subdir = paste0(my_dir,'/',pattern_name,'/')
resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
resultdir1 <- paste0(my_dir,'/split_mixture_patch/',
                     'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
if(!file.exists(resultdir)) dir.create(path = resultdir, recursive = TRUE)
save_name = 'probs'
selected_genes <- c(read.table(file = paste0(resultdir1, save_name,'_0.9_selected_genes.txt'),
                               sep = '\t'))[[1]]
NonSVGs <- c(read.table(file = paste0(resultdir1, save_name,'_0.9_NonSVGs.txt'),
                        sep = '\t'))[[1]]
SVGs <- c(read.table(file = paste0(resultdir1, save_name,'_0.9_SVGs.txt'),
                     sep = '\t'))[[1]]
colData(object_subset) <- DataFrame(CD)
mini_object = object_subset[selected_genes,]
save_dir = my_dir
colData(mini_object)$"cell_ID" <- rownames(colData(mini_object))

## start simulating
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
                     platform = platform,
                     SvGene = SVGs,
                     NonSvGene = NonSVGs)
}

################## original clusters 
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
  
  spatial_gene_results <- SV_edgeR_counts(simobject,
                                          layer_names = "pattern",
                                          original_layer_names = "pattern",
                                          dir = resultdir,           
                                          pattern_name = pattern_name,
                                          cluster_method='original',
                                          platform = "Visium",
                                          return_object = 'data.table')
  spatial_gene_results[, method := 'Original']
  save(spatial_gene_results, file = paste0(resultdir,"result_SV_edgeR_original_clusters.rda"))
  print(i)
}
