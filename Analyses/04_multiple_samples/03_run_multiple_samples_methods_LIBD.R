rm(list = ls())
pattern_type = "Manual_clusters_patch"
pattern_name = "Manual_clusters_pattern"
spatial_probs = c(0.5,0.7)
source("./Analysis/04_multiple_samples/01_multiple_samples_additional_methods.R")
source("./Analysis/04_multiple_samples/01_multiple_samples_simulation_function.R")
dir = paste0("./DESpace_data/Simulation/Output/MultiSample/LIBD/",pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
path = paste0(dir,"probs_0.5_0.7_final_object151507.rda")
load(path)
sce1 = final_object; rm(final_object)
combined_cluster <- c(1,4,5)
coordinate_name = c("imagerow","imagecol")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce1)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce1)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce1))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
pattern = covariates[as.numeric(colData(sce1)$layer_guess_reordered) %in% combined_cluster,]
pattern_ids = pattern$cell_ID
pattern_name = 'Manual_clusters_pattern'

cell_meta = as.data.table(colData(sce1))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce1) = DataFrame(cell_meta)
#sce1 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.7_final_object151669.rda")
load(path)
sce2 = final_object; rm(final_object)
combined_cluster <- c(4,7)
coordinate_name = c("imagerow","imagecol")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce2)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce2)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce2))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
pattern = covariates[as.numeric(colData(sce2)$layer_guess_reordered) %in% combined_cluster,]
pattern_ids = pattern$cell_ID
pattern_name = 'Manual_clusters_pattern'

cell_meta = as.data.table(colData(sce2))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce2) = DataFrame(cell_meta)
#sce2 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.7_final_object151673.rda")
load(path)
sce3 = final_object; rm(final_object)
combined_cluster <- c(2,3,4)
coordinate_name = c("imagerow","imagecol")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce3)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce3)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce3))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
pattern = covariates[as.numeric(colData(sce3)$layer_guess_reordered) %in% combined_cluster,]
pattern_ids = pattern$cell_ID
pattern_name = 'Manual_clusters_pattern'

cell_meta = as.data.table(colData(sce3))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce3) = DataFrame(cell_meta)
#sce3 = combined_object; rm(combined_object)
sce.combined = SingleCellExperiment::cbind(sce1, sce2, sce3, deparse.level = 1)

#colData(sce) <- cbind(colData(combined_object), joint_cluster)
sample_id = data.frame(colData(sce.combined)$sample_name)
layer_names = 'JointCluster'

joint_cluster <- colData(sce.combined)[["Manual_clusters_pattern"]]
colData(sce.combined)$JointCluster <- as.factor(joint_cluster)

#colData(sce) <- cbind(colData(combined_object), joint_cluster)
sample_id = data.frame(colData(sce.combined)$sample_name)
layer_names = 'JointCluster'
save(sce.combined, file = paste0(dir, 'sce_original_cluster.rda'))
set.seed(123)
spatial_gene_results <- Single_multi_edgeR(sce_object = sce.combined,
                                           #coordinate_name = c("row", "col"),
                                           layer_names = layer_names,
                                           default = F,
                                           save = T,
                                           edgeR_method = c('multi','single'),
                                           cluster_name = 'Original')
set.seed(123)
spatial_gene_results <- run_LIBD_registration1(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_spatialLIBD_registration_woFilter.rda"))

set.seed(123)
spatial_gene_results <- run_scran_findMarkers(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_scran_findMarkers.rda"))

set.seed(123)
spatial_gene_results <- run_seurat_FindAllMarkers(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_seurat_FindAllMarkers.rda"))

######################### bottom pattern
######################### 
rm(list = ls())
pattern_type = "bottom_patch"
pattern_name = "bottom_pattern"
spatial_probs = c(0.5,0.7)
source("./Analysis/04_multiple_samples/01_multiple_samples_additional_methods.R")
source("./Analysis/04_multiple_samples/01_multiple_samples_simulation_function.R")
dir = paste0("./DESpace_data/Simulation/Output/MultiSample/LIBD/",pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
path = paste0(dir,"probs_0.5_0.7_final_object151507.rda")
load(path)
sce1 = final_object; rm(final_object)
combined_cluster <- c(1,4,5)
coordinate_name = c("imagerow","imagecol")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce1)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce1)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce1))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
pattern = covariates %>% filter(sdimx > 0 & sdimy < mean(covariates$sdimy))
pattern_ids = pattern$cell_ID
pattern_name = 'bottom_pattern'

cell_meta = as.data.table(colData(sce1))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce1) = DataFrame(cell_meta)
#sce1 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.7_final_object151669.rda")
load(path)
sce2 = final_object; rm(final_object)
combined_cluster <- c(4,7)
coordinate_name = c("imagerow","imagecol")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce2)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce2)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce2))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
pattern = covariates %>% filter(sdimx > 0 & sdimy < mean(covariates$sdimy))
pattern_ids = pattern$cell_ID
pattern_name = 'bottom_pattern'

cell_meta = as.data.table(colData(sce2))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce2) = DataFrame(cell_meta)
#sce2 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.7_final_object151673.rda")
load(path)
sce3 = final_object; rm(final_object)
combined_cluster <- c(2,3,4)
coordinate_name = c("imagerow","imagecol")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce3)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce3)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce3))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
pattern = covariates %>% filter(sdimx > 0 & sdimy < mean(covariates$sdimy))
pattern_ids = pattern$cell_ID
pattern_name = 'bottom_pattern'

cell_meta = as.data.table(colData(sce3))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce3) = DataFrame(cell_meta)
#sce3 = combined_object; rm(combined_object)
sce.combined = SingleCellExperiment::cbind(sce1, sce2, sce3, deparse.level = 1)
joint_cluster <- colData(sce.combined)[["bottom_pattern"]]
colData(sce.combined)$JointCluster <- as.factor(joint_cluster)

#colData(sce) <- cbind(colData(combined_object), joint_cluster)
sample_id = data.frame(colData(sce.combined)$sample_name)
layer_names = 'JointCluster'
save(sce.combined, file = paste0(dir, 'sce_original_cluster.rda'))

path = paste0(dir,'sce_original_cluster.rda')
load(path)
#colData(sce) <- cbind(colData(combined_object), joint_cluster)
sample_id = data.frame(colData(sce.combined)$sample_name)
layer_names = 'JointCluster'
set.seed(123)
spatial_gene_results <- Single_multi_edgeR(sce_object = sce.combined,
                                           #coordinate_name = c("row", "col"),
                                           layer_names = layer_names,
                                           default = F,
                                           save = T,
                                           edgeR_method = c('multi','single'),
                                           cluster_name = 'Original')
save(spatial_gene_results, file = paste0(dir,"result_SV_edgeR_counts.rda"))
set.seed(123)
spatial_gene_results <- run_LIBD_registration1(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_spatialLIBD_registration.rda"))

set.seed(123)
spatial_gene_results <- run_scran_findMarkers(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_scran_findMarkers.rda"))

set.seed(123)
spatial_gene_results <- run_seurat_FindAllMarkers(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_seurat_FindAllMarkers.rda"))

######################### circle pattern
######################### 
rm(list = ls())
pattern_type = "circle_patch"
pattern_name = "circle_pattern"
spatial_probs = c(0.5,0.7)
source("./Analysis/04_multiple_samples/01_multiple_samples_additional_methods.R")
source("./Analysis/04_multiple_samples/01_multiple_samples_simulation_function.R")
dir = paste0("./DESpace_data/Simulation/Output/MultiSample/LIBD/",pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
path = paste0(dir,"probs_0.5_0.7_final_object151507.rda")
load(path)
sce1 = final_object; rm(final_object)
combined_cluster <- c(1,4,5)
coordinate_name = c("imagerow","imagecol")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce1)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce1)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce1))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
mx <- round(mean(covariates$sdimx))
my <- round(mean(covariates$sdimy))
dist_x <- (covariates$sdimx - mx)^2 
dist_y <-  (covariates$sdimy - my)^2
pattern = covariates[dist_x + dist_y < 10000,]
pattern_ids = pattern$cell_ID
pattern_name = 'circle_pattern'

cell_meta = as.data.table(colData(sce1))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce1) = DataFrame(cell_meta)
#sce1 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.7_final_object151669.rda")
load(path)
sce2 = final_object; rm(final_object)
combined_cluster <- c(4,7)
coordinate_name = c("imagerow","imagecol")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce2)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce2)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce2))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
mx <- round(mean(covariates$sdimx))
my <- round(mean(covariates$sdimy))
dist_x <- (covariates$sdimx - mx)^2 
dist_y <-  (covariates$sdimy - my)^2
pattern = covariates[dist_x + dist_y < 10000,]
pattern_ids = pattern$cell_ID
pattern_name = 'circle_pattern'

cell_meta = as.data.table(colData(sce2))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce2) = DataFrame(cell_meta)
#sce2 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.7_final_object151673.rda")
load(path)
sce3 = final_object; rm(final_object)
combined_cluster <- c(2,3,4)
coordinate_name = c("imagerow","imagecol")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce3)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce3)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce3))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
mx <- round(mean(covariates$sdimx))
my <- round(mean(covariates$sdimy))
dist_x <- (covariates$sdimx - mx)^2 
dist_y <-  (covariates$sdimy - my)^2
pattern = covariates[dist_x + dist_y < 10000,]
pattern_ids = pattern$cell_ID
pattern_name = 'circle_pattern'

cell_meta = as.data.table(colData(sce3))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce3) = DataFrame(cell_meta)
#sce3 = combined_object; rm(combined_object)
sce.combined = SingleCellExperiment::cbind(sce1, sce2, sce3, deparse.level = 1)
joint_cluster <- colData(sce.combined)[["circle_pattern"]]
colData(sce.combined)$JointCluster <- as.factor(joint_cluster)

#colData(sce) <- cbind(colData(combined_object), joint_cluster)
sample_id = data.frame(colData(sce.combined)$sample_name)
layer_names = 'JointCluster'
save(sce.combined, file = paste0(dir, 'sce_original_cluster.rda'))

path = paste0(dir,'sce_original_cluster.rda')
load(path)
#colData(sce) <- cbind(colData(combined_object), joint_cluster)
sample_id = data.frame(colData(sce.combined)$sample_name)
layer_names = 'JointCluster'
set.seed(123)
spatial_gene_results <- Single_multi_edgeR(sce_object = sce.combined,
                                           #coordinate_name = c("row", "col"),
                                           layer_names = layer_names,
                                           default = F,
                                           save = T,
                                           edgeR_method = c('multi','single'),
                                           cluster_name = 'Original')
save(spatial_gene_results, file = paste0(dir,"result_SV_edgeR_counts.rda"))
set.seed(123)
spatial_gene_results <- run_LIBD_registration1(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_spatialLIBD_registration.rda"))

set.seed(123)
spatial_gene_results <- run_scran_findMarkers(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_scran_findMarkers.rda"))

set.seed(123)
spatial_gene_results <- run_seurat_FindAllMarkers(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_seurat_FindAllMarkers.rda"))

######################### MixClusters_pattern
######################### 
rm(list = ls())
cluster_index = list(3,4,5,6,7)
pattern_type = "mixture_patch"
pattern_name = "MixClusters_pattern"
spatial_probs = c(0.5,0.7)
source("./Analysis/04_multiple_samples/01_multiple_samples_additional_methods.R")
source("./Analysis/04_multiple_samples/01_multiple_samples_simulation_function.R")
dir = paste0("./DESpace_data/Simulation/Output/MultiSample/LIBD/",pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
path = paste0(dir,'sce_original_cluster.rda')
load(path)
#colData(sce) <- cbind(colData(combined_object), joint_cluster)
sample_id = data.frame(colData(sce.combined)$sample_name)
sce.combined$JointCluster <- as.character(sce.combined$layer_guess_reordered)
sce.combined$JointCluster[sce.combined$JointCluster %in% c("Layer1", "Layer2")] <- "OUT"
sce.combined$JointCluster[is.na(sce.combined$JointCluster)] <- "OUT"
table(sce.combined$JointCluster)
layer_names <- "JointCluster"
set.seed(123)
spatial_gene_results <- Single_multi_edgeR(sce_object = sce.combined,
                                           #coordinate_name = c("row", "col"),
                                           layer_names = layer_names,
                                           default = F,
                                           save = T,
                                           edgeR_method = c('multi','single'),
                                           cluster_name = 'Original')
save(spatial_gene_results, file = paste0(dir,"result_SV_edgeR_counts.rda"))

rm(spatial_gene_results)
set.seed(123)
spatial_gene_results <- run_LIBD_registration(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_spatialLIBD_registration.rda"))

rm(spatial_gene_results)
set.seed(123)
spatial_gene_results <- run_scran_findMarkers(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_scran_findMarkers.rda"))

rm(spatial_gene_results)
set.seed(123)
spatial_gene_results <- run_seurat_FindAllMarkers(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_seurat_FindAllMarkers.rda"))

######################### MixClusters_reverse_pattern
######################### 
rm(list = ls())
`%notin%` <- Negate(`%in%`)
cluster_index = list(3,4,5,6,7)
pattern_type = "mixture_reverse_patch"
pattern_name = "MixClusters_reverse_pattern"
spatial_probs = c(0.5,0.7)
source("./Analysis/04_multiple_samples/01_multiple_samples_additional_methods.R")
source("./Analysis/04_multiple_samples/01_multiple_samples_simulation_function.R")
dir = paste0("./DESpace_data/Simulation/Output/MultiSample/LIBD/",pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
path = paste0(dir,'sce_original_cluster.rda')
load(path)
#colData(sce) <- cbind(colData(combined_object), joint_cluster)
sample_id = data.frame(colData(sce.combined)$sample_name)
sce.combined$JointCluster <- as.character(sce.combined$layer_guess_reordered)
sce.combined$JointCluster[sce.combined$JointCluster %in% c("Layer1", "Layer2")] <- "OUT"
sce.combined$JointCluster[is.na(sce.combined$JointCluster)] <- "OUT"
table(sce.combined$JointCluster)
layer_names <- "JointCluster"
set.seed(123)
spatial_gene_results <- Single_multi_edgeR(sce_object = sce.combined,
                                           #coordinate_name = c("row", "col"),
                                           layer_names = layer_names,
                                           default = F,
                                           save = T,
                                           edgeR_method = c('multi','single'),
                                           cluster_name = 'Original')
save(spatial_gene_results, file = paste0(dir,"result_SV_edgeR_counts.rda"))

rm(spatial_gene_results)
source("./Analysis/04_multiple_samples/01_multiple_samples_additional_methods.R")
source("./Analysis/04_multiple_samples/01_multiple_samples_simulation_function.R")
set.seed(123)
spatial_gene_results <- run_LIBD_registration(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_spatialLIBD_registration.rda"))

rm(spatial_gene_results)
set.seed(123)
spatial_gene_results <- run_scran_findMarkers(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_scran_findMarkers.rda"))

rm(spatial_gene_results)
set.seed(123)
spatial_gene_results <- run_seurat_FindAllMarkers(sce.combined)
save(spatial_gene_results, file = paste0(dir,"result_SV_seurat_FindAllMarkers.rda"))