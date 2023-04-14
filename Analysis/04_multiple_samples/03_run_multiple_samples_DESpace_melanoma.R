rm(list = ls())
pattern_type = "BayesSpace_clusters_patch"
pattern_name = "BayesSpace_clusters_pattern"
spatial_probs = c(0.5,0.8)
dir = paste0("./DESpace_data/Simulation/Output/MultiSample/melanoma/",pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
path = paste0(dir,"probs_0.5_0.8_final_objectmel1_rep1.rda")
load(path)
sce1 = final_object; rm(final_object)
combined_cluster <- c(1,4)
coordinate_name = c("row","col")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce1)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce1)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce1))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
pattern = covariates[as.numeric(colData(sce1)$spatial.cluster) %in% combined_cluster,]
pattern_ids = pattern$cell_ID
cell_meta = as.data.table(colData(sce1))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce1) = DataFrame(cell_meta)
path = paste0(dir,"probs_0.5_0.8_final_objectmel2_rep1.rda")
load(path)
sce2 = final_object; rm(final_object)
combined_cluster <- c(1,4)
coordinate_name = c("row","col")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce2)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce2)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce2))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
pattern = covariates[as.numeric(colData(sce2)$spatial.cluster) %in% combined_cluster,]
pattern_ids = pattern$cell_ID
pattern_name = 'BayesSpace_clusters_pattern'

cell_meta = as.data.table(colData(sce2))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce2) = DataFrame(cell_meta)
#sce2 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.8_final_objectmel3_rep1.rda")
load(path)
sce3 = final_object; rm(final_object)
combined_cluster <- c(1,4)
coordinate_name = c("row","col")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce3)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce3)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce3))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
pattern = covariates[as.numeric(colData(sce3)$spatial.cluster) %in% combined_cluster,]
pattern_ids = pattern$cell_ID
pattern_name = 'BayesSpace_clusters_pattern'

cell_meta = as.data.table(colData(sce3))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce3) = DataFrame(cell_meta)
path = paste0(dir,"probs_0.5_0.8_final_objectmel4_rep1.rda")
load(path)
sce4 = final_object; rm(final_object)
combined_cluster <- c(1,2)
coordinate_name = c("row","col")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce4)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce4)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce4))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
pattern = covariates[as.numeric(colData(sce4)$spatial.cluster) %in% combined_cluster,]
pattern_ids = pattern$cell_ID
pattern_name = 'BayesSpace_clusters_pattern'

cell_meta = as.data.table(colData(sce4))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce4) = DataFrame(cell_meta)
sce.combined = SingleCellExperiment::cbind(sce1, sce2, sce3, sce4, deparse.level = 1)

path = paste0("./DESpace_data/Simulation/Output/MultiSample/melanoma/",
              pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
joint_cluster <- colData(sce.combined)[["BayesSpace_clusters_pattern"]]
colData(sce.combined)$JointCluster <- as.factor(joint_cluster)

sample_id = data.frame(colData(sce.combined)$sample_name)
layer_names = 'JointCluster'
save(sce.combined, file = paste0(path, 'sce_original_cluster.rda'))
print("SV_edgeR method: ")
source("./Analysis/03_simulations/02_simulation_function_data_specific.R")
source("./Analysis/04_multiple_samples/01_multiple_samples_simulation_function.R")
set.seed(123)
spatial_gene_results <- Single_multi_edgeR(sce_object = sce.combined,
                                           #coordinate_name = c("row", "col"),
                                           layer_names = layer_names,
                                           default = F,
                                           save = T,
                                           edgeR_method = c('multi','single'),
                                           cluster_name = 'Original')
save(spatial_gene_results, file = paste0(path,"result_SV_edgeR_counts.rda"))


######################### right pattern
######################### 
rm(list = ls())
pattern_type = "right_patch"
pattern_name = "right_pattern"
spatial_probs = c(0.5,0.8)
dir = paste0("./DESpace_data/Simulation/Output/MultiSample/melanoma/",pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
path = paste0(dir,"probs_0.5_0.8_final_objectmel1_rep1.rda")
load(path)
sce1 = final_object; rm(final_object)
combined_cluster <- c(1,4)
coordinate_name = c("row","col")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce1)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce1)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce1))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
pattern = covariates %>% filter(sdimx > 0 & sdimy < mean(covariates$sdimy))
pattern_ids = pattern$cell_ID
pattern_name = 'right_pattern'

cell_meta = as.data.table(colData(sce1))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce1) = DataFrame(cell_meta)
#sce1 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.8_final_objectmel2_rep1.rda")
load(path)
sce2 = final_object; rm(final_object)
combined_cluster <- c(1,4)
coordinate_name = c("row","col")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce2)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce2)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce2))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
pattern = covariates %>% filter(sdimx > 0 & sdimy < mean(covariates$sdimy))
pattern_ids = pattern$cell_ID
pattern_name = 'right_pattern'

cell_meta = as.data.table(colData(sce2))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce2) = DataFrame(cell_meta)
#sce2 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.8_final_objectmel3_rep1.rda")
load(path)
sce3 = final_object; rm(final_object)
combined_cluster <- c(1,4)
coordinate_name = c("row","col")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce3)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce3)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce3))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
pattern = covariates %>% filter(sdimx > 0 & sdimy < mean(covariates$sdimy))
pattern_ids = pattern$cell_ID
pattern_name = 'right_pattern'

cell_meta = as.data.table(colData(sce3))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce3) = DataFrame(cell_meta)
#sce3 = combined_object; rm(combined_object)
path = paste0(dir,"probs_0.5_0.8_final_objectmel4_rep1.rda")
load(path)
sce4 = final_object; rm(final_object)
combined_cluster <- c(1,2)
coordinate_name = c("row","col")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce4)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce4)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce4))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
pattern = covariates %>% filter(sdimx > 0 & sdimy < mean(covariates$sdimy))
pattern_ids = pattern$cell_ID
pattern_name = 'right_pattern'

cell_meta = as.data.table(colData(sce4))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce4) = DataFrame(cell_meta)
sce.combined = SingleCellExperiment::cbind(sce1, sce2, sce3, sce4, deparse.level = 1)

path = paste0("./DESpace_data/Simulation/Output/MultiSample/melanoma/",
              pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
joint_cluster <- colData(sce.combined)[["right_pattern"]]
colData(sce.combined)$JointCluster <- as.factor(joint_cluster)

#colData(sce) <- cbind(colData(combined_object), joint_cluster)
sample_id = data.frame(colData(sce.combined)$sample_name)
layer_names = 'JointCluster'
save(sce.combined, file = paste0(path, 'sce_original_cluster.rda'))
print("SV_edgeR method: ")
source("./Analysis/03_simulations/02_simulation_function_data_specific.R")
source("./Analysis/04_multiple_samples/01_multiple_samples_simulation_function.R")
set.seed(123)
spatial_gene_results <- Single_multi_edgeR(sce_object = sce.combined,
                                           #coordinate_name = c("row", "col"),
                                           layer_names = layer_names,
                                           default = F,
                                           save = T,
                                           edgeR_method = c('multi','single'),
                                           cluster_name = 'Original')
save(spatial_gene_results, file = paste0(path,"result_SV_edgeR_counts.rda"))




######################### circle pattern
######################### 
rm(list = ls())
pattern_type = "circle_patch"
pattern_name = "circle_pattern"
spatial_probs = c(0.5,0.8)
dir = paste0("./DESpace_data/Simulation/Output/MultiSample/melanoma/",pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
path = paste0(dir,"probs_0.5_0.8_final_objectmel1_rep1.rda")
load(path)
sce1 = final_object; rm(final_object)
combined_cluster <- c(1,4)
coordinate_name = c("row","col")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce1)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce1)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce1))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
mx <- round(mean(covariates$sdimx))
my <- round(mean(covariates$sdimy))
dist_x <- (covariates$sdimx - mx)^2 
dist_y <-  (covariates$sdimy - my)^2
pattern = covariates[dist_x + dist_y < 15,]
pattern_ids = pattern$cell_ID
pattern_name = 'circle_pattern'

cell_meta = as.data.table(colData(sce1))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce1) = DataFrame(cell_meta)
#sce1 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.8_final_objectmel2_rep1.rda")
load(path)
sce2 = final_object; rm(final_object)
combined_cluster <- c(1,4)
coordinate_name = c("row","col")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce2)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce2)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce2))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
mx <- round(mean(covariates$sdimx))
my <- round(mean(covariates$sdimy))
dist_x <- (covariates$sdimx - mx)^2 
dist_y <-  (covariates$sdimy - my)^2
pattern = covariates[dist_x + dist_y < 15,]
pattern_ids = pattern$cell_ID
pattern_name = 'circle_pattern'

cell_meta = as.data.table(colData(sce2))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce2) = DataFrame(cell_meta)
#sce2 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.8_final_objectmel3_rep1.rda")
load(path)
sce3 = final_object; rm(final_object)
combined_cluster <- c(1,4)
coordinate_name = c("row","col")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce3)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce3)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce3))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
mx <- round(mean(covariates$sdimx))
my <- round(mean(covariates$sdimy))
dist_x <- (covariates$sdimx - mx)^2 
dist_y <-  (covariates$sdimy - my)^2
pattern = covariates[dist_x + dist_y < 15,]
pattern_ids = pattern$cell_ID
pattern_name = 'circle_pattern'

cell_meta = as.data.table(colData(sce3))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce3) = DataFrame(cell_meta)
#sce3 = combined_object; rm(combined_object)
path = paste0(dir,"probs_0.5_0.8_final_objectmel4_rep1.rda")
load(path)
sce4 = final_object; rm(final_object)
combined_cluster <- c(1,2)
coordinate_name = c("row","col")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce4)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(sce4)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(sce4))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
mx <- round(mean(covariates$sdimx))
my <- round(mean(covariates$sdimy))
dist_x <- (covariates$sdimx - mx)^2 
dist_y <-  (covariates$sdimy - my)^2
pattern = covariates[dist_x + dist_y < 15,]
pattern_ids = pattern$cell_ID
pattern_name = 'circle_pattern'

cell_meta = as.data.table(colData(sce4))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_ids, 'in', 'out')]
cluster_index = NULL
colData(sce4) = DataFrame(cell_meta)
sce.combined = SingleCellExperiment::cbind(sce1, sce2, sce3, sce4, deparse.level = 1)

path = paste0("./DESpace_data/Simulation/Output/MultiSample/melanoma/",
              pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
joint_cluster <- colData(sce.combined)[["circle_pattern"]]
colData(sce.combined)$JointCluster <- as.factor(joint_cluster)

#colData(sce) <- cbind(colData(combined_object), joint_cluster)
sample_id = data.frame(colData(sce.combined)$sample_name)
layer_names = 'JointCluster'
save(sce.combined, file = paste0(path, 'sce_original_cluster.rda'))
print("SV_edgeR method: ")
source("./Analysis/03_simulations/02_simulation_function_data_specific.R")
source("./Analysis/04_multiple_samples/01_multiple_samples_simulation_function.R")
set.seed(123)
spatial_gene_results <- Single_multi_edgeR(sce_object = sce.combined,
                                           #coordinate_name = c("row", "col"),
                                           layer_names = layer_names,
                                           default = F,
                                           save = T,
                                           edgeR_method = c('multi','single'),
                                           cluster_name = 'Original')
save(spatial_gene_results, file = paste0(path,"result_SV_edgeR_counts.rda"))




#######
######################### MixClusters_pattern
######################### 
rm(list = ls())
cluster_index = list(1,c(2,3),4)
pattern_type = "mixture_patch"
pattern_name = "MixClusters_pattern"
spatial_probs = c(0.5,0.8)
dir = paste0("./DESpace_data/Simulation/Output/MultiSample/melanoma/",pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
# if(j==1){combined_cluster <- c(1,4)}
# if(j==2){combined_cluster <- c(1,4)}
# if(j==3){combined_cluster <- c(1,4)}
# if(j==4){combined_cluster <- c(1,2)}
path = paste0(dir,"probs_0.5_0.8_final_objectmel1_rep1.rda")
load(path)
sce1 = combined_object; rm(combined_object)
combined_cluster <- c(1,4)
colData(sce1)$spatial.cluster[colData(sce1)$spatial.cluster == 3] <- 2
colData(sce1)$spatial.cluster[colData(sce1)$spatial.cluster == 4] <- 3
pattern_name = 'MixClusters_pattern'
cell_meta = as.data.table(colData(sce1))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := as.factor(cell_meta$spatial.cluster)]
colData(sce1) = DataFrame(cell_meta)
colData(sce1)$sdimx = as.numeric(as.character(colData(sce1)[["row"]]))
colData(sce1)$sdimy = as.numeric(as.character(colData(sce1)[["col"]]))

#sce1 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.8_final_objectmel2_rep1.rda")
load(path)
sce2 = combined_object; rm(combined_object)
combined_cluster <- c(1,4)
colData(sce2)$spatial.cluster[colData(sce2)$spatial.cluster == 3] <- 2
colData(sce2)$spatial.cluster[colData(sce2)$spatial.cluster == 4] <- 3
pattern_name = 'MixClusters_pattern'
cell_meta = as.data.table(colData(sce2))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := as.factor(cell_meta$spatial.cluster)]
colData(sce2) = DataFrame(cell_meta)
colData(sce2)$sdimx = as.numeric(as.character(colData(sce2)[["row"]]))
colData(sce2)$sdimy = as.numeric(as.character(colData(sce2)[["col"]]))
#sce2 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.8_final_objectmel3_rep1.rda")
load(path)
sce3 = combined_object; rm(combined_object)
combined_cluster <- c(1,4)
colData(sce3)$spatial.cluster[colData(sce3)$spatial.cluster == 3] <- 2
colData(sce3)$spatial.cluster[colData(sce3)$spatial.cluster == 4] <- 3
pattern_name = 'MixClusters_pattern'
cell_meta = as.data.table(colData(sce3))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := as.factor(cell_meta$spatial.cluster)]
colData(sce3) = DataFrame(cell_meta)
colData(sce3)$sdimx = as.numeric(as.character(colData(sce3)[["row"]]))
colData(sce3)$sdimy = as.numeric(as.character(colData(sce3)[["col"]]))

#sce3 = combined_object; rm(combined_object)
path = paste0(dir,"probs_0.5_0.8_final_objectmel4_rep1.rda")
load(path)
sce4 = combined_object; rm(combined_object)
combined_cluster <- c(1,2)
colData(sce4)$spatial.cluster[colData(sce4)$spatial.cluster == 3] <- 2
colData(sce4)$spatial.cluster[colData(sce4)$spatial.cluster == 4] <- 3
pattern_name = 'MixClusters_pattern'
cell_meta = as.data.table(colData(sce4))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := as.factor(cell_meta$spatial.cluster)]
colData(sce4) = DataFrame(cell_meta)
colData(sce4)$sdimx = as.numeric(as.character(colData(sce4)[["row"]]))
colData(sce4)$sdimy = as.numeric(as.character(colData(sce4)[["col"]]))

sce.combined = SingleCellExperiment::cbind(sce1, sce2, sce3, sce4, deparse.level = 1)

path = paste0("./DESpace_data/Simulation/Output/MultiSample/melanoma/",
              pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
joint_cluster <- colData(sce.combined)[["MixClusters_pattern"]]
colData(sce.combined)$JointCluster <- as.factor(joint_cluster)

#colData(sce) <- cbind(colData(combined_object), joint_cluster)
sample_id = data.frame(colData(sce.combined)$sample_name)
layer_names = 'JointCluster'
save(sce.combined, file = paste0(path, 'sce_original_cluster.rda'))
print("SV_edgeR method: ")
source("./Analysis/03_simulations/02_simulation_function_data_specific.R")
source("./Analysis/04_multiple_samples/01_multiple_samples_simulation_function.R")
set.seed(123)
spatial_gene_results <- Single_multi_edgeR(sce_object = sce.combined,
                                           #coordinate_name = c("row", "col"),
                                           layer_names = layer_names,
                                           default = F,
                                           save = T,
                                           edgeR_method = c('multi','single'),
                                           cluster_name = 'Original')
save(spatial_gene_results, file = paste0(path,"result_SV_edgeR_counts.rda"))






######################### MixClusters_reverse_pattern
######################### 
rm(list = ls())
cluster_index = list(1,c(2,3),4)
pattern_type = "mixture_reverse_patch"
pattern_name = "MixClusters_reverse_pattern"
spatial_probs = c(0.5,0.8)
dir = paste0("./DESpace_data/Simulation/Output/MultiSample/melanoma/",pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
# if(j==1){combined_cluster <- c(1,4)}
# if(j==2){combined_cluster <- c(1,4)}
# if(j==3){combined_cluster <- c(1,4)}
# if(j==4){combined_cluster <- c(1,2)}
path = paste0(dir,"probs_0.5_0.8_final_objectmel1_rep1.rda")
load(path)
sce1 = combined_object; rm(combined_object)
combined_cluster <- c(1,4)
colData(sce1)$spatial.cluster[colData(sce1)$spatial.cluster == 3] <- 2
colData(sce1)$spatial.cluster[colData(sce1)$spatial.cluster == 4] <- 3
pattern_name = 'MixClusters_reverse_pattern'
cell_meta = as.data.table(colData(sce1))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := as.factor(cell_meta$spatial.cluster)]
colData(sce1) = DataFrame(cell_meta)
colData(sce1)$sdimx = as.numeric(as.character(colData(sce1)[["row"]]))
colData(sce1)$sdimy = as.numeric(as.character(colData(sce1)[["col"]]))

#sce1 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.8_final_objectmel2_rep1.rda")
load(path)
sce2 = combined_object; rm(combined_object)
combined_cluster <- c(1,4)
colData(sce2)$spatial.cluster[colData(sce2)$spatial.cluster == 3] <- 2
colData(sce2)$spatial.cluster[colData(sce2)$spatial.cluster == 4] <- 3
pattern_name = 'MixClusters_reverse_pattern'
cell_meta = as.data.table(colData(sce2))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := as.factor(cell_meta$spatial.cluster)]
colData(sce2) = DataFrame(cell_meta)
colData(sce2)$sdimx = as.numeric(as.character(colData(sce2)[["row"]]))
colData(sce2)$sdimy = as.numeric(as.character(colData(sce2)[["col"]]))
#sce2 = combined_object; rm(combined_object)

path = paste0(dir,"probs_0.5_0.8_final_objectmel3_rep1.rda")
load(path)
sce3 = combined_object; rm(combined_object)
combined_cluster <- c(1,4)
colData(sce3)$spatial.cluster[colData(sce3)$spatial.cluster == 3] <- 2
colData(sce3)$spatial.cluster[colData(sce3)$spatial.cluster == 4] <- 3
pattern_name = 'MixClusters_reverse_pattern'
cell_meta = as.data.table(colData(sce3))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := as.factor(cell_meta$spatial.cluster)]
colData(sce3) = DataFrame(cell_meta)
colData(sce3)$sdimx = as.numeric(as.character(colData(sce3)[["row"]]))
colData(sce3)$sdimy = as.numeric(as.character(colData(sce3)[["col"]]))

#sce3 = combined_object; rm(combined_object)
path = paste0(dir,"probs_0.5_0.8_final_objectmel4_rep1.rda")
load(path)
sce4 = combined_object; rm(combined_object)
combined_cluster <- c(1,2)
colData(sce4)$spatial.cluster[colData(sce4)$spatial.cluster == 3] <- 2
colData(sce4)$spatial.cluster[colData(sce4)$spatial.cluster == 4] <- 3
pattern_name = 'MixClusters_reverse_pattern'
cell_meta = as.data.table(colData(sce4))
#cell_meta = cbind(cell_meta, cell_ID = rownames(colData(object_subset)))
cell_meta[, (pattern_name) := as.factor(cell_meta$spatial.cluster)]
colData(sce4) = DataFrame(cell_meta)
colData(sce4)$sdimx = as.numeric(as.character(colData(sce4)[["row"]]))
colData(sce4)$sdimy = as.numeric(as.character(colData(sce4)[["col"]]))

sce.combined = SingleCellExperiment::cbind(sce1, sce2, sce3, sce4, deparse.level = 1)

path = paste0("./DESpace_data/Simulation/Output/MultiSample/melanoma/",
              pattern_name,"/probs_",spatial_probs[1],"_",spatial_probs[2],"FALSE/")
joint_cluster <- colData(sce.combined)[["MixClusters_reverse_pattern"]]
colData(sce.combined)$JointCluster <- as.factor(joint_cluster)

#colData(sce) <- cbind(colData(combined_object), joint_cluster)
sample_id = data.frame(colData(sce.combined)$sample_name)
layer_names = 'JointCluster'
save(sce.combined, file = paste0(path, 'sce_original_cluster.rda'))
print("SV_edgeR method: ")
source("./Analysis/03_simulations/02_simulation_function_data_specific.R")
source("./Analysis/04_multiple_samples/01_multiple_samples_simulation_function.R")
set.seed(123)
spatial_gene_results <- Single_multi_edgeR(sce_object = sce.combined,
                                           #coordinate_name = c("row", "col"),
                                           layer_names = layer_names,
                                           default = F,
                                           save = T,
                                           edgeR_method = c('multi','single'),
                                           cluster_name = 'Original')
save(spatial_gene_results, file = paste0(path,"result_SV_edgeR_counts.rda"))


#################################################################################################################