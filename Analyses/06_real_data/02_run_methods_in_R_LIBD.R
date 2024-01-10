source("./Analyses/06_real_data/01_methods_in_R_RealData.R")

## LIBD data
# load filtered data
LIBD_files <- list.files(path = "./DESpace_data/Data/LIBD",
                         pattern = "^15 *.*rda$", full.names = TRUE)

sample_names <- substr(gsub(".rda.*$", "", gsub("^.*/", "", LIBD_files)),1,6)
List_sce_one = list()

# load filtered melanoma data
set.seed(169612)
setwd("./DESpace_data/Data/LIBD")
times_all = NULL
for(ii in c(1:12)){
  times = c()  
  load(LIBD_files[ii])
  sample <- sample_names[ii]
  sce_one = logNormCounts(sce_one)
  
  print("spark")
  a = system.time({
    SV_spark(sce=sce_one,
             coordinate_name = c("imagerow","imagecol"),
             num_core = 1,
             path = './DESpace_data/Real/LIBD/',
             sample_name=sample)
  })
  times[1] = a[3]
  names(times[1]) = "SPARK"
  
  print("manual_edgeR_counts")
  b=system.time({
    SV_edgeR_counts(filtered_sce = sce_one,
                    num_core = 1,
                    covariates = NULL,
                    layer_names = "layer_guess_reordered",
                    path = './DESpace_data/Real/LIBD/',
                    sample_name=sample)
  })
  times[2] = b[3]
  names(times[2]) = "Manual_edgeR"
  print("BayesSpace_edgeR_counts")
  tt = system.time({
    sce_bayes <- clustes_BayesSpace(sce_one, q=7, platform='Visium')
  })
  save(sce_bayes, file = paste0(dir, 'BayesSpace_clusters.rda'))
  c = system.time({
    SV_edgeR_counts(filtered_sce = sce_bayes,
                    num_core = 1,
                    covariates = NULL,
                    layer_names = "spatial.cluster",
                    path = './DESpace_data/Real/LIBD/BayesSpace_',
                    sample_name=sample)
  })
  
  times[3] = c[3]
  names(times[3]) = "BayesSpace_edgeR"
  
  print("sparkx")
  d = system.time({
    SV_sparkx(sce=sce_one,
              coordinate_name = c("imagerow","imagecol"),
              num_core = 1,
              path = './DESpace_data/Real/LIBD/',
              sample_name=sample)
  })
  times[4] = d[3]
  names(times[4]) = "SPARK-X"
  
  print("MERINGUE")
  e = system.time({
    SV_meringue(sce=sce_one,
                coordinate_name = c("imagerow","imagecol"),
                path = './DESpace_data/Real/LIBD/',
                sample_name=sample)
  })
  times[5] = e[3]
  names(times[5]) = "MERINGUE"
  
  f = system.time({
    SV_nnSVG(sce = sce_one,
             coordinate_name = c("imagerow","imagecol"),
             num_core = 1,
             path = './DESpace_data/Real/LIBD/',
             sample_name=sample)
  })
  times[6] = f[3]
  names(times[6]) = "nnSVG"
  
  ## Run stLearn before
  st_results_path <- paste0("./DESpace_data/Real/LIBD/",
                            sample, "_stLearn_results.csv")
  stLearn_results <- read.csv(st_results_path, sep = ',', header = TRUE)
  colnames(stLearn_results) <- c("cell_ID","imagecol","imagerow","stLearn_pca_kmeans")
  
  CD <- colData(sce_one)
  colData(sce_one) <- cbind(CD, stLearn_results$stLearn_pca_kmeans)
  len = length(colnames(colData(sce_one)))
  colnames(colData(sce_one))[len] <- "stLearn_pca_kmeans"
  print("stLearn_edgeR_counts")
  set.seed(123) 
  
  g = system.time({
    SV_edgeR_counts(filtered_sce = sce_one,
                    num_core = 1,
                    layer_names = "stLearn_pca_kmeans",
                    path = './DESpace_data/Real/LIBD/stLearn_',
                    sample_name=sample)
  })
  times[7] = g[3]
  names(times[7]) = "stLearn_edgeR"
  
  h=system.time({
    set.seed(123)
    sce = sce_one
    layer = colData(sce)[[layer_names]] 
    sce = sce[, !is.na(layer)]
    layer = as.factor(colData(sce)[[layer_names]])
    scranResults <- scran::findMarkers(
      sce, groups = layer, 
      pval.type = "all",
    )
  })
  save(scranResults, file = paste0(dir, 'Manual_scran_results.rda'))
  times[8] = h[3]
  names(times[8]) = "Manual_scran"
  
  layer_names = "layer_guess_reordered"
  print("manual_seurat")
  j=system.time({
    set.seed(123)
    sce = sce_one
    layer = colData(sce)[[layer_names]] 
    sce = sce[, !is.na(layer)]
    layer = as.factor(colData(sce)[[layer_names]])
    counts <- counts(sce)
    seurat <- CreateSeuratObject(counts = counts)
    t2 <- seurat@meta.data
    t <- data.frame(colData(sce))
    t3 <- data.frame(t,t2)
    seurat@meta.data <- t3
    
    Idents(seurat) <- layer
    # find markers for every cluster compared to all remaining cells
    seuratResults <- Seurat::FindAllMarkers(seurat,
                                            logfc.threshold = -Inf,
                                            min.pct = 0,
                                            test.use = "wilcox",
                                            min.diff.pt = -Inf,
                                            min.cells.feature = 0,
                                            min.cells.group = 0,
                                            return.thresh = 2#, min.pct = 0.25, logfc.threshold = 0.25
    )
  })
  save(seuratResults, file = paste0(dir, 'Manual_seurat_results.rda'))
  times[9] = j[3]
  names(times[9]) = "Manual_seurat"
  
  st_results_path <- paste0("./DESpace_data/Real/LIBD/",
                            sample, "_stLearn_results.csv")
  stLearn_results <- read.csv(st_results_path, sep = ',', header = TRUE)
  colnames(stLearn_results) <- c("cell_ID","imagecol","imagerow","stLearn_pca_kmeans")
  
  CD <- colData(sce_one)
  colData(sce_one) <- cbind(CD, stLearn_results$stLearn_pca_kmeans)
  len = length(colnames(colData(sce_one)))
  colnames(colData(sce_one))[len] <- "stLearn_pca_kmeans"
  print("stLearn_scran")
  layer_names = "stLearn_pca_kmeans"
  
  k = system.time({
    set.seed(123)
    sce = sce_one
    layer = colData(sce)[[layer_names]] 
    sce = sce[, !is.na(layer)]
    layer = as.factor(colData(sce)[[layer_names]])
    scranResults <- scran::findMarkers(
      sce, groups = layer, 
      pval.type = "all",
    )
  })
  
  save(scranResults, file = paste0(dir, 'StLearn_scran_results.rda'))
  times[10] = k[3]
  names(times[10]) = "StLearn_scran"
  
  layer_names = "stLearn_pca_kmeans"
  l=system.time({
    set.seed(123)
    sce = sce_one
    layer = colData(sce)[[layer_names]] 
    sce = sce[, !is.na(layer)]
    layer = as.factor(colData(sce)[[layer_names]])
    counts <- counts(sce)
    seurat <- CreateSeuratObject(counts = counts)
    t2 <- seurat@meta.data
    t <- data.frame(colData(sce))
    t3 <- data.frame(t,t2)
    seurat@meta.data <- t3
    
    Idents(seurat) <- layer
    # find markers for every cluster compared to all remaining cells
    seuratResults <- Seurat::FindAllMarkers(seurat,
                                            logfc.threshold = -Inf,
                                            min.pct = 0,
                                            test.use = "wilcox",
                                            min.diff.pt = -Inf,
                                            min.cells.feature = 0,
                                            min.cells.group = 0,
                                            return.thresh = 2#, min.pct = 0.25, logfc.threshold = 0.25
    )
  })
  save(seuratResults, file = paste0(dir, 'StLearn_seurat_results.rda'))
  times[11] = l[3]
  names(times[11]) = "StLearn_seurat"
  
  layer_names = "spatial.cluster"
  i = system.time({
    set.seed(123)
    sce = sce_bayes
    layer = colData(sce)[[layer_names]] 
    sce = sce[, !is.na(layer)]
    layer = as.factor(colData(sce)[[layer_names]])
    scranResults <- scran::findMarkers(
      sce, groups = layer, 
      pval.type = "all",
    )
  })
  save(scranResults, file = paste0(dir, 'BayesSpace_scran_results.rda'))
  times[12] = i[3]
  names(times[12]) = "BayesSpace_scran"
  
  layer_names = "spatial.cluster"
  q = system.time({
    set.seed(123)
    sce = sce_bayes
    layer = colData(sce)[[layer_names]] 
    sce = sce[, !is.na(layer)]
    layer = as.factor(colData(sce)[[layer_names]])
    counts <- counts(sce)
    seurat <- CreateSeuratObject(counts = counts)
    t2 <- seurat@meta.data
    t <- data.frame(colData(sce))
    t3 <- data.frame(t,t2)
    seurat@meta.data <- t3
    
    Idents(seurat) <- layer
    # find markers for every cluster compared to all remaining cells
    seuratResults <- Seurat::FindAllMarkers(seurat,
                                            logfc.threshold = -Inf,
                                            min.pct = 0,
                                            test.use = "wilcox",
                                            min.diff.pt = -Inf,
                                            min.cells.feature = 0,
                                            min.cells.group = 0,
                                            return.thresh = 2#, min.pct = 0.25, logfc.threshold = 0.25
    )
  })
  save(seuratResults, file = paste0(dir, 'BayesSpace_seurat_results.rda'))
  times[13] = q[3]
  names(times[13]) = "BayesSpace_seurat"
  
  times[14] = sample
  names(times[14]) = "sample_id"
  
  times[15] = tt[3]
  names(times[15]) = "BayesSpace_cluster"
  
  times_all = rbind(times_all,times)
  print(paste0("finish:",sample))
}
colnames(times_all) <- c("SPARK","Manual_edgeR_counts","BayesSpace_edgeR",
                         "SPARK-X","MERINGUE","nnSVG","stLearn_edgeR",
                         "Manual_scran",
                         "Manual_seurat","StLearn_scran",
                         "StLearn_seurat","BayesSpace_scran","BayesSpace_seurat","sample_id",
                         "BayesSpace_cluster")

write.csv(times_all, paste0("./DESpace_data/Real/LIBD/computational_cost.csv"))
