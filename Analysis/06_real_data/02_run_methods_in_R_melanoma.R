source("./Analysis/06_real_data/01_methods_in_R_RealData.R")

melanoma_files <- list.files(path = "./DESpace_data/Data/melanoma",
                             pattern = "^mel *.*rda$", full.names = TRUE)
sample_names <- substr(gsub(".rda.*$", "", gsub("^.*/", "", melanoma_files)),1,9)

times_all = NULL
for(ii in c(1:8)){
  times = c()  
  load(melanoma_files[ii])
  sample <- sample_names[ii]
  sce_one = logNormCounts(sce_one)
  print("sparkx")
  set.seed(123)
  a = system.time({
    SV_sparkx(sce=sce_one,
              coordinate_name = c("row","col"),
              num_core = 1,
              path = './DESpace_data/Real/melanoma/',
              sample_name=sample)
  })
  times[1] = a[3]
  names(times[1]) = "SPARK-X"

  print("BayesSpace_edgeR_counts")
  set.seed(123)
  b = system.time({
    sce_bayes <- clusters_BayesSpace(sce_one, q=4,  platform='ST')
    SV_edgeR_counts(filtered_sce = sce_bayes,
                    num_core = 1,
                    covariates = NULL,
                    layer_names = "spatial.cluster",
                    path = './DESpace_data/Real/melanoma/BayesSpace_',
                    sample_name=sample)
  })
  times[2] = b[3]
  names(times[2]) = "BayesSpace_edgeR"

  print("nnSVG")
  set.seed(123)
  c = system.time({
    SV_nnSVG(sce = sce_one,
             coordinate_name = c("row","col"),
             num_core = 1,
             path = './DESpace_data/Real/melanoma/',
             sample_name=sample)
  })
  times[3] = c[3]
  names(times[3]) = "nnSVG"

  print("spark")
  set.seed(123)
  d = system.time({
    SV_spark(sce=sce_one,
             coordinate_name = c("row","col"),
             num_core = 1,
             path = './DESpace_data/Real/melanoma/',
             sample_name=sample)
  })
  times[4] = d[3]
  names(times[4]) = "SPARK"
  
  print("MERINGUE")
  set.seed(123)
  e = system.time({
    SV_meringue(sce=sce_one,
                coordinate_name = c("row","col"),
                path = './DESpace_data/Real/melanoma/',
                sample_name=sample)
  })
  times[5] = e[3]
  names(times[5]) = "MERINGUE"
  
  ## Run stLearn before
  st_results_path <- paste0("./DESpace_data/Real/melanoma/",
                            sample, "_stLearn_results.csv")
  stLearn_results <- read.csv(st_results_path, sep = ',', header = TRUE)
  colnames(stLearn_results) <- c("cell_ID","imagecol","imagerow","stLearn_pca_kmeans")
  
  CD <- colData(sce_one)
  colData(sce_one) <- cbind(CD, stLearn_results$stLearn_pca_kmeans)
  len = length(colnames(colData(sce_one)))
  colnames(colData(sce_one))[len] <- "stLearn_pca_kmeans"
  print("stLearn_edgeR_counts")
  set.seed(123) 
  
  f = system.time({
    SV_edgeR_counts(filtered_sce = sce_one,
                    num_core = 1,
                    layer_names = "stLearn_pca_kmeans",
                    path = './DESpace_data/Real/melanoma/stLearn_',
                    sample_name=sample)
  })
  times[6] = f[3]
  names(times[6]) = "stLearn_edgeR"
  
  
  
  times[7] = sample
  names(times[7]) = "sample_id"
  
  times_all = rbind(times_all,times)
  print(paste0("finish:",sample))
}
colnames(times_all) <- c("SPARK-X","BayesSpace_edgeR","nnSVG",
                         "SPARK","MERINGUE","stLearn_edgeR","sample_id")
write.csv(times_all, paste0("./DESpace_data/Real/melanoma/computational_cost.csv"))


