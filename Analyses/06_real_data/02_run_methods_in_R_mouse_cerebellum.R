source("./Analyses/06_real_data/01_methods_in_R_RealData.R")

load("./DESpace_data/Data/mouse_cerebellum/cerebellum_filtered.rda")
set.seed(123)
sce_one = logNormCounts(sce)
exprs <- assays(sce_one)$counts
sce_one <- addPerCellQCMetrics(sce_one)
sce_one <- logNormCounts(sce_one)
CD <- colData(sce_one)
cell_ID = rownames(colData(sce_one))
colData(sce_one) <- cbind(CD, cell_ID)
st_result_path <- paste0('./DESpace_data/Data/mouse_cerebellum/100_stLearn_results.csv')
stLearn_results <- read.csv(st_result_path, sep = ',', header = TRUE)
colnames(stLearn_results) <- c("cell_ID","imagecol","imagerow","stLearn_pca_kmeans")

CD <- colData(sce_one)
colData(sce_one) <- cbind(CD, stLearn_results$stLearn_pca_kmeans)
sce_one <- scuttle::logNormCounts(sce_one)
times_all = NULL
List_sce_one = list()
`%notin%` <- Negate(`%in%`)
sample <- "SlideSeq2"
times = c()  

print("stLearn_edgeR_counts")
set.seed(123)
a=system.time({
  SV_edgeR_counts(filtered_sce = sce_one,
                  num_core = 1,
                  covariates = NULL,
                  layer_names = "stLearn_results$stLearn_pca_kmeans",
                  path = './DESpace_data/Real/mouse_cerebellum/StLearn_',
                  sample_name=sample)
})
times[1] = a[3]
names(times[1]) = "StLearn_edgeR"

print("BayesSpace_edgeR_counts")
set.seed(123)
b = system.time({
  dec <- scran::modelGeneVar(sce_one)
  top <- scran::getTopHVGs(dec, n = 2000)
  sce <- scater::runPCA(sce_one, subset_row=top)
  sce <- spatialPreprocess(sce, n.PCs=7, log.normalize=FALSE)
  q = 5
  rm(dec)
  sce <- spatialCluster(sce, q=q, 
                        gamma=5, save.chain=FALSE)
  SV_edgeR_counts(filtered_sce = sce,
                  num_core = 1,
                  covariates = NULL,
                  layer_names = "spatial.cluster",
                  path = './DESpace_data/Real/mouse_cerebellum/BayesSpace_',
                  sample_name=sample)
})
rm(sce)

times[2] = b[3]
names(times[2]) = "BayesSpace_edgeR"


print("sparkx")
set.seed(123)
c = system.time({
  SV_sparkx(sce=sce_one,
            coordinate_name = c("row","col"),
            num_core = 1,
            path = './DESpace_data/Real/mouse_cerebellum/',
            sample_name=sample)
})
times[3] = c[3]
names(times[3]) = "SPARK-X"

print("MERINGUE")
set.seed(123)
d = system.time({
  SV_meringue(sce=sce_one,
              coordinate_name = c("row","col"),
              path = './DESpace_data/Real/mouse_cerebellum/',
              sample_name=sample)
})
times[4] = d[3]
names(times[4]) = "MERINGUE"
print("spark")
library(reticulate)
e = system.time({
  SV_spark(sce=sce_one,
           coordinate_name = c("row","col"),
           num_core = 1,
           path = './DESpace_data/Real/mouse_cerebellum/',
           sample_name=sample)
})
times[5] = e[3]
names(times[5]) = "SPARK"

print("nnSVG")
set.seed(123)
f = system.time({
  SV_nnSVG(sce = sce_one,
           coordinate_name = c("row","col"),
           num_core = 1,
           path = './DESpace_data/Real/mouse_cerebellum/',
           sample_name=sample)
})
times[6] = f[3]
names(times[6]) = "nnSVG"

times[7] = sample
names(times[7]) = "sample_id"

times_all = rbind(times_all,times)
print(paste0("finish:",sample))

colnames(times_all) <- c("StLearn_edgeR","BayesSpace_edgeR",
                         "SPARK-X","MERINGUE","SPARK","nnSVG","sample_id")
write.csv(times_all, paste0("./DESpace_data/Real/mouse_cerebellum/computational_cost.csv"))
