rm(list=ls())
library(SpatialExperiment)
library(BiocFileCache)
library(scater)
`%notin%` <- Negate(`%in%`)
setwd("./Data/melanoma")

###############################################################################
####### melanoma: do BayesSpace with original dataset (without filtering)######
##############################################################################
sample_names <- c("ST_mel1_rep1","ST_mel1_rep2",
                  "ST_mel2_rep1","ST_mel2_rep2",
                  "ST_mel3_rep1","ST_mel3_rep2",
                  "ST_mel4_rep1","ST_mel4_rep2")
List_sce_one <- list()

for (i in seq_along(sample_names)){
  melanoma_counts <- read.csv(paste0(sample_names[i],"_counts.tsv"), sep = '\t', header = TRUE)
  
  gene_names <- sub(" .*", "", melanoma_counts[,1])
  gene_id <- sub(".* ", "", melanoma_counts[,1])
  rownames(melanoma_counts) <- gene_id
  melanoma_counts <- melanoma_counts[,-1]
  
  coord <- colnames(melanoma_counts)
  raw <- sub(".*X", "", coord)
  row <- as.data.frame(as.factor(sub("x.*", "", raw)))
  col <- as.data.frame(as.factor(sub(".*x", "", coord)))
  df <- cbind(row,col)
  rownames(df) <- raw
  colnames(df) <- c("row","col")
  sce = SingleCellExperiment(assays=list(counts=melanoma_counts), 
                             colData=df)
  sce <- logNormCounts(sce)
  assays(sce)$cpm <- calculateCPM(sce)
  
  colData(sce) = cbind(colData(sce), perCellQCMetrics(sce))
  rowData(sce) = cbind(rowData(sce), perFeatureQCMetrics(sce))
  colData(sce)
  melanoma <- sce
  set.seed(100)
  system.time({
    dec <- scran::modelGeneVar(melanoma)
    top <- scran::getTopHVGs(dec, n = 2000)
    
    set.seed(101)
    melanoma <- scater::runPCA(melanoma, subset_row = top)
    
    ## Add BayesSpace metadata
    melanoma <- spatialPreprocess(melanoma, platform="ST", skip.PCA=TRUE)
    
    q <- 4  # Number of clusters
    d <- 7  # Number of PCs
    
    set.seed(100)
    melanoma <- spatialCluster(melanoma, q=q, d=d, platform='ST',
                               nrep=50000, gamma=2)
  })
  save(melanoma, file = paste0(sample_names[i],".rda"))
}

###############################################################################
########################### QC and filtering ##############################
##############################################################################
rm(list = ls())
library(SpatialExperiment)
library(BiocFileCache)
library(scater)
`%notin%` <- Negate(`%in%`)
setwd("./Data/melanoma")

melanoma_files <- list.files(path = "./",
                             pattern = "^ST_mel *.*rda$", full.names = TRUE)
sample_names <- gsub("ST_","",gsub(".rda.*$", "", gsub("^.*/", "", melanoma_files)))

qc_lib_size_all <- c(250, 500, 1000, 850, 
                     750, 300, 1000, 3500)
qc_detected_all <- c(300, 350, 800, 650, 
                     400, 200, 500, 2000)

for(i in seq_along(sample_names)){
  load(melanoma_files[i])
  sce <- melanoma
  (spe <- SpatialExperiment(
    assays = list(counts = assay(sce)),
    colData = colData(sce)
  ))
  # subset to keep only spots over tissue
  spatialCoords(spe) <- as.matrix(colData(sce)[,1:2])
  spatialCoordsNames(spe) <- c("x_coord","y_coord")
  
  spatialData(spe) <- cbind(spatialData(spe), 
                            in_tissue = rep(sample_names[i],dim(colData(sce))[1]))  
  # calculate per-spot QC metrics and store in colData
  spe <- addPerCellQC(spe,)

  #### sample size:
  sample_id <- sample_names[i]
  spe1 <- spe[, spatialData(spe)$in_tissue == sample_id]
  dim(spe1)
  
  ## Library size
  # select QC threshold for library size
  qc_lib_size <- colData(spe1)$sum < qc_lib_size_all[i]
  colData(spe1)$qc_lib_size <- qc_lib_size
  
  ## Number of expressed features
  # select QC threshold for expressed genes
  qc_detected <- colData(spe1)$detected < qc_detected_all[i]
  colData(spe1)$qc_detected <- qc_detected
  
  ## Remove low-quality spots
  # combined set of discarded spots
  discard <- qc_lib_size | qc_detected
  
  # store in object
  colData(spe1)$discard <- discard
  
  # remove combined set of low-quality spots
  spe1 <- spe1[, !colData(spe1)$discard]
  
  ## Gene level
  ### basic filter: remove undetected genes
  # calculate logcounts (log-transformed normalized counts) and store in object
  spe1 <- logNormCounts(spe1)
  qc_undetected_gene <- rowSums(assays(spe1)$logcounts > 0) > 0
  
  ### basic filter: remove lowly expressed genes: at least xx non-zero cells:
  qc_low_gene <- rowSums(assays(spe1)$counts > 0) >=20
  ## remove combined set of low-quality spots
  
  spe1 <- spe1[qc_low_gene & qc_undetected_gene,]
  dim(spe1)
  
  ## store the filtered dataset
  sce_one = spe1
  path = './Data/melanoma/'
  save(sce_one, file = paste0(path,sample_id,"_melanoma_clean.rda"))
}