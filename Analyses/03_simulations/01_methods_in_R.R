### Run methods in R: 
### tool name (function name)
### BayesSpace (clusters_BayesSpace), DESpace (SV_edgeR_counts), nnSVG (nnSVG.test), MERINGUE (meringue), SPARK (spark), SPARK-X (spark_x)

suppressMessages({library(SingleCellExperiment)
  library(BiocParallel)
  library(scater)
  library(limma)
  library(scuttle)
  library(edgeR)
  library(ggplot2)
  library(DEXSeq)
  library(DESeq2)
  library(MERINGUE)
  library(BayesSpace)
  library(devtools)  # if not installed: install.packages('devtools')
  library(remotes)  # if not installed: install.packages('remotes')
  # remotes::install_github("RubD/Giotto") 
  # compilation problems (gfortran)?
  # this version does not require C compilation
  # remotes::install_github("RubD/Giotto@cless") 
  library(data.table)
  library(dplyr)
  library(rlang)
  #BiocManager::install("qvalue",force = TRUE)
  library(qvalue)
  library(Seurat)  
  #library(SeuratData)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(patchwork)
  #install.packages('jackstraw')
  #library(Giotto)
  library(SPARK)
  #installGiottoEnvironment()
  #BiocManager::install('MAST')
  library(MAST)
  library(nnSVG)
  library(scran)
})
`%notin%` <- Negate(`%in%`)
addPerCellQCMetrics <- function(x, ...) {
  colData(x) <- cbind(colData(x), scuttle::perCellQCMetrics(x, ...))
  x
}

clusters_BayesSpace <- function(sce, q, platform){
  message("Start clustering (BayesSpace)...")
  set.seed(123)
  dec <- scran::modelGeneVar(sce)
  #if(platform == "Visium"){n = 1000} else{n = 2000}
  n=2000
  top <- scran::getTopHVGs(dec, n=n)
  
  sce <- scater::runPCA(sce, subset_row=top)
  print(platform)
  if(platform %notin% c("ST","Visium")){stop("Execution stopped. The platform should be ST or Visium.")}
  sce <- spatialPreprocess(sce, platform=platform, skip.PCA=TRUE)
  sce <- spatialCluster(sce, q=q,platform=platform,init.method = c("mclust"),
                        model = c("t"), save.chain=TRUE)
  return(sce)
}

SV_edgeR_counts = function(sce_object,
                           min_count = 50,
                           coordinate_name = c("row", "col"),
                           num_core = 4,
                           prior.count = 10, # edgeR:cpm
                           covariates = NULL,
                           return_object = c('data.table'),
                           layer_names = "layer_guess_reordered",
                           original_layer_names = "spatial.cluster",
                           default = T,
                           save = T,
                           dir = "~",
                           pattern_name = 'Manual_clusters_patch',
                           cluster_method='BayesSpace',
                           platform,
                           ...) {
  
  set.seed(123)
  # determine parameter
  # return_object = match.arg(return_object, c('data.table'))
  
  # data.table variables
  genes = FDR = PValue = NULL
  
  # print message with information #
  message("using 'SV_edgeR_counts' for spatial gene/pattern detection. ")
  
  
  ## extract expression values from gobject
  expr = assays(sce_object)$counts
  
  ## extract metadata from gobject
  metadata = colData(sce_object)
  covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce_object)[[coordinate_name[1]]])),
                                y=as.numeric(as.character(colData(sce_object)[[coordinate_name[2]]])))
  colnames(covariates) <- c("image_x","image_y")
  metadata <- cbind(metadata, covariates)
  ## create SingleCellExperiment object
  sce = SingleCellExperiment(assays=list(counts=expr), 
                             colData=metadata)
  
  libsizes <- colSums(expr)
  size.factors <- libsizes/mean(libsizes)
  logcounts(sce) <- log2(t(t(expr)/size.factors) + 1)
  metadata[[original_layer_names]] <- as.factor(metadata[[original_layer_names]])
  print(metadata)
  if(default == T){
    q = 10
    
  }else if(default == F){
    if(pattern_name %in% c('circle_patch','right_patch','bottom_patch','StLearn_clusters_patch','BayesSpace_clusters_patch')){q=2}
    else {q = length(unique(metadata[[original_layer_names]]))}
  }
  
  print(paste0("num clusters: ", q))
  if (q == 1){
    stop("Error: the number of cluster must be at least 2.")
  }
  ## clusters
  print(platform)
  if (layer_names == "spatial.cluster"){
    sce <- clusters_BayesSpace(sce, q = q, platform = platform)}
  
  layer = colData(sce)[[layer_names]] 
  sce = sce[, !is.na(layer)]
  
  if(cluster_method == 'BayesSpace'){save(sce, file = paste0(dir, 'sce_edgeR.rda'))
  }else {save(sce, file = paste0(dir, cluster_method,'_sce_edgeR.rda'))}
  # store again after removing NA layers:
  layer = as.factor(colData(sce)[[layer_names]])
  q = nlevels(layer)
  print(paste0("new num clusters: ", q))
  if (q == 1){
    stop("Error: BayesSpace identifies one cluster.")
  }
  # use "cluster" in the design: we look for differences in expression btw clusters:
  design = data.frame(condition = layer)
  design$condition <- droplevels(design$condition)
  y <- DGEList(counts=assays(sce)$counts, 
               genes=rownames(assays(sce)$counts))
  
  y$samples$lib.size <- colSums(y$counts)
  y <- calcNormFactors(y)
  
  design_model <- model.matrix(~design$condition)
  rownames(design_model) <- colnames(y)
  #is.fullrank(design_model)
  BPPARAM = MulticoreParam(4)
  y <- estimateDisp(y, design_model,# robust=TRUE, 
                    BPPARAM = BPPARAM)
  ##y <- estimateGLMCommonDisp(y, design_model)
  ##y <- estimateGLMTrendedDisp(y, design_model)
  ##y <- estimateGLMTagwiseDisp(y, design_model)
  
  fit <- glmFit(y, design_model, BPPARAM = BPPARAM)
  q = length(unique(colData(sce)[[layer_names]]))
  print(paste0("new num clusters: ", q))
  # Note that glmLRT has conducted a test for the last coefficient in the linear model:
  lrt <- glmLRT(fit, coef = 2:q)
  res_edgeR = topTags(lrt, n = Inf)
  
  ## return results ##
  if(return_object == 'data.table'){
    DT_results = data.table::as.data.table(res_edgeR[[1]])
    gene_names = rownames(res_edgeR[[1]])
    DT_results[, genes := gene_names]
    data.table::setorder(DT_results, FDR, PValue)
    return(DT_results)
  }
}

nnSVG.test = function(sce_object,
                      coordinate_name = c("row", "col"),
                      num_core = 4,
                      return_object = c('data.table', 'nnSVG'),
                      ...) {
  set.seed(123)
  
  # determine parameter
  return_object = match.arg(return_object, c('data.table', 'nnSVG'))
  
  # data.table variables
  genes =  pval = padj = NULL
  
  ## test if nnSVG is installed ##
  package_check(pkg_name = 'nnSVG',
                repository = c('Bioc'))
  #repository = c('github'),
  #github_repo = 'lmweber/nnSVG')
  
  
  # print message with information #
  message("using 'nnSVG' for spatial gene/pattern detection. If used in published research, please cite:
   Weber L.M. et al. (2022), 'nnSVG: scalable identification of spatially variable genes using nearest-neighbor Gaussian processes', bioRxiv")
  
  
  ## extract expression values from gobject
  expr = assays(sce_object)$counts
  
  ## extract metadata from gobject
  metadata = colData(sce_object)
  
  ## extract metadata from gobject
  rowdata = rowData(sce_object)
  
  ## extract coordinates from gobject
  locs = cbind.data.frame(x=as.numeric(as.character(colData(sce_object)[[coordinate_name[1]]])),
                          y=as.numeric(as.character(colData(sce_object)[[coordinate_name[2]]])))
  colnames(locs) <- c("image_x","image_y")
  
  rownames(locs) = colnames(expr)
  
  # Load reducedDim
  reducedDimNames <-
    SingleCellExperiment::reducedDims(sce_object)
  ## create SpatialExperiment object for analysis
  spe <- SpatialExperiment::SpatialExperiment(
    rowData = rowdata,
    colData = cbind(metadata,locs),
    assays = expr,
    reducedDims = reducedDimNames,
    spatialCoordsNames = c("image_x","image_y")
  )
  names(assays(spe)) <- "counts"
  # calculate log-transformed normalized counts using scran package
  # (alternatively could calculate deviance residuals using scry package)
  #set.seed(123)
  qclus <- scran::quickCluster(spe)
  spe <- scran::computeSumFactors(spe, cluster = qclus)
  spe <- logNormCounts(spe)
  
  #set.seed(123)
  
  spe <- nnSVG::nnSVG(spe, BPPARAM = MulticoreParam(num_core))
  
  ## return results ##
  if(return_object == 'nnSVG'){
    return(spe)
  }else if(return_object == 'data.table'){
    DT_results = data.table::as.data.table(rowData(spe))
    gene_names = rownames(rowData(spe))
    DT_results[, genes := gene_names]
    data.table::setorder(DT_results, padj, pval)
    return(DT_results)
  }
}

meringue = function(sce_object,raw_sce,
                    min.reads = 0, 
                    min.lib.size = 0,
                    coordinate_name = c("row", "col"),
                    covariates = NULL,
                    return_object = c('data.table'),
                    ...) {
  
  set.seed(123)
  # determine parameter
  #return_object = match.arg(return_object, c('data.table'))
  
  # data.table variables
  genes =  p.adj = p.value = NULL
  
  ## test if SPARK is installed ##
  package_check(pkg_name = 'MERINGUE',
                repository = c('github'),
                github_repo = 'JEFworks-Lab/MERINGUE')
  
  
  # print message with information #
  message("using 'MERINGUE' for spatial gene/pattern detection. If used in published research, please cite:
  Miller, B., Bambah-Mukku, D., Dulac, C., Zhuang, X. and Fan, J. 'Characterizing spatial gene expression heterogeneity in spatially resolved single-cell transcriptomics data with nonuniform cellular densities.' 
  Genome Research. May 2021.
doi: 10.1101/gr.271288.120 ")
  
  
  ## extract expression values from gobject
  expr = assays(sce_object)$counts
  
  ## extract metadata from gobject
  metadata = colData(sce_object)
  ## extract coordinates from gobject
  covariates = cbind.data.frame(x=as.numeric(as.character(colData(sce_object)[[coordinate_name[1]]])),
                                y=as.numeric(as.character(colData(sce_object)[[coordinate_name[2]]])))
  colnames(covariates) <- c("image_x","image_y")
  
  ## extract coordinates from gobject
  rownames(covariates) = colnames(expr)
  
  ## create SPARK object for analysis and filter out lowly expressed genes
  ## Remove poor datasets and genes
  #counts <- cleanCounts(counts = as.matrix(expr), 
  #                      min.reads = min.reads, 
  #                      min.lib.size = min.lib.size, 
  #                      plot=FALSE,
  #                      verbose=TRUE)
  counts <- as.matrix(expr)
  ## CPM normalize
  mat <- normalizeCounts(counts = counts, 
                         log=FALSE,
                         verbose=TRUE)
  
  
  ## Spatially-aware analysis
  ## Get neighbor-relationships
  pos <- as.data.frame(covariates)[,1:2]
  #print(head(pos))
  
  colnames(pos) <- c("x","y")
  print(length(colnames(counts)))
  print(length( rownames(pos)))
  rownames(pos) = colnames(counts)
  
  print("Get Spatial Neighbors")
  w <- getSpatialNeighbors(pos, filterDist = 10)
  
  ## Identify sigificantly spatially auto-correlated genes
  print("Identify sigificantly spatially auto-correlated genes")
  I <- getSpatialPatterns(mat, w)
  
  ## return results ##
  if(return_object == 'data.table'){
    DT_results = data.table::as.data.table(I)
    gene_names = rownames(I)
    DT_results[, genes := gene_names]
    data.table::setorder(DT_results, p.adj, p.value)
    return(DT_results)
  }
}

spark = function(sce_object,
                 percentage = 0,
                 min_count = 0,
                 coordinate_name = c("row", "col"),
                 num_core = 5,
                 covariates = NULL,
                 return_object = c('data.table', 'spark'),
                 ...) {
  set.seed(123)
  
  # determine parameter
  return_object = match.arg(return_object, c('data.table', 'spark'))
  
  # data.table variables
  genes =  adjusted_pvalue = combined_pvalue = NULL
  
  ## test if SPARK is installed ##
  package_check(pkg_name = 'SPARK',
                repository = c('github'),
                github_repo = 'xzhoulab/SPARK')
  
  
  # print message with information #
  message("using 'SPARK' for spatial gene/pattern detection. If used in published research, please cite:
  Sun, Shiquan, Jiaqiang Zhu, and Xiang Zhou. 'Statistical Analysis of Spatial Expression Pattern for Spatially Resolved Transcriptomic Studies.'
          BioRxiv, October 21, 2019, 810903. https://doi.org/10.1101/810903.")
  
  
  ## extract expression values from gobject
  expr = assays(sce_object)$counts
  
  ## extract metadata from gobject
  metadata = colData(sce_object)
  ## extract coordinates from gobject
  locs = cbind.data.frame(x=as.numeric(as.character(colData(sce_object)[[coordinate_name[1]]])),
                          y=as.numeric(as.character(colData(sce_object)[[coordinate_name[2]]])))
  colnames(locs) <- c("image_x","image_y")
  
  rownames(locs) = colnames(expr)
  
  ## create SPARK object for analysis and filter out lowly expressed genes
  sobject = SPARK::CreateSPARKObject(counts = expr,
                                     location = locs[,1:2],
                                     percentage = percentage,
                                     min_total_counts = min_count)
  
  ## total counts for each cell
  sobject@lib_size = apply(sobject@counts, 2, sum)
  
  ## extract covariates ##
  if(!is.null(covariates)) {
    
    # first filter giotto object based on spark object
    filter_cell_ids = colnames(sobject@counts)
    filter_gene_ids = rownames(sobject@counts)
    #tempgobject = subsetGiotto(gobject, cell_ids = filter_cell_ids, gene_ids = filter_gene_ids)
    tempobject = sce_object[filter_gene_ids, filter_cell_ids]
    metadata = colData(tempobject)
    
    if(!covariates %in% colnames(metadata)) {
      warning(covariates, ' was not found in the cell metadata of the SingleCellExperiment object, will be set to NULL \n')
      covariates = NULL
    } else {
      covariates = metadata[[covariates]]
    }
  }
  
  ## Fit statistical model under null hypothesis
  sobject = SPARK::spark.vc(sobject,
                            covariates = covariates,
                            lib_size = sobject@lib_size,
                            num_core = num_core,
                            verbose = F,
                            ...)
  
  ## test spatially expressed pattern genes
  ## calculating pval
  sobject = SPARK::spark.test(sobject,
                              check_positive = T,
                              verbose = F)
  
  ## return results ##
  if(return_object == 'spark'){
    return(sobject)
  }else if(return_object == 'data.table'){
    DT_results = data.table::as.data.table(sobject@res_mtest)
    gene_names = rownames(sobject@counts)
    DT_results[, genes := gene_names]
    data.table::setorder(DT_results, adjusted_pvalue, combined_pvalue)
    return(DT_results)
  }
}

spark_x = function(sce_object,
                   #percentage = 0.1,
                   #min_count = 10,
                   coordinate_name = c("row", "col"),
                   num_core = 5,
                   #covariates = NULL,
                   return_object = c('data.table', 'spark_x'),
                   ...) {
  set.seed(123)
  
  # determine parameter
  return_object = match.arg(return_object, c('data.table', 'spark_x'))
  
  # data.table variables
  genes =  adjusted_pvalue = combined_pvalue = NULL
  
  ## test if SPARK is installed ##
  package_check(pkg_name = 'SPARK',
                repository = c('github'),
                github_repo = 'xzhoulab/SPARK')
  
  
  # print message with information #
  message("using 'SPARK' for spatial gene/pattern detection. If used in published research, please cite:
  Sun, Shiquan, Jiaqiang Zhu, and Xiang Zhou. 'Statistical Analysis of Spatial Expression Pattern for Spatially Resolved Transcriptomic Studies.'
          BioRxiv, October 21, 2019, 810903. https://doi.org/10.1101/810903.")
  
  
  ## extract expression values from gobject
  expr = assays(sce_object)$counts
  
  ## extract metadata from gobject
  metadata = colData(sce_object)
  ## extract coordinates from gobject
  locs = cbind.data.frame(x=as.numeric(as.character(colData(sce_object)[[coordinate_name[1]]])),
                          y=as.numeric(as.character(colData(sce_object)[[coordinate_name[2]]])))
  colnames(locs) <- c("image_x","image_y")
  
  rownames(locs) = colnames(expr)
  location <- as.matrix(locs[,1:2])
  
  expr = as.matrix(expr)
  ## Analyze the data with SPARK-X
  sparkX <- sparkx(expr,location,numCores=num_core,option="mixture")
  ## return results ##
  if(return_object == 'sparkx'){
    return(sparkX)
  }else if(return_object == 'data.table'){
    DT_results = data.table::as.data.table(sparkX$res_mtest)
    gene_names = rownames(sparkX$res_mtest)
    DT_results[, genes := gene_names]
    data.table::setorder(DT_results, adjustedPval, combinedPval)
    return(DT_results)
  }
}
