library('SPARK')
library(SingleCellExperiment)
library(BiocParallel)
library(scater)
library(limma)
library(edgeR)
library(ggplot2)
library(DEXSeq)
library(DESeq2)
library(reticulate)
library(Matrix)

spe_to_sce <-  function(spe){
  expr = assays(spe)$counts
  # extract metadata from gobject
  metadata = colData(spe)
  sce = SingleCellExperiment(assays=list(counts=expr), 
                             colData=metadata)
  libsizes <- colSums(expr)
  size.factors <- libsizes/mean(libsizes)
  logcounts(sce) <- log2(t(t(expr)/size.factors) + 1)
  return(sce)
}

package_check = function(pkg_name,
                         repository = c('CRAN', 'Bioc', 'github'),
                         github_repo = NULL) {
  
  repository = match.arg(repository, choices = c('CRAN', 'Bioc', 'github'))
  
  if(repository == 'CRAN') {
    
    if(!requireNamespace(pkg_name, quietly = TRUE)) {
      stop("\n package ", pkg_name ," is not yet installed \n",
           "To install: \n",
           "install.packages('",pkg_name,"')",
           call. = FALSE)
    } else {
      return(TRUE)
    }
    
    
  } else if(repository == 'Bioc') {
    
    if(!requireNamespace(pkg_name, quietly = TRUE)) {
      stop("\n package ", pkg_name ," is not yet installed \n",
           "To install: \n",
           "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager');
         BiocManager::install('",pkg_name,"')",
           call. = FALSE)
    } else {
      return(TRUE)
    }
    
  } else if(repository == 'github') {
    
    if(is.null(github_repo)) stop("provide the github repo of package, e.g. 'johndoe/cooltool' ")
    
    if(!requireNamespace(pkg_name, quietly = TRUE)) {
      stop("\n package ", pkg_name ," is not yet installed \n",
           "To install: \n",
           "devtools::install_github('",github_repo,"')",
           call. = FALSE)
    } else {
      return(TRUE)
    }
    
  }
  
}

SV_spark <- function(sce_object,
                     coordinate_name = c("imagerow","imagecol"),
                     num_core = 5,
                     path = '~',
                     sample_name,
                     percentage = 0,
                     min_count = 0,
                     covariates = NULL,
                     ...
){
  set.seed(123)
  
  #   # determine parameter
  #   return_object = match.arg(return_object, c('data.table', 'spark'))
  
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
  write.csv(sobject@res_mtest, paste0(path,"spark_results_",sample_name,".csv"))
  return(sobject)
}

SV_sparkx <- function(sce_object,
                      coordinate_name = c("imagerow","imagecol"),
                      num_core = 5,
                      path = '~',
                      sample_name,
                      ...
){
  set.seed(123)
  
  #   # determine parameter
  #   return_object = match.arg(return_object, c('data.table', 'spark_x'))
  
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
  
  write.csv(sparkX$res_mtest, paste0(path,"sparkx_results_",sample_name,".csv"))
  return(sparkX)
}



suppressMessages(library(MERINGUE))
SV_meringue <- function(sce_object,
                        coordinate_name = c("imagerow","imagecol"),
                        path = './',
                        sample_name = sample,
                        ...
){
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
  start_time <- Sys.time()
  w <- getSpatialNeighbors(pos, filterDist = 10)
  
  ## Identify sigificantly spatially auto-correlated genes
  print("Identify sigificantly spatially auto-correlated genes")
  I <- getSpatialPatterns(mat, w)
  end_time <- Sys.time()
  print(end_time - start_time)
  write.csv(I, paste0(path,"meringue_results_",sample_name,".csv"))
  return(I)
}

SV_edgeR_counts = function(filtered_sce,
                           num_core = 4,
                           layer_names = "layer_guess_reordered",
                           path = '~',
                           sample_name = NULL,
                           save = T,
                           ...) {
  set.seed(123)
  BPPARAM = MulticoreParam(num_core)
  # print message with information #
  message("using 'SV_edgeR' for spatial gene detection. ")
  
  sce = filtered_sce
  layer = colData(sce)[[layer_names]] 
  sce = sce[, !is.na(layer)]
  
  if(length(assayNames(sce)) == 1){
    counts <- assays(sce)$counts
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    logcounts(sce) <- log2(t(t(counts)/size.factors) + 1)}
  
  # store again after removing NA layers:
  layer = as.factor(colData(sce)[[layer_names]])
  # use "cluster" in the design: we look for differences in expression btw clusters:
  design = data.frame(condition = layer)
  design$condition <- droplevels(design$condition)
  y <- DGEList(counts=assays(sce)$counts, 
               genes=rownames(assays(sce)$counts))
  y$samples$lib.size <- colSums(y$counts)
  y <- calcNormFactors(y)
  design_model <- model.matrix(~design$condition)
  rownames(design_model) <- colnames(y)
  y <- estimateDisp(y, design_model, robust=TRUE, BPPARAM = BPPARAM)
  fit <- glmFit(y, design_model, BPPARAM = BPPARAM)
  # Note that glmLRT has conducted a test for the last coefficient in the linear model:
  lrt <- glmLRT(fit, coef = 2:length(unique(layer)))
  # test all 6 layers together (coef = 2:7) or specify an individual one.
  res_edgeR = topTags(lrt, n = Inf)
  
  ## return results ##
  if(save == T){
    write.csv(res_edgeR[[1]], paste0(path,"edgeR_results_",sample_name,".csv"))
  }
  return(list(results = res_edgeR[[1]],
              estimated_y = y))
  
}

library(BayesSpace)
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

#' @title nnSVG
#' @name nnSVG
#' @description Compute spatially expressed genes with nnSVG method
#' @param sce_object SingleCellExperiment object
#' @param num_core number of cores to use
#' @param return_object type of result to return (data.table or SpatialExperiment object)
#' @return data.table with nnSVG spatial genes results or the SpatialExperiment object
#' @details This function is a wrapper for the method implemented in the nnSVG package:
#' @export
SV_nnSVG = function(sce_object,
                    coordinate_name = c("row", "col"),
                    num_core = 4,
                    path = '~',
                    sample_name = NULL,
                    ...) {
  set.seed(123)
  
  #   # determine parameter
  #   return_object = match.arg(return_object, c('data.table', 'nnSVG'))
  
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
  write.csv(rowData(spe), paste0(path,"nnSVG_results_",sample_name,".csv"))
  return(rowData(spe))
}
