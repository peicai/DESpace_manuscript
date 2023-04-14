######################################################################################
############################### packages ##############################################
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scater)
  #library(harmony)
  library(BayesSpace)
  library(BiocFileCache)
  library(edgeR)
  library(BiocParallel)
  library(doParallel)
  library(parallel)
  library(MASS)
  library(foreach)
  library(ggpubr)
  library(tidyverse)
})

`%notin%` <- Negate(`%in%`)

################################ Processing ###########################################
#' Convert SpatialExperiment to SingleCellExperiment
#' 
#' @param spe SpatialExperiment
#'   
#' @return sce SingleCellExperiment
#' 
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
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

#####################################################################################
############################### Gene level test #####################################
#' @title SV_edgeR
#' @name SV_edgeR
#' @description Creates a known spatial pattern for selected genes one-by-one and runs the different spatial gene detection tests
#' @param sce_object SingleCellExperiment
#' @param layer_col Column name of spatial clusters in colData(sce)
#' @param sample_col Column name of sample ids in colData(sce)
#' @param replicates Single sample or multi-sample test
#' @param num_core Number of cores
#' @return List of results (gene_results and estimated_y): edgeR test results and estimated dispersion
#' @export

SV_edgeR = function(sce_object,
                    layer_col = 'joint_cluster',
                    sample_col = NULL,
                    replicates = F,
                    num_core = 4,
                    verbose = TRUE,
                    ...) {

 # data.table variables
  genes = FDR = PValue = NULL
  
  # print message with information #
  message("using 'SV_edgeR' for spatial gene/pattern detection. ")
  
  # if replicates == T, run a multi-sample edgeR test
  # results: a list; test results (DT_results_multi) and estimated dispersion (estimated_y_multi)
  if(!is.null(sample_col)){
    num_sample <- length(unique(colData(sce_object)[[sample_col]]))
  }else(num_sample = NULL)
  
  if(replicates == T){
    if( sample_col %notin% colnames(colData(sce_object))){
      sprintf("'sample_col' %s  not in colData(sce)", sample_col)
      return(NULL)
    }
    stopifnot(num_sample > 1)
    #nlevels(as.factor(colData(sce_object)[[sample_col]])) > 1
    if( layer_col %notin% colnames(colData(sce_object))){
      sprintf("'layer_col' %s  not in colData(sce)", layer_col)
      return(NULL)
    }

    message("multi-sample test")
    if(verbose){
      list[gene_results, estimated_y, glmLRT, glmFit] = .multi_edgeR_test(sce = sce_object, 
                                                                          layer_col = layer_col, 
                                                                          sample_col = sample_col, 
                                                                          num_core = num_core,
                                                                          verbose = verbose)
      }else{
      list[gene_results, estimated_y] = .multi_edgeR_test(sce = sce_object, 
                                                          layer_col = layer_col, 
                                                          sample_col = sample_col, 
                                                          num_core = num_core,
                                                          verbose = verbose)
      }
    
    
    DT_results = gene_results
    data.table::setorder(DT_results, FDR, PValue)
    DT_results_multi = DT_results
    estimated_y_multi = estimated_y
  }
  
  
  # if replicates == F, run single-sample edgeR tests
  # if sce_object is a list of sce objects, run single-sample edgeR test for each sample via the function ".multi_single_edgeR_test"
  # results: a list; test results (DT_results_single) and estimated dispersion (estimated_y_single)

  if(replicates == F && !is.null(num_sample)){
    if( layer_col %notin% colnames(colData(sce_object))){
      sprintf("'layer_col' %s  not in colData(sce)", layer_col)
      return(NULL)
    }
    if( sample_col %notin% colnames(colData(sce_object))){
      sprintf("'sample_col' %s  not in colData(sce)", sample_col)
      return(NULL)
    }
    message("using 'single' for spatial gene/pattern detection. ")
    
    sample_names = unique(colData(sce_object)[[sample_col]])
    n_sample = length(sample_names)
    results = data.frame()
    results = lapply(1:n_sample,.multi_single_edgeR_test,
                     sce_object=sce_object,
                     sample_names=sample_names,
                     layer_col=layer_col,
                     verbose = verbose
    )
    
    results_all = do.call("rbind", results)
    DT_results_single = results_all$DT_results1
    estimated_y_single = results_all$estimated_y
  }
  
  # if replicates == F, run a single-sample edgeR test
  # if sce_object is a sce object, run a single-sample edgeR test 
  # results: a list; test results (DT_results_single) and estimated dispersion (estimated_y_single)
  
  if(replicates == F && is.null(num_sample)){
    if( layer_col %notin% colnames(colData(sce_object))){
      sprintf("'layer_col' %s  not in colData(sce)", layer_col)
      return(NULL)
    }

    message("single sample test")
    # if "logcounts" not included in assayNames(sce), calculate logcounts
    if(length(assayNames(sce_object)) == 1){
      counts <- assays(sce_object)$counts
      libsizes <- colSums(counts)
      size.factors <- libsizes/mean(libsizes)
      logcounts(sce_object) <- log2(t(t(counts)/size.factors) + 1)}
    
    if(verbose){
      list[gene_results, estimated_y, glmLRT, glmFit] = .single_edgeR_test(sce = sce_object, 
                                                                           layer_col = layer_col, 
                                                                           num_core = num_core,
                                                                           verbose = verbose)
    }else{
      list[gene_results, estimated_y] = .single_edgeR_test(sce = sce_object, 
                                                           layer_col = layer_col, 
                                                           num_core = num_core,
                                                           verbose = verbose)
    }
    
    DT_results_single = gene_results
    estimated_y_single = estimated_y
  }
  
  # if replicates == T, gene_results estimated_y -> from multi-sample results
  if(replicates == T){
    gene_results = DT_results_multi
    estimated_y = estimated_y_multi
  }
  # if replicates == F, gene_results, estimated_y -> from single-sample results
  if(replicates == F){
    gene_results = DT_results_single
    estimated_y = estimated_y_single
  }
  
  if(verbose){
    results_list <- list(gene_results, estimated_y, glmLRT, glmFit)
  }else{
    results_list <- list(gene_results, estimated_y)
  }

  ## return results ##
  #if(save == T){
  #  write.csv(results_list, paste0(path,"multi_edgeR_results.csv"))
  #}
  return(results_list)
}

#' Run single sample edgeR test for multiple samples together 
#' 
#' @param i i-th sample
#' @param sce_object SingleCellExperiment with sample_name in colData
#' @param sample_names Vector of all sample names (e.g., sample_names <- c("151507","151508","151509"))
#' @param layer_col Column name of spatial clusters in colData(sce)
#' @param num_core Number of cores
#' @return List of results (DT_results1 and estimated_y): edgeR test results and estimated dispersion
#' 
#' @keywords internal
.multi_single_edgeR_test <- function(i,
                                     sce_object,
                                     sample_names,
                                     layer_col,
                                     num_core = 4,
                                     verbose = TRUE){
  sce = subset(sce_object,,sample_name == sample_names[i])
  if(verbose){
    list[gene_results, estimated_y, glmLRT, glmFit] = .single_edgeR_test(sce, layer_col, num_core, verbose)
  }else{
    list[gene_results, estimated_y] = .single_edgeR_test(sce, layer_col, num_core, verbose)
  }
   
  DT_results1 = data.table::as.data.table(gene_results)
  DT_results1[, sample := sample_names[i]]
  data.table::setorder(DT_results1, FDR, PValue)
  if(verbose){
    return(list(DT_results1 = DT_results1,
              estimated_y = estimated_y,
              glmLRT = lrt,
              glmFit = fit))
  }else{
    return(list(DT_results1 = DT_results1,
                estimated_y = estimated_y))
  }
}
#' Run single sample edgeR test for one sample  
#' 
#' @param sce SingleCellExperiment
#' @param layer_col Column name of spatial clusters in colData(sce)
#' @param num_core Number of cores
#' @return List of results (gene_results and estimated_y): edgeR test results and estimated dispersion
#'
#' @keywords internal
.single_edgeR_test <- function(sce, 
                               layer_col, 
                               num_core,
                               verbose = TRUE){
  layer = as.factor(colData(sce)[[layer_col]] )
  sce = sce[, !is.na(layer)]
  design = data.frame(condition = layer)
  
  design$condition <- droplevels(design$condition)
  y <- DGEList(counts=assays(sce)$counts, 
               genes=rownames(assays(sce)$counts))
  y$samples$lib.size <- colSums(y$counts)
  y <- calcNormFactors(y)
  design_model <- model.matrix(~design$condition)
  rownames(design_model) <- colnames(y)
   
  y <- estimateDisp(y, design_model, robust=TRUE)
  fit <- glmFit(y, design_model)
  q = nlevels(factor(colData(sce)[[layer_col]]))
  
  lrt <- glmLRT(fit, coef = 2:q)
  res_edgeR = topTags(lrt, n = Inf)
  results = data.table::as.data.table(res_edgeR[[1]][c("genes","LR","logCPM", "PValue", "FDR")])
  colnames(results)[1] <- "gene_id"
  if(verbose){
    return(list(gene_results = results,
                estimated_y = y,
                glmLRT = lrt,
                glmFit = fit))
  }else{
    return(list(gene_results = results,
                estimated_y = y))
  }
  
}

#' Run multi-sample edgeR test for multiple samples 
#' 
#' @param sce SingleCellExperiment 
#' @param layer_col Column name of spatial clusters in colData(sce)
#' @param sample_col column name of sample ids in colData(sce) 
#' @param num_core Number of cores
#' @return List of results (gene_results and estimated_y): edgeR test results and estimated dispersion
#'
#'@keywords internal
.multi_edgeR_test <- function(sce, 
                              layer_col, 
                              sample_col, 
                              num_core,
                              verbose = TRUE){
  layer = as.factor(colData(sce)[[layer_col]] )
  sce = sce[, !is.na(layer)]
  design = data.frame(condition = factor(colData(sce)[[layer_col]]),
                      sample_id = factor(colData(sce)[[sample_col]]))
  
  design$condition <- droplevels(design$condition)
  y <- DGEList(counts=assays(sce)$counts, 
               genes=rownames(assays(sce)$counts))
  y$samples$lib.size <- colSums(y$counts)
  y <- calcNormFactors(y)
  design_model <- model.matrix(~design$condition + design$sample_id)
  rownames(design_model) <- colnames(y)
   
  y <- estimateDisp(y, design_model, robust=TRUE)
  fit <- glmFit(y, design_model)
  q = nlevels(factor(colData(sce)[[layer_col]]))
  
  lrt <- glmLRT(fit, coef = 2:q)
  res_edgeR = topTags(lrt, n = Inf)
  results = data.table::as.data.table(res_edgeR[[1]][c("genes","LR", "logCPM", "PValue", "FDR")])
  colnames(results)[1] <- "gene_id"
  if(verbose){
    return(list(gene_results = results,
                estimated_y = y,
                glmLRT = lrt,
                glmFit = fit))
  }else{
    return(list(gene_results = results,
                estimated_y = y))
  }
  
}

############################### Individual cluster ###############################################

#' @title individual_test
#' @name individual_test
#' @description Creates a known spatial pattern for selected genes one-by-one and runs the different spatial gene detection tests
#' @param sce_object SingleCellExperiment
#' @param layer_col Column name of spatial clusters in colData(sce)
#' @param sample_col Column name of sample ids in colData(sce)
#' @param edgeR_y Pre-estimated dispersion; if it's null, compute dispersion
#' @param replicates Single sample or multi-sample test
#' @param num_core Number of cores
#' @return List of results (gene_results and estimated_y): edgeR test results and estimated dispersion
#' @export

individual_test <- function(sce_object,
                            layer_col = "layer_guess_reordered",
                            sample_col = "sample_name",
                            edgeR_y=NULL, 
                            replicates = FALSE, 
                            num_core=4,
                            ... ){
  if( layer_col %notin% colnames(colData(sce_object))){
    sprintf("'layer_col' %s  not in colData(sce)", layer_col)
    return(NULL)
  }

  layer = colData(sce_object)[[layer_col]] 
  layer  <- droplevels(layer)
  sce_object = sce_object[, !is.na(layer)]
  # Re-label the cluster (e.g., cluster 1 vs. rest)

  print("Pre-processing")
  layer_list <- foreach (i = 1: nlevels(as.factor(layer))) %dopar% {
      cluster = levels(as.factor(layer))[i]
      .cluster_label(sce_object, cluster_list = cluster, layer_col = layer_col)
    
  }

   
  print("Start modeling")
  
  # if edgeR_y != NULL -> use dispersion taken from gene-level test (with all clusters)
  if(!is.null(edgeR_y)){
    if(replicates == T){
      if( sample_col %notin% colnames(colData(sce_object))){
        sprintf("'sample_col' %s  not in colData(sce)", sample_col)
        return(NULL)
      }
      sample_id = factor(colData(sce_object)[[sample_col]])
      result_list <- lapply(layer_list, .layer_test1_multi, 
                            y = edgeR_y, sample_id = sample_id)
      #single_cluster_results = lapply(result_list, `[[`, 1)
    }else{
      
      result_list <- lapply(layer_list, .layer_test1, y = edgeR_y)
      #single_cluster_results = lapply(result_list, `[[`, 1)
    }
  }else if(is.null(edgeR_y)){
    # if edgeR_y == NULL -> dispersion computed each time
    print("Input data (estimated dispersion) is missing.")
    print("Re-compute dispersion.")
    edgeR_y <- DGEList(counts=assays(sce_object)$counts, 
                 genes=rownames(assays(sce_object)$counts))
    edgeR_y$samples$lib.size <- colSums(edgeR_y$counts)
    edgeR_y <- calcNormFactors(edgeR_y)
    
    if(replicates == T){
      sample_id = factor(colData(sce_object)[[sample_col]])
      result_list <- lapply(layer_list, .layer_test2_multi, 
                            y = edgeR_y, sample_id = sample_id)
      #single_cluster_results = lapply(result_list, `[[`, 1)
    }else{
      result_list <- lapply(layer_list, .layer_test2, y = edgeR_y)
      #single_cluster_results = lapply(result_list, `[[`, 1)
    }
    }
  #head(single_cluster_results)

  names(result_list) <- levels(as.factor(layer))
  message("Returning results")
  return(cluster_results = result_list)
  
}


############################### layer test ################
#' Re-label clusters 
#' 
#' @param sce SingleCellExperiment 
#' @param cluster_list The cluster name  
#' @param layer_col Column name of spatial clusters in colData(sce)
#' @return Vector of cluster labels (only keep \code{cluster_list} cluster, set other clusters as "Other")
#'
#'@keywords internal
.cluster_label <- function(sce, 
                           cluster_list, 
                           layer_col = "layer_guess_reordered"){
  #print(cluster)
  metadata = as.data.frame(colData(sce))
  layer_rename = c()
  layer_rename = as.character(metadata[[layer_col]])
  layer_rename[layer_rename != cluster_list] <- 'Other'
  layer_rename = as.factor(layer_rename)
  one_layer = layer_rename
  return(one_layer)
}

#' Layer test 1
#' We use the dispersion estimate (estimateDisp) from the gene-level model, and simply run the cluster-specific test with all clusters
#' Much faster (estimateDisp is not re-run); but it's slightly wrong: use the dispersion estimated from the wrong design (not the cluster-specific design we use of testing)
#' 
#' @param layers Vector of cluster labels 
#' @param y Estimated dispersion  
#' @return Table of edgeR test results (genes, LR, logFC, PValue and FDR)
#'
#'@keywords internal
.layer_test1 <- function(layers, 
                         y){
  one_layer = layers
  design_model = model.matrix(~one_layer)
  fit <- glmFit(y, design_model)
  lrt <- glmLRT(fit, coef = 2)
  results <- topTags(lrt, n = Inf)
  results = data.table::as.data.table(results[[1]][c("genes","LR", "logCPM", "logFC", "PValue", "FDR")])
  colnames(results)[1] <- "gene_id"
  return(results)
}

#' Layer test 1 [multi-sample version]
#' 
#' @param layers Vector of cluster labels 
#' @param sample_id Vector of sample ids (\code{factor(colData(sce_object)[[sample_col]])})
#' @param y Estimated dispersion 
#' @return Table of edgeR test results (genes, LR, logFC, PValue and FDR)
#'
#'@keywords internal
.layer_test1_multi <- function(layers, 
                               sample_id, 
                               y){
  one_layer = layers
  design_model = model.matrix(~one_layer+sample_id)
  fit <- glmFit(y, design_model)
  lrt <- glmLRT(fit, coef = 2)
  results <- topTags(lrt, n = Inf)
  results = data.table::as.data.table(results[[1]][c("genes", "LR", "logCPM", "logFC", "PValue", "FDR")])
  colnames(results)[1] <- "gene_id"
  return(results)
}

#' Layer test 2
#' We re-compute the dispersion estimate (estimateDisp) for each cluster-specific test with the new design
#' More accurate, but much much slower
#' 
#' @param layers Vector of cluster labels 
#' @param y Estimated dispersion 
#' @return Table of edgeR test results (genes, LR, logFC, PValue and FDR)
#'
#'@keywords internal
.layer_test2 <- function(layers, 
                         y){
  # layers = factor(one_layer, levels  = c('Other', 'Layer3'))
  one_layer = layers #== layer
  design_model = model.matrix(~one_layer)
  limma::is.fullrank(design_model)
  rownames(design_model) <- colnames(y)
  y <- estimateDisp(y, robust=TRUE, design = design_model)
  fit <- glmFit(y, design_model)
  lrt <- glmLRT(fit, coef = 2)
  results <- topTags(lrt, n = Inf)
  results = data.table::as.data.table(results[[1]][c("genes", "LR", "logCPM", "logFC", "PValue", "FDR")])
  colnames(results)[1] <- "gene_id"
  return(results)
}
#' Layer test 2 [multi-sample version]
#' 
#' @param layers Vector of cluster labels 
#' @param sample_id Vector of sample ids (\code{factor(colData(sce_object)[[sample_col]])})
#' @param y Estimated dispersion  
#' @return Table of edgeR test results (genes, LR, logFC, PValue and FDR)
#'
#'@keywords internal
.layer_test2_multi <- function(layers, 
                               sample_id, 
                               y){
  # layers = factor(one_layer, levels  = c('Other', 'Layer3'))
  one_layer = layers #== layer
  design_model = model.matrix(~one_layer)
  limma::is.fullrank(design_model)
  rownames(design_model) <- colnames(y)
  y <- estimateDisp(y, robust=TRUE, design = design_model)
  fit <- glmFit(y, design_model)
  lrt <- glmLRT(fit, coef = 2)
  results <- topTags(lrt, n = Inf)
  results = data.table::as.data.table(results[[1]][c("genes", "LR", "logCPM", "logFC", "PValue", "FDR")])
  colnames(results)[1] <- "gene_id"
  return(results)
}


############################## top results  ###############################################################

#' @title top_results
#' @name top_results
#' @description Filter significant results. \code{top_results} returns the significant results obtained via \code{\link{SV_edgeR}} and \code{\link{individual_test}}
#' @param gene_results  A \code{\linkS4class{data.frame}} with results as returned from \code{\link{SV_edgeR}}
#' @param cluster_results A list of \code{\linkS4class{data.frame}} with results as returned from \code{\link{individual_test}}
#' @param cluster Specify the cluster in interest. A character indicating the cluster(s) whose results have to be returned
#' Results from all clusters are returned by default ("NULL")
#' @param select A character indicating which results should be extracted
#' Results of adjusted p-value and log fold-change ("both", default choice) are returned. Only adjusted p-value (FDR) or log fold-change (logFC) is returned
#' @param high_low A character indicating whether to return: all results ("both" or "BOTH"), only highly abundant cluster SVGs ("high" or "HIGH") or lowly abundant cluster SVGs ("low" or "LOW")
#' In `cluster_results`, a logFC > 0 indicates highly abundant cluster1 SVGs (compared to rest clusters); while a logFC < 0 indicates lowly abundant cluster1 SVGs (compared to rest clusters)
#' @return A \code{\linkS4class{data.frame}} object or a list of \code{\linkS4class{data.frame}} with results
#' Given \code{gene_results}, results combining gene-and cluster-level are returned
#' If \code{high_low} is specified, a list of \code{\linkS4class{data.frame}} containing highly and/or lowly abundant cluster SVG results is returned
#' If \code{cluster} is specified, a subset \code{\linkS4class{data.frame}} only containing SVG results of the cluster is returned
#' Otherwise, a \code{\linkS4class{data.frame}} containing selected results (by default `select = "both"`; `FDR` and `logFC`) of all clusters is returned
#' @export
top_results = function(gene_results = NULL,
                       cluster_results, 
                       cluster = NULL,
                       select = "both",
                       high_low = NULL){
  stopifnot(
    is.list(cluster_results),
    #is.numeric(cluster_order), length(cluster_order) == 1L,
    is.character(select), length(select) == 1L
  )
  
  cluster_results <- lapply(cluster_results, as.data.frame)

  # If there is gene level results (stored in gene_results), return merged results (merge gene and cluster level results)
  if(!is.null(gene_results)){
    gene_results <- as.data.frame(gene_results)
    merge_results <- .merge_results(gene_results, cluster_results, select)
    merge_results <- as.data.frame(merge_results) %>% arrange(gene_FDR)
    return(res = merge_results)
  }
  
  if(!is.null(high_low)){
    stopifnot(is.character(high_low), length(high_low) == 1L)
    if( !(high_low %in% c("both", "BOTH", "high", "HIGH", "low", "LOW")) ){
      message("'high_low' should be one of: 'both', 'BOTH', 'high', 'HIGH', 'low', 'LOW'")
      return(NULL)
    }
    }
  if(!is.null(cluster)){
    stopifnot(is.character(cluster), length(cluster) == 1L)
  }
  #if( !(select %in% c("p_adj", "FDR", "p_val", "logFC")) ){
  #  message("'select' should be one of: 'p_adj', 'FDR', 'p_val', 'logFC'")
  #  return(NULL)
  #}
  #if(select == 'p_adj'){select = "FDR"}
  #if(select == 'p_val'){select = "PValue"}
  
  # Extract and combine results of the selected column ('select') from each sample
  com_results <- .getValueRes(cluster_results = cluster_results, select = 'FDR')
  # Top genes -> for each gene/row, min p-values/FDR; for each gene, return a vector of layer names
  #rank_results_list <- apply(com_results,1,function(x) which(x == min(x)))
  #rank_results <- as.data.frame(paste(lapply(rank_results_list,names),sep = ","))
  sel_layer_min <- apply(com_results,1,which.min)
  rank_results<-as.data.frame(colnames(com_results)[sel_layer_min])
  colnames(rank_results) <- "Layer"
  # For each gene, extract LR, PValue and FDR based on the layer name (i.e., the layer has min FDR)
  rank_results <- cbind(rank_results, 
                        lapply(1:nrow(com_results), function(x) cluster_results[[sel_layer_min[x]]][x,]) %>% 
                                                                                             bind_rows()) 
  # Rename columns; 3:7 -> "LR", "logCPM", "logFC", "PValue", "FDR"
  colnames(rank_results)[3:7] <- paste0("Layer_", colnames(rank_results)[3:7])
  colnames(com_results) <- paste0(names(cluster_results), "_FDR")
  # Add the "gene id" column
  com_results['gene_id'] <- rownames(com_results)
  rank_results['gene_id'] <- rownames(com_results)
  rownames(rank_results) <- rank_results$gene_id
  # Add the results of identified layers (i.e., rank_results) into com_results
  com_results <- merge(rank_results, com_results, by = "gene_id") 
  # Check genes: p-values/FDR for all layers are same
  # gene_all_layer <- as.data.frame(com_results==apply(com_results,1, min))
  # gene_all_layer <- rownames(gene_all_layer %>% filter(rowSums(gene_all_layer) != 1))
  rownames(com_results) <- com_results$gene_id
  # Sort results based on FDR
  com_results <- as.data.frame(com_results) %>% arrange(Layer_FDR)
  # Sort results based on the specific cluster i; return top genes for cluster i
  if(!is.null(cluster)){
    gene_id <- com_results[order(com_results[, paste0((cluster), "_FDR")], decreasing = FALSE), 'gene_id']
    # Reorder rows of rank_results based on FDR of the specific cluster i
    top_results_cluster <- rank_results[match(gene_id, rank_results$gene_id),]
    # Subset Layer == cluster i
    top_genes_cluster <- subset(top_results_cluster, Layer == (cluster))
  }
  
  # Sort results based on the specific cluster i; return top genes for cluster i
  if(!is.null(cluster) & !is.null(high_low)){
    sel_FC <- .getValueRes(cluster_results = cluster_results, select = 'logFC')
    
    if(high_low %in% c("HIGH","high")){
      FC_high_genes <- rownames(sel_FC[sel_FC[,(cluster)] >= 0,])
      high_genes <- as.data.frame(top_genes_cluster[top_genes_cluste$gene_id %in% FC_high_genes,])
      low_genes <- NULL
    }
    if(high_low %in% c("LOW","low")){
      FC_low_genes <- rownames(sel_FC[sel_FC[,(cluster)] < 0,])
      low_genes <- as.data.frame(top_genes_cluster[top_genes_cluster$gene_id %in% FC_low_genes,])
      high_genes <- NULL
    }
    if(high_low %in% c("BOTH","both")){
      FC_high_genes <- rownames(sel_FC[sel_FC[,(cluster)] >= 0,])
      high_genes <- as.data.frame(top_genes_cluster[top_genes_cluster$gene_id %in% FC_high_genes,])
      
      FC_low_genes <- rownames(sel_FC[sel_FC[,(cluster)] < 0,])
      low_genes <- as.data.frame(top_genes_cluster[top_genes_cluster$gene_id %in% FC_low_genes,])
    }
    
  }
  if(is.null(high_low)){
    if(is.null(cluster)){
      return(res = com_results)
    }else{
      return(top_genes = top_genes_cluster)
    }
  }else{
    return(list(#res = com_results, 
                #top_genes = top_genes_cluster,
                high_genes = high_genes,
                low_genes = low_genes
                ))
    }
}

#' Merge gene-and cluster-level results
#' 
#' @param gene_results  A \code{\linkS4class{data.frame}} with results as returned from \code{\link{SV_edgeR}}
#' @param cluster_results A list of \code{\linkS4class{data.frame}} with results as returned from \code{\link{individual_test}}
#' @param select A character indicating which results should be extracted
#' @return A table of adjusted p-value and/or log fold-change ("both", default choice) for each cluster is returned
#'
#'@keywords internal
.merge_results <- function(gene_results,
                           cluster_results, 
                           select = "both"){
  stopifnot(
    is.list(cluster_results),
    is.data.frame(gene_results),
    is.character(select), length(select) == 1L
  )
  
  if( !(select %in% c("p_adj", "FDR", "logFC","both")) ){
    message("'select' should be one of: 'p_adj', 'FDR', 'logFC','both'")
    return(NULL)
  }
  #if(dim(gene_results)[2] == 5){
 #   colnames(gene_results) <- c("gene_id", "gene_id", "gene_LR","gene_Pvalue","gene_FDR")
  #}else(colnames(gene_results) <- c("gene_id", "gene_LR","gene_Pvalue","gene_FDR"))
  colnames(gene_results) <- c("gene_id", "gene_LR", "gene_logCPM", "gene_Pvalue","gene_FDR")
  
  if(select != 'both'){
    com_results <- .getValueRes(cluster_results = cluster_results, select = select)
    colnames(com_results) <- paste0(colnames(com_results), "_", select)
    com_results['gene_id'] <- rownames(com_results)
    # Merge gene-level and cluster-level results
    merge_results <- merge(gene_results, com_results, by = "gene_id")
  }
  
  if(select == 'both'){
    com_results1 <- .getValueRes(cluster_results = cluster_results, select = c('FDR'))
    colnames(com_results1) <- paste0(colnames(com_results1), "_", 'FDR')
    
    com_results2 <- .getValueRes(cluster_results = cluster_results, select = c('logFC'))
    colnames(com_results2) <- paste0(colnames(com_results2), "_", 'logFC')
    
    com_results <- cbind(com_results1, com_results2)
    com_results['gene_id'] <- rownames(com_results)
    # Merge gene-level and cluster-level results
    merge_results <- merge(gene_results, com_results, by = "gene_id")
  }
  return(merge_results)
}

#' Extract 'FDR' (select = "FDR') results for each sample from a list (\code{cluster_results}; obtained via \code{\link{individual_test}})
#' 
#' @param cluster_results A list of \code{\linkS4class{data.frame}} with results as returned from \code{\link{individual_test}}
#' @param select A character indicating which results should be extracted
#' @return A \code{\linkS4class{data.frame}} containing selected results ("FDR", default choice) for each cluster is returned
#'
#'@keywords internal
.getValueRes <- function(cluster_results = single_cluster_results,
                        select = "FDR"){
  ll <- lapply(1:length(cluster_results), function(i){
    data.frame(gene_id = cluster_results[[i]]$gene_id, 
               value  = get(select, cluster_results[[i]]))})
  names(ll) <- names(cluster_results)
 
  ll2 <- lapply(ll, function(LL) lapply(LL, `[`, order(LL$gene_id)) )
  genesID <- ll2[[1]]$gene_id
  ll2 <-lapply(ll2, "[",  "value")
  
  com_results <- as.data.frame(do.call(cbind, lapply(ll2, as.data.frame)))
  rownames(com_results) <- genesID
  colnames(com_results) <- names(cluster_results)
  
  return(com_results)
}

##################################################
################# top results #####################
.top_results <- function(gene_results = NULL,
                         cluster_results, 
                         cluster = NULL,
                         select = "both",
                         high_low = NULL){
  stopifnot(
    is.list(cluster_results),
    #is.numeric(cluster_order), length(cluster_order) == 1L,
    is.character(select), length(select) == 1L
  )
  
  cluster_results <- lapply(cluster_results, as.data.frame)
  
  # If there is gene level results (stored in gene_results), return merged results (merge gene and cluster level results)
  if(!is.null(gene_results)){
    gene_results <- as.data.frame(gene_results)
    merge_results <- .merge_results(gene_results, cluster_results, select)
    merge_results <- as.data.frame(merge_results) %>% arrange(gene_FDR)
    return(res = merge_results)
  }
  
  if(!is.null(high_low)){
    stopifnot(is.character(high_low), length(high_low) == 1L)
    if( !(high_low %in% c("both", "BOTH", "high", "HIGH", "low", "LOW")) ){
      message("'high_low' should be one of: 'both', 'BOTH', 'high', 'HIGH', 'low', 'LOW'")
      return(NULL)
    }
  }
  if(!is.null(cluster)){
    stopifnot(is.character(cluster), length(cluster) == 1L)
  }
  if( !(select %in% c("p_adj", "FDR", "p_val", "logFC")) ){
    message("'select' should be one of: 'p_adj', 'FDR', 'p_val', 'logFC'")
    return(NULL)
  }
  if(select == 'p_adj'){select = "FDR"}
  if(select == 'p_val'){select = "PValue"}
  
  # Extract and combine results of the selected column ('select') from each sample
  com_results <- .getValueRes(cluster_results = cluster_results, select = select)
  # Top genes -> for each gene/row, min p-values/FDR; for each gene, return a vector of layer names
  #rank_results_list <- apply(com_results,1,function(x) which(x == min(x)))
  #rank_results <- as.data.frame(paste(lapply(rank_results_list,names),sep = ","))
  sel_layer_min <- apply(com_results,1,which.min)
  rank_results<-as.data.frame(colnames(com_results)[sel_layer_min])
  colnames(rank_results) <- "Layer"
  # For each gene, extract LR, PValue and FDR based on the layer name (i.e., the layer has min FDR)
  rank_results <- cbind(rank_results, 
                        lapply(1:nrow(com_results), function(x) cluster_results[[sel_layer_min[x]]][x,]) %>% 
                          bind_rows()) 
  # Rename columns; 3:7 -> "LR", "logCPM", "logFC", "PValue", "FDR"
  colnames(rank_results)[3:7] <- paste0("Layer_", colnames(rank_results)[3:7])
  colnames(com_results) <- paste0(names(cluster_results), "_", select)
  # Add the "gene id" column
  com_results['gene_id'] <- rownames(com_results)
  rank_results['gene_id'] <- rownames(com_results)
  rownames(rank_results) <- rank_results$gene_id
  # Add the results of identified layers (i.e., rank_results) into com_results
  com_results <- merge(rank_results, com_results, by = "gene_id") 
  # Check genes: p-values/FDR for all layers are same
  # gene_all_layer <- as.data.frame(com_results==apply(com_results,1, min))
  # gene_all_layer <- rownames(gene_all_layer %>% filter(rowSums(gene_all_layer) != 1))
  rownames(com_results) <- com_results$gene_id
  # Sort results based on FDR
  com_results <- as.data.frame(com_results) %>% arrange(Layer_FDR)
  # Sort results based on the specific cluster i; return top genes for cluster i
  if(!is.null(cluster)){
    gene_name <- com_results[order(com_results[, paste0((cluster), "_FDR")], decreasing = FALSE), 'gene_name']
    # Reorder rows of rank_results based on FDR of the specific cluster i
    top_results_cluster <- rank_results[match(gene_name, rank_results$gene_name),]
    # Subset Layer == cluster i
    top_genes_cluster <- subset(top_results_cluster, Layer == (cluster))
  }
  
  # Sort results based on the specific cluster i; return top genes for cluster i
  if(!is.null(cluster) & !is.null(high_low)){
    sel_FC <- .getValueRes(cluster_results = cluster_results, select = 'logFC')
    
    if(high_low %in% c("HIGH","high")){
      FC_high_genes <- rownames(sel_FC[sel_FC[,(cluster)] >= 0,])
      high_genes <- as.data.frame(top_genes_cluster[top_genes_cluste$gene_name %in% FC_high_genes,])
      low_genes <- NULL
    }
    if(high_low %in% c("LOW","low")){
      FC_low_genes <- rownames(sel_FC[sel_FC[,(cluster)] < 0,])
      low_genes <- as.data.frame(top_genes_cluster[top_genes_cluster$gene_name %in% FC_low_genes,])
      high_genes <- NULL
    }
    if(high_low %in% c("BOTH","both")){
      FC_high_genes <- rownames(sel_FC[sel_FC[,(cluster)] >= 0,])
      high_genes <- as.data.frame(top_genes_cluster[top_genes_cluster$gene_name %in% FC_high_genes,])
      
      FC_low_genes <- rownames(sel_FC[sel_FC[,(cluster)] < 0,])
      low_genes <- as.data.frame(top_genes_cluster[top_genes_cluster$gene_name %in% FC_low_genes,])
    }
    
  }
  if(is.null(high_low)){
    if(is.null(cluster)){
      return(res = com_results)
    }else{
      return(top_genes = top_genes_cluster)
    }
  }else{
    return(list(#res = com_results, 
      #top_genes = top_genes_cluster,
      high_genes = high_genes,
      low_genes = low_genes
    ))
  }
}
