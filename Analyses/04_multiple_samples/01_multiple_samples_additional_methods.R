## Get Enrichment statistics for snRNA-seq data
library("spatialLIBD")
library("SingleCellExperiment")
library(Seurat)
library(data.table)
registration_pseudobulk_modified <-
  function(
    sce,
    var_registration,
    var_sample_id,
    covars = NULL,
    min_ncells = 10,
    pseudobulk_rds_file = NULL) {
    ## Check that inputs are correct
    stopifnot(is(sce, "SingleCellExperiment"))
    stopifnot(var_registration %in% colnames(colData(sce)))
    stopifnot(var_sample_id %in% colnames(colData(sce)))
    stopifnot(all(
      !c("registration_sample_id", "registration_variable") %in% colnames(colData(sce))
    ))
    
    ## Avoid any incorrect inputs that are otherwise hard to detect
    stopifnot(!var_registration %in% covars)
    stopifnot(!var_sample_id %in% covars)
    stopifnot(var_registration != var_sample_id)
    
    ## Check that the values in the registration variable are ok
    uniq_var_regis <- unique(sce[[var_registration]])
    if (any(grepl("\\+|\\-", uniq_var_regis))) {
      stop(
        "Remove the + and - signs in colData(sce)[, '",
        var_registration,
        "'] to avoid downstream issues.",
        call. = FALSE
      )
    }
    
    ## Pseudo-bulk for our current BayesSpace cluster results
    message(Sys.time(), " make pseudobulk object")
    ## I think this needs counts assay
    sce_pseudo <- scuttle::aggregateAcrossCells(
      sce,
      DataFrame(
        registration_variable = sce[[var_registration]],
        registration_sample_id = sce[[var_sample_id]]
      )
    )
    colnames(sce_pseudo) <-
      paste0(
        sce_pseudo$registration_sample_id,
        "_",
        sce_pseudo$registration_variable
      )
    
    ## Check that the covariates are present
    if (!is.null(covars)) {
      for (covariate_i in covars) {
        if (sum(is.na(sce_pseudo[[covariate_i]])) == ncol(sce_pseudo)) {
          stop(
            "Covariate '",
            covariate_i,
            "' has all NAs after pseudo-bulking. Might be due to not being a sample-level covariate.",
            call. = FALSE
          )
        }
      }
    }
    
    ## Drop pseudo-bulked samples that had low initial contribution
    ## of raw-samples. That is, pseudo-bulked samples that are not
    ## benefiting from the pseudo-bulking process to obtain higher counts.
    if (!is.null(min_ncells)) {
      message(
        Sys.time(),
        " dropping ",
        sum(sce_pseudo$ncells < min_ncells),
        " pseudo-bulked samples that are below 'min_ncells'."
      )
      sce_pseudo <- sce_pseudo[, sce_pseudo$ncells >= min_ncells]
    }
    
    if (is.factor(sce_pseudo$registration_variable)) {
      ## Drop unused var_registration levels if we had to drop some due
      ## to min_ncells
      sce_pseudo$registration_variable <- droplevels(sce_pseudo$registration_variable)
    }
    
    ## Drop lowly-expressed genes
    # message(Sys.time(), " drop lowly expressed genes")
    # keep_expr <-
    #   edgeR::filterByExpr(sce_pseudo, group = sce_pseudo$registration_variable)
    # sce_pseudo <- sce_pseudo[which(keep_expr), ]
    # 
    ## Compute the logcounts
    message(Sys.time(), " normalize expression")
    logcounts(sce_pseudo) <-
      edgeR::cpm(edgeR::calcNormFactors(sce_pseudo),
                 log = TRUE,
                 prior.count = 1
      )
    
    if (is(sce_pseudo, "SpatialExperiment")) {
      ## Drop things we don't need
      spatialCoords(sce_pseudo) <- NULL
      imgData(sce_pseudo) <- NULL
    }
    if (!is.null(pseudobulk_rds_file)) {
      message(Sys.time(), " saving sce_pseudo to ", pseudobulk_rds_file)
      saveRDS(sce_pseudo, file = pseudobulk_rds_file)
    }
    
    ## Done!
    return(sce_pseudo)
  }
run_LIBD_registration1 <- function(sce.combined){
  ## Add gene-level information
  rowData(sce.combined)$gene_id <- rownames(sce.combined)
  rowData(sce.combined)$gene_name <- rownames(sce.combined)
  sce.combined$Cluster <- as.character(sce.combined$JointCluster)
  sce.combined$Cluster[sce.combined$JointCluster == "in"] <- "up"
  sce.combined$Cluster[sce.combined$JointCluster == "out"] <- "down"
  ## Perform the spatial registration
  sce_pseudo <- registration_pseudobulk_modified(
    sce = sce.combined,
    var_registration = "Cluster", 
    var_sample_id = "sample_name", 
    min_ncells = 0)
  
  registration_mod <- registration_model(sce_pseudo)
  head(registration_mod)
  block_cor <- registration_block_cor(sce_pseudo, registration_mod)
  
  ## compute ANOVA statistics
  # This function computes the gene ANOVA F-statistics 
  # (at least one group is different from the rest). 
  # These F-statistics can be used for spatial registration with layer_stat_cor() and related functions. 
  # Although, they are more typically used for identifying ANOVA-marker genes.
  # results_anova_nan <- registration_stats_anova(sce_pseudo,
  #                                               block_cor = NaN,
  #                                               gene_name = "gene_name",
  #                                               gene_ensembl = "gene_name"
  # )
  # # 
  # results_anova <- registration_stats_anova(sce_pseudo,
  #                                           block_cor,
  #                                           gene_name = "gene_name",
  #                                           gene_ensembl = "gene_name"
  # )
  # 
  # results_anova_merged <- merge(results_anova_nan, results_anova_nocovar)
  # Error: You need 'var_registration' to have at least 3 different values to compute an F-statistic.
  
  
  ## compute enrichment statistics
  ## This function computes the gene enrichment t-statistics (one group > the rest). 
  ## These t-statistics are the ones typically used for spatial registration with 
  ## layer_stat_cor() and related functions.
  results_enrichment <- registration_stats_enrichment(sce_pseudo,
                                                      block_cor,
                                                      gene_name = "gene_name",
                                                      gene_ensembl = "gene_name"
  )
  # Specifying `block_cor = NaN` then ignores the correlation structure
  results_enrichment_nan <- registration_stats_enrichment(sce_pseudo,
                                                          block_cor = NaN,
                                                          gene_name = "gene_name",
                                                          gene_ensembl = "gene_name"
                                                          
  )
  
  ## compute pairwise statistics
  ## This function computes the gene pairwise t-statistics (one group > another, for all combinations). 
  ## These t-statistics can be used for spatial registration with layer_stat_cor() and related functions. 
  ## Although, they are more typically used for identifying pairwise-marker genes.
  
  results_pairwise <- registration_stats_pairwise(sce_pseudo,
                                                  registration_mod,
                                                  block_cor,
                                                  gene_name = "gene_name",
                                                  gene_ensembl = "gene_name"
  )
  
  results_pairwise_nan <- registration_stats_pairwise(sce_pseudo,
                                                      registration_mod,
                                                      block_cor = NaN,
                                                      gene_name = "gene_name",
                                                      gene_ensembl = "gene_name"
  )
  
  return(list(    
    results_enrichment, results_enrichment_nan,
    results_pairwise, results_pairwise_nan))
}

run_LIBD_registration <- function(sce.combined){
  ## Add gene-level information
  rowData(sce.combined)$gene_id <- rownames(sce.combined)
  rowData(sce.combined)$gene_name <- rownames(sce.combined)
  sce.combined$Cluster <- as.character(sce.combined$JointCluster)
  # sce.combined$Cluster[sce.combined$JointCluster == "in"] <- "up"
  # sce.combined$Cluster[sce.combined$JointCluster == "out"] <- "down"
  ## Perform the spatial registration
  sce_pseudo <- registration_pseudobulk_modified(
    sce = sce.combined,
    var_registration = "Cluster", 
    var_sample_id = "sample_name", 
    min_ncells = 0)
  
  registration_mod <- registration_model(sce_pseudo)
  head(registration_mod)
  block_cor <- registration_block_cor(sce_pseudo, registration_mod)
  
  ## compute ANOVA statistics
  # This function computes the gene ANOVA F-statistics 
  # (at least one group is different from the rest). 
  # These F-statistics can be used for spatial registration with layer_stat_cor() and related functions. 
  # Although, they are more typically used for identifying ANOVA-marker genes.
  results_anova_nan <- registration_stats_anova(sce_pseudo,
                                                block_cor = NaN,
                                                gene_name = "gene_name",
                                                gene_ensembl = "gene_name"
  )
  # 
  results_anova <- registration_stats_anova(sce_pseudo,
                                            block_cor,
                                            gene_name = "gene_name",
                                            gene_ensembl = "gene_name"
  )
  # 
  # results_anova_merged <- merge(results_anova_nan, results_anova_nocovar)
  # Error: You need 'var_registration' to have at least 3 different values to compute an F-statistic.
  
  
  ## compute enrichment statistics
  ## This function computes the gene enrichment t-statistics (one group > the rest). 
  ## These t-statistics are the ones typically used for spatial registration with 
  ## layer_stat_cor() and related functions.
  results_enrichment <- registration_stats_enrichment(sce_pseudo,
                                                      block_cor,
                                                      gene_name = "gene_name",
                                                      gene_ensembl = "gene_name"
  )
  # Specifying `block_cor = NaN` then ignores the correlation structure
  results_enrichment_nan <- registration_stats_enrichment(sce_pseudo,
                                                          block_cor = NaN,
                                                          gene_name = "gene_name",
                                                          gene_ensembl = "gene_name"
                                                          
  )
  
  ## compute pairwise statistics
  ## This function computes the gene pairwise t-statistics (one group > another, for all combinations). 
  ## These t-statistics can be used for spatial registration with layer_stat_cor() and related functions. 
  ## Although, they are more typically used for identifying pairwise-marker genes.
  
  results_pairwise <- registration_stats_pairwise(sce_pseudo,
                                                  registration_mod,
                                                  block_cor,
                                                  gene_name = "gene_name",
                                                  gene_ensembl = "gene_name"
  )
  
  results_pairwise_nan <- registration_stats_pairwise(sce_pseudo,
                                                      registration_mod,
                                                      block_cor = NaN,
                                                      gene_name = "gene_name",
                                                      gene_ensembl = "gene_name"
  )
  
  return(list(    results_anova, results_anova_nan,
                  results_enrichment, results_enrichment_nan,
                  results_pairwise, results_pairwise_nan))
}

run_scran_findMarkers <- function(sce.combined){
  sce.combined$Cluster <- as.character(sce.combined$JointCluster)
  # sce.combined$Cluster[sce.combined$JointCluster == "in"] <- "up"
  # sce.combined$Cluster[sce.combined$JointCluster == "out"] <- "down"
  layer = as.factor(colData(sce.combined)[["Cluster"]])
  markers_block <- scran::findMarkers(
    sce.combined, groups = layer,
    block = sce.combined$sample_name, 
    pval.type = "all"
  )
  return(markers_block)
}

run_seurat_FindAllMarkers <- function(sce.combined){
  sce.combined$Cluster <- as.character(sce.combined$JointCluster)
  # sce.combined$Cluster[sce.combined$JointCluster == "in"] <- "up"
  # sce.combined$Cluster[sce.combined$JointCluster == "out"] <- "down"
  counts <- assays(sce.combined)[[1]]
  colnames(counts) <- paste0(colData(sce.combined)$cell_ID, "_", colData(sce.combined)$sample_name)
  seurat <- CreateSeuratObject(counts = counts)
  t2 <- seurat@meta.data
  t <- data.frame(colData(sce.combined))
  t3 <- data.frame(t,t2)
  seurat@meta.data <- t3
  layer = as.factor(colData(sce.combined)[["Cluster"]])
  Idents(seurat) <- layer
  # find markers for every cluster compared to all remaining cells
  seuratResults1 <- Seurat::FindAllMarkers(seurat,
                                           grouping.var = "sample_name",
                                           logfc.threshold = -Inf,
                                           min.pct = 0,
                                           test.use = "wilcox",
                                           min.diff.pt = -Inf,
                                           min.cells.feature = 0,
                                           min.cells.group = 0,
                                           return.thresh = 2#, min.pct = 0.25, logfc.threshold = 0.25
  )
  seuratResults2 <- Seurat::FindAllMarkers(seurat, 
                                           recorrect_umi = FALSE, test.use="MAST", 
                                           latent.vars="sample_name",
                                           logfc.threshold = -Inf,
                                           min.diff.pt = -Inf,
                                           min.cells.feature = 0,
                                           min.cells.group = 0,
                                           return.thresh = 2,
                                           min.pct = 0)
  return(list(seuratResults1, seuratResults2))
}
# seuratResults3 <- FindConservedMarkers(seurat, 
#                                        ident.1 = "down", ident.2 = NULL, 
#                                        grouping.var="sample_name",
#                      verbose = TRUE)
## only return a ranked list of putative conserved markers

# https://github.com/satijalab/seurat/issues/3970
# Running FindAllMarkers on the integrated data is not recommended. I would run it within each of your datasets and then compare or run FindConservedMarkers and have you biological replicates marked as variable for iding the two groups. 
# I would use this for reference https://www.biostars.org/p/409790/.

# https://github.com/satijalab/seurat/issues/6302
# since MAST DE test can model donor_ID, it will be able to mitigate such issue. 
# The MAST test will identify the DEGs that are consistently elevated/down-regulated across all 
# the samples in your case group. As for sample-specific "DEGs", 
# they will be captured by the Donor_ID term in the MAST model. 
# So personally, MAST model might provide a more stringent list of DEGs than Wilcox test.