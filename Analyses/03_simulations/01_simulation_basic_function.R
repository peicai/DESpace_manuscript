source("./Analysis/03_simulations/01_methods_in_R.R")

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

#' @title addCellMetadata
#' @description adds cell metadata to the giotto object
#' @param gobject giotto object
#' @param new_metadata new cell metadata to use (data.table, data.frame, ...)
#' @param vector_name (optional) custom name if you provide a single vector
#' @param by_column merge metadata based on cell_ID column in pDataDT (default = FALSE)
#' @param column_cell_ID column name of new metadata to use if by_column = TRUE
#' @return giotto object
#' @details You can add additional cell metadata in two manners:
#' \itemize{
#'   \item{1. Provide a data.table or data.frame with cell annotations in the same order as the cell_ID column in pDataDT(gobject) }
#'   \item{2. Provide a data.table or data.frame with cell annotations and specificy which column contains the cell IDs, these cell IDs need to match with the cell_ID column in pDataDT(gobject)}
#' }
#' @export
addMetadata <- function(sce_object,
                        new_metadata,
                        vector_name = NULL,
                        by_column = FALSE,
                        column_cell_ID = NULL) {
  # data.table variables
  cell_ID = NULL
  
  cell_metadata = colData(sce_object)
  ordered_cell_IDs = rownames(colData(sce_object))
  
  if(is.vector(new_metadata) | is.factor(new_metadata)) {
    original_name = deparse(substitute(new_metadata))
    new_metadata = data.table::as.data.table(new_metadata)
    
    if(!is.null(vector_name) & is.character(vector_name)) {
      colnames(new_metadata) = vector_name
    } else {
      colnames(new_metadata) = original_name
    }
    
  } else {
    new_metadata = data.table::as.data.table(new_metadata)
  }
  
  if(is.null(column_cell_ID)) {
    column_cell_ID = 'cell_ID'
  }
  
  # overwrite columns with same name
  new_col_names = colnames(new_metadata)
  new_col_names = new_col_names[new_col_names != column_cell_ID]
  old_col_names = colnames(cell_metadata)
  old_col_names = old_col_names[old_col_names != 'cell_ID']
  same_col_names = new_col_names[new_col_names %in% old_col_names]
  
  
  if(length(same_col_names) >= 1) {
    cat('\n these column names were already used: ', same_col_names, '\n',
        'and will be overwritten \n')
    cell_metadata = as.data.table(cell_metadata)
    cell_metadata[, (same_col_names) := NULL]
  }
  
  
  
  if(by_column == FALSE) {
    cell_metadata = cbind(cell_metadata, new_metadata)
  } else {
    if(is.null(column_cell_ID)) stop('You need to provide cell_ID column')
    cell_metadata <- data.table::merge.data.table(as.data.table(cell_metadata), by.x = 'cell_ID',
                                                  new_metadata, by.y = column_cell_ID,
                                                  all.x = T)
  }
  
  # reorder
  cell_metadata = cell_metadata[match(ordered_cell_IDs, cell_ID)]
  
  rownames(cell_metadata) <- rownames(colData(sce_object))
  sce_object = SingleCellExperiment(assays=list(counts=counts(sce_object)), 
                                    colData=cell_metadata)
  return(sce_object)
}


## Simulate single-gene spatial patterns ####
## --------------------------------------- ##
#' @title simulateOneGenePatternGiottoObject
#' @name simulateOneGenePatternGiottoObject
#' @description Create a simulated spatial pattern for one selected gnee
#' @param gobject giotto object
#' @param pattern_name name of spatial pattern
#' @param pattern_cell_ids cell ids that make up the spatial pattern
#' @param gene_name selected gene
#' @param spatial_prob probability for a high expressing gene value to be part of the spatial pattern
#' @param gradient_direction direction of gradient
#' @param show_pattern show the discrete spatial pattern
#' @param pattern_colors 2 color vector for the spatial pattern
#' @param \dots additional parameters for (re-)normalizing
#' @return Reprocessed Giotto object for which one gene has a forced spatial pattern
#' @export
simulateOneGenePatternObject = function(sce_object,
                                        coordinate_name = c("row", "col"),
                                        pattern_name = 'pattern',
                                        cell_meta,
                                        #pattern_cell_ids = pattern_cell_ids,
                                        gene_name = NULL,
                                        spatial_prob = 0.95,
                                        gradient_direction = NULL,
                                        show_pattern = TRUE,
                                        pattern_colors = c('in' = 'green', 'out' = 'red'),
                                        layer_names = layer_names,
                                        simobject,
                                        save_raw = TRUE,
                                        save_dir = '~',
                                        save_name = 'counts',
                                        prob_i = 0.5,
                                        ...) {
  # data.table variables
  cell_ID = sdimx_y = sdimx = sdimy = NULL
  
  #newgobject = addCellMetadata(sce_object,
  #                             new_metadata = cell_meta[,c('cell_ID', pattern_name), with = F],
  #                             by_column = T,
  #                             column_cell_ID = 'cell_ID')
  newobject = addMetadata(sce_object,
                          new_metadata = cell_meta[,c('cell_ID', pattern_name), with = F],
                          by_column = T,
                          column_cell_ID = 'cell_ID')  
  
  ## merge cell metadata and cell coordinate data
  cell_meta = data.table::as.data.table(colData(newobject))
  cell_coord = data.table::as.data.table(cbind.data.frame(x=as.numeric(as.character(colData(sce_object)[[coordinate_name[1]]])),
                                                          y=as.numeric(as.character(colData(sce_object)[[coordinate_name[2]]])),
                                                          cell_ID = colData(sce_object)[['cell_ID']]))
  cell_meta = data.table::merge.data.table(cell_meta, cell_coord, by = 'cell_ID')
  
  ## get number of cells within pattern
  cell_number = nrow(cell_meta[get(pattern_name) == 'in'])
  
  
  ## normalized expression
  counts <- assay(newobject, "counts")
  libsizes <- colSums(counts)
  size.factors <- libsizes/mean(libsizes)
  logcounts(newobject) <- log2(t(t(counts)/size.factors) + 1)
  expr_data = assay(newobject, "logcounts")
  result_list = list()
  
  ## raw expression
  raw_expr_data = assays(newobject)$counts
  raw_result_list = list()
  
  
  ## create the spatial expression pattern for the specified gene
  # 1. rank all gene values from the cells from high to low
  # 2. move the highest expressing values to the spatial pattern using a probability
  #     - 0.5 is the control = random
  #     - 1 is perfection: all the highest values go to the pattern
  #     - 0.5 to 1 is decreasing noise levels
  
  if(is.null(gene_name)) stop('a gene name needs to be provided')
  
  # rank genes
  gene_vector = expr_data[rownames(expr_data) == gene_name, ]
  sort_expr_gene = sort(gene_vector, decreasing = T)
  
  # number of cells in and out the pattern
  total_cell_number = length(sort_expr_gene)
  remaining_cell_number = total_cell_number - cell_number
  print(paste0("total_cell_number: ", total_cell_number))
  print(paste0("cell_number: ", cell_number))
  print(paste0("remaining cell number: ", remaining_cell_number))
  # calculate outside probability
  outside_prob = 1 - spatial_prob
  prob_vector = c(rep(spatial_prob, cell_number), rep(outside_prob, remaining_cell_number))
  
  # first get the 'in' pattern sample values randomly
  sample_values = sample(sort_expr_gene, replace = F, size = cell_number, prob = prob_vector)
  
  # then take the remaining 'out' pattern values randomly
  remain_values = sort_expr_gene[!names(sort_expr_gene) %in% names(sample_values)]
  remain_values = sample(remain_values, size = length(remain_values))
  
  
  
  
  ## A. within pattern ##
  # ------------------- #
  in_cell_meta = cell_meta[get(pattern_name) == 'in']
  
  # if gradient is wanted
  # does not work with 0.5!! is not random!!
  if(!is.null(gradient_direction)) {
    # sort in_ids according to x, y or  xy coordinates to create gradient
    in_cell_meta[, sdimx_y := abs(sdimx)+ abs(sdimy)]
    # order according to gradient direction
    in_cell_meta = in_cell_meta[order(get(gradient_direction))]
  }
  in_ids = in_cell_meta$cell_ID
  
  # preparation for raw matrix
  sample_values_id_vector = names(sample_values)
  names(sample_values_id_vector) = in_ids
  
  
  ## B. outside pattern ##
  # -------------------- #
  out_ids = cell_meta[get(pattern_name) == 'out']$cell_ID
  
  # preparation for raw matrix
  remain_values_id_vector = names(remain_values)
  names(remain_values_id_vector) = out_ids
  
  
  
  
  ## raw matrix
  # swap the cell ids #
  raw_gene_vector = raw_expr_data[rownames(raw_expr_data) == gene_name,]
  
  raw_new_sample_vector = raw_gene_vector[sample_values_id_vector]
  names(raw_new_sample_vector) = names(sample_values_id_vector)
  
  raw_new_remain_vector = raw_gene_vector[remain_values_id_vector]
  names(raw_new_remain_vector) = names(remain_values_id_vector)
  
  new_sim_raw_values = c(raw_new_sample_vector, raw_new_remain_vector)
  new_sim_raw_values = new_sim_raw_values[names(raw_gene_vector)]
  
  # change the original matrices
  raw_expr_data[rownames(raw_expr_data) == gene_name,] = new_sim_raw_values
  assays(newobject)$counts = raw_expr_data
  
  # add counts for selected gene to simulate gobject
  new_raw_expr_data = assays(simobject)$counts
  new_raw_expr_data[rownames(new_raw_expr_data) == gene_name,] = new_sim_raw_values
  assays(simobject)$counts = new_raw_expr_data
  
  # recalculate normalized values
  #simgobject <- normalizeGiotto(gobject = simgobject, ...)
  #simgobject <- addStatistics(gobject = simgobject)
  simobject <- scuttle::logNormCounts(simobject)
  return(simobject)
  
}



#' @title CombineSimulatedGiottoObject
#' @name CombineSimulatedGiottoObject
#' @description Creates final giotto object (SV_percent \* spaital_probs[1] + (1-SV_percent)\*spatial_probs[2])
#' @keywords internal
CombineSimulatedObject = function(sce_object,
                                  pattern_name = 'pattern',
                                  coordinate_name = c("row", "col"),
                                  pattern_cell_ids = NULL,
                                  gene_name = NULL,
                                  spatial_probs = c(0.5, 1),
                                  rep_i = 1,
                                  SV_percent = 0.33,
                                  layer_names=layer_names,
                                  save_plot = FALSE,
                                  save_raw = TRUE,
                                  save_norm = TRUE,
                                  save_dir = '~',
                                  verbose = TRUE,
                                  run_simulations = TRUE,
                                  simobject,
                                  save_name = 'probs',
                                  default = T,
                                  ...) {
  set.seed(123)
  if(is.null(pattern_cell_ids)) {
    stop('pattern_cell_ids can not be NULL \n')
  }
  prob_list = list()
  cell_meta_original = data.table::as.data.table(colData(sce_object))
  
  if(length(spatial_probs) == 2){
    final_object = sce_object
    GroudTruth = data.frame(matrix(ncol = 3, nrow = dim(sce_object)[2]))
    colnames(GroudTruth) = c("high_SVG", "low_SVG", "non_SVG")
    ## SV genes
    sample_genes = sample(as.vector(rownames(assays(sce_object)$counts)), round(length(gene_name)*SV_percent))
    
    ## For 3 “easy” patterns, simulate 50% of SV patterns with in = highly abundant area, and 50% of SV patterns with in = lowly abundant area.
    high_SVG = sample(as.vector(sample_genes), round(length(sample_genes)*1/2))
    low_SVG = sample_genes[sample_genes %notin% high_SVG]
    non_SVG = rownames(assays(sce_object)$counts)[rownames(assays(sce_object)$counts) %notin% sample_genes]
    
    ## Save SVGs
    subdir = paste0(save_dir,'/',pattern_name,'/')
    resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
    
    status <- ifelse(sample_genes %in% high_SVG, "high", "low")
    SVG <- as.data.frame(cbind(gene_ids = sample_genes, status))
    
    write.table(SVG,
                file = paste0(resultdir,  'SV_genes.txt'),
                sep = '\t')
    
    gene_split <- list(non_SVG, sample_genes)
    GroudTruth$high_SVG <- with(cell_meta_original,ifelse(cell_ID %in% pattern_cell_ids, 'in', 'out'))
    GroudTruth$non_SVG <- with(cell_meta_original,ifelse(cell_ID %in% pattern_cell_ids, 'in', 'out'))
    GroudTruth$low_SVG <- with(cell_meta_original,ifelse(cell_ID %in% pattern_cell_ids, 'out', 'in'))
    write.table(as.data.frame(GroudTruth),
                file = paste0(resultdir, 'Ground_Truth.csv'),
                sep = '\t')
    # for first spatial probability
    for(p in c(1:2)){
      prob_i = spatial_probs[p]
      gene_list = gene_split[[p]]
      
      if(verbose) cat('\n \n start with ', prob_i, '\n \n')
      if(verbose) cat('\n \n repetition = ', rep_i, '\n \n')
      
      for(gene_ind in 1:length(gene_list)) {
        gene = gene_list[gene_ind]
        print(gene)
        print(paste0(round(gene_ind/length(gene_list)*100),"%"))
        plot_name = paste0('plot_',gene,'_prob', prob_i, '_rep', rep_i)
        cell_meta_original = data.table::as.data.table(colData(sce_object))
        
        ## create and add annotation for pattern
        if(gene %in% c(non_SVG, high_SVG)){
          cell_meta_original[, (pattern_name) := ifelse(cell_ID %in% pattern_cell_ids, 'in', 'out')]
        }else if(gene %in% low_SVG){
          cell_meta_original[, (pattern_name) := ifelse(cell_ID %in% pattern_cell_ids, 'out', 'in')]
        }
        
        gg = simulateOneGenePatternObject(sce_object,
                                          coordinate_name = coordinate_name,
                                          pattern_name = pattern_name,
                                          cell_meta = cell_meta_original,
                                          #pattern_cell_ids = pattern_cell_ids,
                                          gene_name = gene,
                                          spatial_prob = prob_i,
                                          gradient_direction = NULL,
                                          show_pattern = FALSE,
                                          pattern_colors = c('in' = 'green', 'out' = 'red'),
                                          layer_names = layer_names,
                                          simobject,
                                          prob_i = prob_i,
                                          save_dir = save_dir
        )
        simobject = gg
      }
      # save raw data
      if(save_raw == TRUE) {
        
        folder_path = paste0(save_dir, '/', pattern_name)
        if(!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
        subdir = paste0(save_dir,'/',pattern_name,'/')
        resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
        write.table(x = as.matrix(counts(gg)),
                    file = paste0(resultdir,  save_name,'_',prob_i, '_raw_data.txt'),
                    sep = '\t')
        save(gg, file = paste0(resultdir, save_name,'_',prob_i,'_final_object.rda'))
        
      }
      # save normalized data
      if(save_norm == TRUE) {
        
        folder_path = paste0(save_dir, '/', pattern_name)
        if(!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
        subdir = paste0(save_dir,'/',pattern_name,'/')
        resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
        write.table(x = as.matrix(assays(gg)$logcounts),
                    file = paste0(resultdir, save_name,'_',prob_i, '_norm_data.txt'),
                    sep = '\t')
      }
      
      raw_expr_data = assays(gg)$counts
      new_sim_raw_values = raw_expr_data[rownames(raw_expr_data) %in% gene_list,]
      raw_expr_data_2 = assays(final_object)$counts
      raw_expr_data_2[rownames(raw_expr_data_2) %in% gene_list,] = new_sim_raw_values
      assays(final_object)$counts = raw_expr_data_2
    }
    
    write.table(x = as.matrix(assays(final_object)$counts),
                file = paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_raw_data.txt'),
                sep = '\t')
    
    save(final_object, file = paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))
    return(final_object)
  }else {
    genes_split = split(rownames(counts(sce_object)),sample(1:length(SV_percent),length(rownames(counts(sce_object))),replace=TRUE,prob=SV_percent))
    combined_object = sce_object
    GroudTruth = data.frame(matrix(ncol = length(spatial_probs), nrow = dim(sce_object)[2]))
    colnames(GroudTruth) = paste0(spatial_probs,"_",seq(spatial_probs))
    for(i in c(1:length(spatial_probs))){
      prob_i = spatial_probs[i]
      gene_list = genes_split[[i]]
      if(verbose) cat('\n \n start with ', prob_i, '_',i,'\n \n')
      simobject = sce_object
      
      if(spatial_probs[1] == 0.5){
        ifelse(i==1, cell_ids <-  colnames(assays(sce_object)$counts), cell_ids <-  pattern_cell_ids[[i-1]])
      }else { #ifelse(i==1, cell_ids <-  colnames(assays(sce_object)$counts), 
        cell_ids <- pattern_cell_ids[[i]]}
      
      cell_meta_original[, (pattern_name) := ifelse(cell_ID %in% cell_ids, 'in', 'out')]
      
      #GroudTruth[,i] <- with(cell_meta_original,ifelse(cell_ID %in% pattern_cell_ids, 'in', 'out'))
      GroudTruth[,i] <- ifelse(cell_meta_original$cell_ID %in% cell_ids, 'in', 'out')
      
      for(gene_ind in 1:length(gene_list)) {
        gene = gene_list[gene_ind]
        print(gene)
        print(paste0(round(gene_ind/length(gene_list)*100),"%"))
        plot_name = paste0('plot_',gene,'_prob', prob_i, '_rep', rep_i)
        
        gg <-   simulateOneGenePatternObject(sce_object,
                                             coordinate_name = coordinate_name,
                                             pattern_name = pattern_name,
                                             cell_meta = cell_meta_original,
                                             #pattern_cell_ids = cell_ids,
                                             gene_name = gene,
                                             spatial_prob = prob_i,
                                             gradient_direction = NULL,
                                             show_pattern = FALSE,
                                             pattern_colors = c('in' = 'green', 'out' = 'red'),
                                             layer_names = layer_names,
                                             simobject,
                                             prob_i = prob_i,
                                             save_dir = save_dir
        )
        #simgobject <- eval(parse(text = paste0("final_gobject_",i)))
        simobject <- gg
        
      }
      
      folder_path = paste0(save_dir, '/', pattern_name)
      if(!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
      subdir = paste0(save_dir,'/',pattern_name,'/')
      resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
      write.table(x = as.matrix(assays(gg)$counts),
                  file = paste0(resultdir,  save_name,'_',prob_i, '_',i,'_raw_data.txt'),
                  sep = '\t')
      save(gg, file = paste0(resultdir, save_name,'_',prob_i,'_',i,'_final_object.rda'))
      
      raw_expr_data = assays(gg)$counts
      new_sim_raw_values = raw_expr_data[rownames(raw_expr_data) %in% gene_list,]
      raw_expr_data_2 = assays(combined_object)$counts
      raw_expr_data_2[rownames(raw_expr_data_2) %in% gene_list,] = new_sim_raw_values
      assays(combined_object)$counts = raw_expr_data_2
      
      
    }
    write.table(as.data.frame(GroudTruth),
                file = paste0(resultdir, 'Ground_Truth.csv'),
                sep = '\t')   
    ## combine gobjects of diff spatial_probs/patterns
    subdir = paste0(save_dir,'/',pattern_name,'/')
    resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
    names(genes_split) <- paste0(rep(c("pattern_"),length(spatial_probs)), 
                                 rep(1:length(spatial_probs)))
    save(genes_split, file = paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_all_genes.rda'))
    write.table(x = as.matrix(unlist(genes_split[-1])),
                file = paste0(resultdir, save_name,'_',spatial_probs[2], '_selected_genes.txt'),
                sep = '\t')
    write.table(x = as.matrix(unlist(genes_split[-1])),
                file = paste0(resultdir,  'SV_genes.txt'),
                sep = '\t')
    write.table(x = as.matrix(assays(combined_object)$counts),
                file = paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_raw_data.txt'),
                sep = '\t')
    
    save(combined_object, file = paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))
    return(combined_object)
  }
  
}


run_spatial_sim_tests_methods = function(sce_object = final_object,raw_sce = raw_object,rep_i = rep_i,
                                         pattern_name = 'pattern',
                                         coordinate_name = coordinate_name,
                                         #pattern_cell_ids = NULL,
                                         gene_name = NULL,
                                         spatial_prob = 0.95,
                                         show_pattern = FALSE,
                                         spat_methods = c('SV_edgeR_counts',
                                                          'SV_edgeR_CPMs','meringue','spatialDE','spark'
                                         ),
                                         spat_methods_params = list(NA,NA, NA,NA,NA),
                                         spat_methods_names = c('SV_edgeR_counts',
                                                                'SV_edgeR_CPMs','meringue','spatialDE','spark'
                                         ),
                                         save_plot = FALSE,
                                         save_raw = FALSE,
                                         save_norm = FALSE,
                                         save_dir = '~',
                                         dir = '~',
                                         save_name = 'plot',
                                         run_simulations = TRUE,
                                         layer_names=layer_names,
                                         original_layer_names = "spatial.cluster",
                                         default = T,
                                         cluster_method = 'BayesSpace',
                                         platform = platform,
                                         ...){
  genes = prob = time = adj.p.value = method = p.val = sd = qval = pval = g = adjusted_pvalue = NULL
  
  ## do simulations ##
  if(run_simulations == TRUE) {
    
    result_list = list()
    for(test in 1:length(spat_methods)) {
      
      # method
      selected_method = spat_methods[test]
      if(!selected_method %in% c(#'spatialDE,
        'nnSVG2','spark_x2',
        'spark', 'spark_x',#'SV_DESeq2,
        'SV_edgeR_counts',#'SV_edgeR_CPMs,
        'meringue', 'nnSVG')) {
        stop(selected_method, ' is not a know spatial method \n')
      }
      print(selected_method)
      # params
      selected_params = spat_methods_params[[test]]
      
      if(length(selected_params) == 1) {
        
        if(is.na(selected_params)) {
          
          
          
          if(selected_method == 'spark') {
            selected_params = list(sce_object = sce_object,
                                   percentage = 0,
                                   min_count = 0,
                                   coordinate_name = coordinate_name,
                                   num_core = 5,
                                   covariates = NULL,
                                   return_object = c('data.table')
            )
            
          } else if(selected_method == 'spark_x') {
            selected_params = list(sce_object = sce_object,
                                   coordinate_name = coordinate_name,
                                   num_core = 5,
                                   return_object = c('data.table'))
          } else if(selected_method == 'spark_x2') {
            selected_params = list(sce_object = sce_object,
                                   coordinate_name = coordinate_name,
                                   num_core = 3,
                                   layer_names = layer_names,
                                   dir = dir,
                                   return_object = c('data.table'))
            
            
            
          } else if(selected_method %in% c('SV_edgeR_counts')) {
            selected_params = list(sce_object = sce_object,
                                   min_count = 0,
                                   coordinate_name = coordinate_name,
                                   num_core = 4,
                                   prior.count = 10, # edgeR:cpm
                                   covariates = NULL,
                                   return_object = c('data.table'),
                                   layer_names = layer_names,
                                   original_layer_names = original_layer_names,
                                   default = default,
                                   save = F,
                                   dir = dir,
                                   pattern_name = pattern_name,
                                   cluster_method = cluster_method,
                                   platform = platform)
            
            
          } else if(selected_method == 'meringue') {
            selected_params = list(sce_object = sce_object,
                                   raw_sce = raw_sce,
                                   min.reads = 0, 
                                   min.lib.size = 0,
                                   coordinate_name= coordinate_name,
                                   covariates = NULL,
                                   return_object = c('data.table')
            )
            
            
          } else if(selected_method == 'nnSVG') {
            selected_params = list(sce_object = sce_object,
                                   coordinate_name = coordinate_name,
                                   num_core = 4,
                                   return_object = c('data.table')
            )
            
          } else if(selected_method == 'nnSVG2') {
            selected_params = list(sce_object = sce_object,
                                   coordinate_name = coordinate_name,
                                   num_core = 3,
                                   layer_names = layer_names,
                                   dir = dir,
                                   return_object = c('data.table'))
          }
          # end if selected_method
        } # end if is.na(selected_params)
        
      } # end if length(selected_params) == 1
      
      # name
      selected_name = spat_methods_names[test]
      
      
      ## RUN Spatial Analysis ##
      message("Run spatial analysis")
      if(selected_method == 'SV_edgeR_counts') 
      {
        message('Method: SV_edgeR_counts')
        ## spark
        start = proc.time()
        SV_edgeR_spatialgenes_sim = do.call('SV_edgeR_counts', c(
          selected_params))
        #print(class(SV_edgeR_spatialgenes_sim))
        #print(str(SV_edgeR_spatialgenes_sim))
        #print(gene_name)
        SV_edgeR_result = SV_edgeR_spatialgenes_sim#[genes == gene_name]
        SV_edgeR_time = proc.time() - start
        
        SV_edgeR_result[, pres := rep_i]
        SV_edgeR_result[, time := SV_edgeR_time[['elapsed']] ]
        
        spatial_gene_results = SV_edgeR_result[,.(genes, FDR,PValue, pres, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value','p.value', 'repetition', 'time')
        spatial_gene_results[, method := 'SV_edgeR_counts'] 
        #subdir = paste0(save_dir,'/',pattern_name,'/')
        if(cluster_method == 'BayesSpace'){
          save(spatial_gene_results, file = paste0(dir,"result_SV_edgeR_counts.rda"))
        }else {save(spatial_gene_results, file = paste0(dir,"result_SV_edgeR_",cluster_method,"counts.rda"))}
        
        
      } else if(selected_method == 'spark') 
      {
        message('Method: spark')
        ## spark
        start = proc.time()
        print(sce_object)
        spark_spatialgenes_sim = do.call('spark', c(
          selected_params))
        #print(class(spark_spatialgenes_sim))
        #print(str(spark_spatialgenes_sim))
        #print(gene_name)
        spark_result = spark_spatialgenes_sim#[genes == gene_name]
        spark_time = proc.time() - start
        
        spark_result[, pres := rep_i]
        spark_result[, time := spark_time[['elapsed']] ]
        
        spatial_gene_results = spark_result[,.(genes, adjusted_pvalue,combined_pvalue, pres, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value','p.value', 'repetition', 'time')
        spatial_gene_results[, method := 'spark']
        
        #subdir = paste0(save_dir,'/',pattern_name,'/')
        save(spatial_gene_results, file = paste0(dir,"result_SV_spark.rda"))
        
        
      } else if(selected_method == 'spark_x') 
      {
        message('Method: spark_x')
        ## spark
        start = proc.time()
        sparkx_spatialgenes_sim = do.call('spark_x', c(selected_params))
        #print(class(spark_spatialgenes_sim))
        #print(str(spark_spatialgenes_sim))
        #print(gene_name)
        sparkx_result = sparkx_spatialgenes_sim#[genes == gene_name]
        sparkx_time = proc.time() - start
        
        sparkx_result[, pres := rep_i]
        sparkx_result[, time := sparkx_time[['elapsed']] ]
        
        spatial_gene_results = sparkx_result[,.(genes, adjustedPval,combinedPval, pres, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value','p.value', 'repetition', 'time')
        spatial_gene_results[, method := 'spark_x']
        
        #subdir = paste0(save_dir,'/',pattern_name,'/')
        save(spatial_gene_results, file = paste0(dir,"result_SV_spark_x.rda"))
        
      } else if(selected_method == 'spark_x2') 
      {
        message('Method: spark_x with covariates (BayesSpace clusters)')
        ## spark
        start = proc.time()
        sparkx_spatialgenes_sim = do.call('spark_x_WithClus', c(selected_params))
        #print(class(spark_spatialgenes_sim))
        #print(str(spark_spatialgenes_sim))
        #print(gene_name)
        sparkx_result = sparkx_spatialgenes_sim#[genes == gene_name]
        sparkx_time = proc.time() - start
        
        sparkx_result[, pres := rep_i]
        sparkx_result[, time := sparkx_time[['elapsed']] ]
        
        spatial_gene_results = sparkx_result[,.(genes, adjustedPval,combinedPval, pres, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value','p.value', 'repetition', 'time')
        spatial_gene_results[, method := 'spark_x_with_clusters']
        
        #subdir = paste0(save_dir,'/',pattern_name,'/')
        save(spatial_gene_results, file = paste0(dir,"result_SV_spark_x_with_clusters_new.rda"))
        
        
      } else if(selected_method == 'meringue') 
      {
        message('Method: MERINGUE')
        ## meringue
        start = proc.time()
        meringue_spatialgenes_sim = do.call('meringue', c(
          selected_params))
        #print(class(meringue_spatialgenes_sim))
        #print(str(meringue_spatialgenes_sim))
        #print(gene_name)
        meringue_result = meringue_spatialgenes_sim#[genes == gene_name]
        meringue_time = proc.time() - start
        
        meringue_result[, pres := rep_i]
        meringue_result[, time := meringue_time[['elapsed']] ]
        
        spatial_gene_results = meringue_result[,.(genes, p.adj,p.value, pres, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value', 'p.value','repetition', 'time')
        spatial_gene_results[, method := 'meringue'] 
        #subdir = paste0(save_dir,'/',pattern_name,'/')
        save(spatial_gene_results, file = paste0(dir,"result_SV_meringue.rda"))
      } else if(selected_method == 'nnSVG') 
      {
        message('Method: nnSVG')
        ## spark
        start = proc.time()
        print(sce_object)
        nnSVG_spatialgenes_sim = do.call('nnSVG.test', c(
          selected_params))
        #print(class(spark_spatialgenes_sim))
        #print(str(spark_spatialgenes_sim))
        #print(gene_name)
        nnSVG_result = nnSVG_spatialgenes_sim#[genes == gene_name]
        nnSVG_time = proc.time() - start
        
        nnSVG_result[, pres := rep_i]
        nnSVG_result[, time := nnSVG_time[['elapsed']] ]
        
        spatial_gene_results = nnSVG_result[,.(genes, padj, pval, pres, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value', 'p.value', 'repetition', 'time')
        spatial_gene_results[, method := 'nnSVG']
        
        #subdir = paste0(save_dir,'/',pattern_name,'/')
        save(spatial_gene_results, file = paste0(dir,"result_SV_nnSVG.rda"))
      } else if(selected_method == 'nnSVG2') 
      {
        message('Method: nnSVG with BayesSpace clusters')
        ## spark
        start = proc.time()
        print(sce_object)
        nnSVG_spatialgenes_sim = do.call('nnSVG.testWithClus', c(
          selected_params))
        #print(class(spark_spatialgenes_sim))
        #print(str(spark_spatialgenes_sim))
        #print(gene_name)
        nnSVG_result = nnSVG_spatialgenes_sim#[genes == gene_name]
        nnSVG_time = proc.time() - start
        
        nnSVG_result[, pres := rep_i]
        nnSVG_result[, time := nnSVG_time[['elapsed']] ]
        
        spatial_gene_results = nnSVG_result[,.(genes, padj, pval, pres, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value', 'p.value', 'repetition', 'time')
        spatial_gene_results[, method := 'nnSVG_with_clusters']
        
        #subdir = paste0(save_dir,'/',pattern_name,'/')
        save(spatial_gene_results, file = paste0(dir,"result_SV_nnSVG_with_clusters_new.rda"))
        
      }
      
      result_list[[test]] = spatial_gene_results
    }
    
    save(result_list, file = paste0(dir,"/result_list_all.rda"))
    results = data.table::rbindlist(l = result_list)
    return(results)
    
  } else {
    return(NULL)
  }
}


#' @title runPatternSimulation
#' @name runPatternSimulation
#' @description Creates a known spatial pattern for selected genes one-by-one and runs the different spatial gene detection tests
#' @param gobject giotto object
#' @param pattern_name name of spatial pattern
#' @param pattern_colors 2 color vector for the spatial pattern
#' @param pattern_cell_ids cell ids that make up the spatial pattern
#' @param gene_names selected genes
#' @param spatial_probs probabilities to test for a high expressing gene value to be part of the spatial pattern
#' @param reps number of random simulation repetitions
#' @param spatial_network_name which spatial network to use for binSpectSingle
#' @param spat_methods vector of spatial methods to test
#' @param spat_methods_params list of parameters list for each element in the vector of spatial methods to test
#' @param spat_methods_names name for each element in the vector of spatial elements to test
#' @param scalefactor library size scaling factor when re-normalizing dataset
#' @param save_plot save intermediate random simulation plots or not
#' @param save_raw save the raw expression matrix of the simulation
#' @param save_norm save the normalized expression matrix of the simulation
#' @param save_dir directory to save results to
#' @param max_col maximum number of columns for final plots
#' @param height height of final plots
#' @param width width of final plots
#' @param run_simulations run simulations (default = TRUE)
#' @param \dots additional parameters for renormalization
#' @return data.table with results
#' @export
runPatternSimulation = function(sce_object, rawobject = object_subset,
                                coordinate_name = coordinate_name,
                                pattern_name = 'pattern',
                                pattern_colors = c('in' = 'green', 'out' = 'red'),
                                pattern_cell_ids = NULL,
                                gene_names = NULL,
                                spatial_probs = c(0.5, 1),
                                reps = 2,
                                SV_percent = 0.33,
                                spatial_network_name = 'kNN_network',
                                spat_methods = c( 'spatialDE', 'spark'),
                                spat_methods_params = list( NA, NA),
                                spat_methods_names = c( 'spatialDE', 'spark'),
                                layer_names=layer_names,
                                scalefactor = 6000,
                                save_plot = T,
                                save_data = TRUE,
                                save_raw = T,
                                save_norm = T,
                                save_dir = '~',
                                max_col = 4,
                                height = 7,
                                width = 7,
                                run_simulations = TRUE,
                                verbose = TRUE,
                                original_layer_names = "spatial.cluster",
                                default = T,
                                cluster_index=NULL,
                                platform = platform,
                                ...) {
  set.seed(123)
  if(length(SV_percent)>1 && sum(SV_percent)!=1 )stop('sum of SV_percent needs to be one')
  
  if(!is.null(cluster_index)){
    if(spatial_probs[1] == 0.5 && 
       (length(spatial_probs) != length(SV_percent) | length(SV_percent) != (length(cluster_index)+1))){
      stop('lengths of spatial probabilities, SV_percent and number of spatial mix patterns need to be the same')
    }
    
    if(spatial_probs[1] != 0.5 && 
       (length(spatial_probs) != length(SV_percent) | length(SV_percent) != (length(cluster_index)) ))
    {stop('lengths of spatial probabilities, SV_percent and number of spatial mix patterns need to be the same')
    }}
  
  if(is.null(cluster_index)){
    if(length(spatial_probs) != (length(SV_percent)+1) )
    {stop('lengths of spatial probabilities, SV_percent and number of spatial mix patterns need to be the same')}
}
  
  # data.table variables
  prob = method = adj.p.value = time = NULL
  simobject = sce_object
  subdir = paste0(save_dir,'/',pattern_name,'/')
  resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
  if(!file.exists(resultdir)) dir.create(path = resultdir, recursive = TRUE)
  all_results = list()
  
  rep_list = list()
  save_name = 'probs'
  for(rep_i in 1:reps) {
    
    
    if(verbose) cat('\n \n repetitiion = ', rep_i, '\n \n')
    if(file.exists(paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))) {
      message("Simulation dataset has been created. 
            Skip data simulation part.
            Load the simulation dataset...")
      load(paste0(resultdir, save_name,'_',spatial_probs[1],'_',spatial_probs[2],'_final_object.rda'))
      if(exists('final_object')) 
        final_object <-  final_object
      else 
        final_object <- combined_object
    }else{
      final_object = CombineSimulatedObject(sce_object,
                                            coordinate_name = coordinate_name,
                                            pattern_name = pattern_name,
                                            pattern_cell_ids = pattern_cell_ids,
                                            gene_name = gene_names,
                                            spatial_probs = spatial_probs,
                                            rep_i = rep_i,
                                            SV_percent = SV_percent,
                                            layer_names=layer_names,
                                            save_plot = FALSE,
                                            save_raw = TRUE,
                                            save_norm = TRUE,
                                            save_dir = save_dir,
                                            verbose = TRUE,
                                            run_simulations = TRUE,
                                            simobject,
                                            default = default)
    }
    print(dim(assays(final_object)$counts))
    print(dim(assays(rawobject)$counts))
    generesults = run_spatial_sim_tests_methods(sce_object = final_object,
                                                coordinate_name = coordinate_name,
                                                rep_i = rep_i,
                                                raw_sce = rawobject,
                                                pattern_name = pattern_name,
                                                #pattern_cell_ids = pattern_cell_ids,
                                                save_dir = save_dir,
                                                gene_name = gene_names,
                                                spatial_prob = spatial_probs,
                                                gradient_direction = NULL,
                                                show_pattern = show_pattern,
                                                layer_names=layer_names,
                                                spat_methods = spat_methods,
                                                spat_methods_params = spat_methods_params,
                                                spat_methods_names = spat_methods_names,
                                                original_layer_names = original_layer_names,
                                                dir = resultdir,
                                                default = default,
                                                platform = platform)
    
    if(run_simulations == TRUE) {
      #generesults[, prob := as.factor(prob)]
      uniq_methods = sort(unique(generesults$method))
      generesults[, method := factor(method, levels = uniq_methods)]
      
      if(save_data == TRUE) {
        
        subdir = paste0(save_dir,'/',pattern_name,'/')
        resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
        if(!file.exists(subdir)) dir.create(path = subdir, recursive = TRUE)
        # write results
        data.table::fwrite(x = generesults, file = paste0(resultdir,rep_i,'_probs_',spatial_probs[1],'_',spatial_probs[2],'_results.txt'), sep = '\t', quote = F)
        
      }
      
      all_results[[rep_i]] = generesults
      
    }
    
  }
  
  
  ## create combined results and visuals
  if(run_simulations == TRUE) {
    
    results = do.call('rbind', all_results)
    
  }}

runStLearnPatternSimulation = function(sce_object, rawobject = object_subset,
                                       coordinate_name = coordinate_name,
                                       st_result_path = paste0('~/Desktop/stLearn/',dataset, '/',pattern_name,'/probs_',
                                                               spatial_probs[1],'_',spatial_probs[2],default,'/',
                                                               sample_names[i],"_stLearn_results.csv"),
                                       pattern_name = 'pattern',
                                       pattern_colors = c('in' = 'green', 'out' = 'red'),
                                       pattern_cell_ids = NULL,
                                       gene_names = NULL,
                                       spatial_probs = c(0.5, 1),
                                       reps = 2,
                                       spatial_network_name = 'kNN_network',
                                       spat_methods = c( 'spatialDE', 'spark'),
                                       spat_methods_params = list( NA, NA),
                                       spat_methods_names = c( 'spatialDE', 'spark'),
                                       layer_names=layer_names,
                                       scalefactor = 6000,
                                       save_plot = T,
                                       save_data = TRUE,
                                       save_raw = T,
                                       save_norm = T,
                                       save_dir = '~',
                                       max_col = 4,
                                       height = 7,
                                       width = 7,
                                       run_simulations = TRUE,
                                       verbose = TRUE,
                                       original_layer_names = "spatial.cluster",
                                       default = T,
                                       cluster_method = 'stLearn',
                                       dataset = 'LIBD',
                                       platform = platform,
                                       ...) {
  
  set.seed(123)
  # data.table variables
  prob = method = adj.p.value = time = NULL
  simobject = sce_object
  #raw_gobject = gobject
  # plot pattern for first gene (the same for all)
  
  subdir = paste0(save_dir,'/',pattern_name,'/')
  resultdir <- paste0(subdir,'probs_',spatial_probs[1],'_',spatial_probs[2],default,'/')
  
  if(!file.exists(resultdir)) dir.create(path = resultdir, recursive = TRUE)
  
  all_results = list()
  
  rep_list = list()
  for(rep_i in 1:reps) {
    
    
    if(verbose) cat('\n \n repetitiion = ', rep_i, '\n \n')
    
    load(paste0(resultdir,"probs_",spatial_probs[1],'_',spatial_probs[2],"_final_object.rda"))
    if(exists("final_object1")) {final_object = final_object1}
    if(exists("final_object")) {final_object = final_object}
    if(exists("combined_object")) {final_object = combined_object}
    print(dim(final_object))
    print(dim(rawobject))
    print(original_layer_names)
    
    stLearn_results <- read.csv(st_result_path, sep = ',', header = TRUE)
    colnames(stLearn_results) <- c("cell_ID","imagecol","imagerow","stLearn_pca_kmeans")
    
    final_object <- addMetadata(
      final_object,
      stLearn_results,
      by_column = TRUE,
      column_cell_ID = "cell_ID"
    )
    
    
    print(head(colData(final_object)))
    generesults = run_spatial_sim_tests_methods(sce_object = final_object,
                                                coordinate_name = coordinate_name,
                                                rep_i = rep_i,
                                                raw_sce = rawobject,
                                                pattern_name = pattern_name,
                                                pattern_cell_ids = pattern_cell_ids,
                                                save_dir = save_dir,
                                                gene_name = gene_names,
                                                spatial_prob = spatial_probs,
                                                gradient_direction = NULL,
                                                show_pattern = show_pattern,
                                                layer_names=layer_names,
                                                spat_methods = spat_methods,
                                                spat_methods_params = spat_methods_params,
                                                spat_methods_names = spat_methods_names,
                                                original_layer_names = original_layer_names,
                                                dir = resultdir,
                                                default = default,
                                                cluster_method = cluster_method,
                                                platform = platform)
    
    if(run_simulations == TRUE) {
      #generesults[, prob := as.factor(prob)]
      uniq_methods = sort(unique(generesults$method))
      generesults[, method := factor(method, levels = uniq_methods)]
      
      if(save_data == TRUE) {
        
        subdir = paste0(save_dir,'/',pattern_name,'/')
        if(!file.exists(subdir)) dir.create(path = subdir, recursive = TRUE)
        # write results
        data.table::fwrite(x = generesults, file = paste0(subdir,'/',rep_i,'_probs_',spatial_probs[1],'_',spatial_probs[2],layer_names,'_results.txt'), sep = '\t', quote = F)
        
      }
      
      all_results[[rep_i]] = generesults
      
    }
    
  }
  
  
  ## create combined results and visuals
  if(run_simulations == TRUE) {
    
    results = do.call('rbind', all_results)
    
  }}
