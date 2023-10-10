library(data.table)
library(iCOBRA)
library(S4Vectors)
library(stringr)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(pROC)
########################## real data analysis #############################
my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
process_data <- function(
    path = '~',
    data = 'melanoma', # or 'LIBD'
    Manual = F, # if there are results from 'Manual clusters'
    threshold = 0.05,
    name_head = -13,
    name_tail = -5,
    ...
){
  file = paste0(path,data)
  edgeR_counts_BayesSpace_files <- list.files(
    file, pattern = "BayesSpace_edgeR_results *.*csv$", full.names = TRUE)
  edgeR_counts_BayesSpace_results <- lapply(edgeR_counts_BayesSpace_files, read.csv)
  
  edgeR_counts_stLearn_files <- list.files(
    file, pattern = "^stLearn_edgeR_results *.*csv$", full.names = TRUE)
  edgeR_counts_stLearn_results <- lapply(edgeR_counts_stLearn_files, read.csv)
  
  scran_BayesSpace_files <- list.files(
    file, pattern =   "BayesSpace_scran_results.rda$", full.names = TRUE)
  scran_BayesSpace_results <- lapply(scran_BayesSpace_files, 
                                     function(xx) {
                                       load(xx, envir=globalenv())
                                       BayesSpace_scran_results <- Reduce(
                                         function(x, y) merge(x, y, by = 'genes', all = T),
                                         lapply(scranResults, function(x) { x$genes <- rownames(x); x }))

                                       # BayesSpace_scran_results <- do.call(cbind, scranResults)
                                       # WRONG: orders of genes of each element in the list are different
                                       rm(scranResults)
                                       ind1 <- which(sub("\\..*", "", colnames(BayesSpace_scran_results)) == "p")
                                       ind2 <- which(sub("\\..*", "", colnames(BayesSpace_scran_results)) == "FDR")
                                       ## Check whether types of p-val and FDR columns are numeric instead of character -> change the minimum value
                                       if(sum(sapply(BayesSpace_scran_results[ind1], class) != 'numeric') != 0){
                                         BayesSpace_scran_results[ind1] <- sapply(BayesSpace_scran_results[ind1], as.numeric)
                                         BayesSpace_scran_results[ind2] <- sapply(BayesSpace_scran_results[ind2], as.numeric)
                                       }
                                       
                                       BayesSpace_scran_pval <- apply(BayesSpace_scran_results[ind1], 1, my.min)
                                       BayesSpace_scran_FDR <- apply(BayesSpace_scran_results[ind2], 1, my.min)
                                       
                                       data.frame(genes = BayesSpace_scran_results$genes,
                                                                   adj.p.value = BayesSpace_scran_FDR,
                                                                   p.value = BayesSpace_scran_pval,
                                                                   method = "BayesSpace_findMarkers")
                                       }
                                     )
  
  scran_stLearn_files <- list.files(
    file, pattern =   "StLearn_scran_results.rda$", full.names = TRUE)
  scran_stLearn_results <- lapply(scran_stLearn_files, 
                                     function(xx) {
                                        load(xx, envir=globalenv())
                                        stLearn_scran_results <- Reduce(
                                          function(x, y) merge(x, y, by = 'genes', all = T),
                                          lapply(scranResults, function(x) { x$genes <- rownames(x); x }))
                                        #do.call(cbind, scranResults)
                                        rm(scranResults)
                                        ind1 <- which(sub("\\..*", "", colnames(stLearn_scran_results)) == "p")
                                        ind2 <- which(sub("\\..*", "", colnames(stLearn_scran_results)) == "FDR")
                                        ## Check whether types of p-val and FDR columns are numeric instead of character -> change the minimum value
                                        if(sum(sapply(stLearn_scran_results[ind1], class) != 'numeric') != 0){
                                          stLearn_scran_results[ind1] <- sapply(stLearn_scran_results[ind1], as.numeric)
                                          stLearn_scran_results[ind2] <- sapply(stLearn_scran_results[ind2], as.numeric)
                                        }
                                        stLearn_scran_pval <- apply(stLearn_scran_results[ind1], 1, my.min)
                                        stLearn_scran_FDR <- apply(stLearn_scran_results[ind2], 1, my.min)
                                        
                                        data.frame(genes = stLearn_scran_results$genes,
                                                                            adj.p.value = stLearn_scran_FDR,
                                                                            p.value = stLearn_scran_pval,
                                                                            method = "StLearn_findMarkers")
                                     })
  

  seurat_BayesSpace_files <- list.files(
    file, pattern =   "BayesSpace_seurat_results.rda$", full.names = TRUE)
  seurat_BayesSpace_results <- lapply(seurat_BayesSpace_files, 
                                     function(xx) {
                                       load(xx, envir=globalenv())
                                       BayesSpace_seurat_results <- as.data.frame(seuratResults)
                                       rm(seuratResults)
                                       BayesSpace_seurat <- stats::reshape(BayesSpace_seurat_results, idvar = "gene", 
                                                                 v.names = c("p_val", "p_val_adj"), 
                                                                 timevar = "cluster", direction = "wide")
                                       
                                       ind1 <- which(sub("\\..*", "", colnames(BayesSpace_seurat)) == "p_val")
                                       ind2 <- which(sub("\\..*", "", colnames(BayesSpace_seurat)) == "p_val_adj")
                                       ## Check whether types of p-val and FDR columns are numeric instead of character -> change the minimum value
                                       if(sum(sapply(BayesSpace_seurat[ind1], class) != 'numeric') != 0){
                                         BayesSpace_seurat[ind1] <- sapply(BayesSpace_seurat[ind1], as.numeric)
                                         BayesSpace_seurat[ind2] <- sapply(BayesSpace_seurat[ind2], as.numeric)
                                       }
                                       BayesSpace_seurat_pval <- apply(BayesSpace_seurat[ind1], 1, my.min)
                                       BayesSpace_seurat_FDR <- apply(BayesSpace_seurat[ind2], 1, my.min)
                                       
                                       data.frame(genes = BayesSpace_seurat$gene,
                                                  adj.p.value = BayesSpace_seurat_FDR,
                                                  p.value = BayesSpace_seurat_pval,
                                                  method = "BayesSpace_FindAllMarkers")
                                     }
  )
  
  seurat_stLearn_files <- list.files(
    file, pattern =   "StLearn_seurat_results.rda$", full.names = TRUE)
  seurat_stLearn_results <- lapply(seurat_stLearn_files, 
                                  function(xx) {
                                    load(xx, envir=globalenv())
                                    stLearn_seurat_results <- as.data.frame( seuratResults)
                                    rm(seuratResults)
                                    stLearn_seurat <- stats::reshape(stLearn_seurat_results, idvar = "gene", 
                                                              v.names = c("p_val", "p_val_adj"), 
                                                              timevar = "cluster", direction = "wide")
                                    
                                    ind1 <- which(sub("\\..*", "", colnames(stLearn_seurat)) == "p_val")
                                    ind2 <- which(sub("\\..*", "", colnames(stLearn_seurat)) == "p_val_adj")
                                    ## Check whether types of p-val and FDR columns are numeric instead of character -> change the minimum value
                                    if(sum(sapply(stLearn_seurat[ind1], class) != 'numeric') != 0){
                                      stLearn_seurat[ind1] <- sapply(stLearn_seurat[ind1], as.numeric)
                                      stLearn_seurat[ind2] <- sapply(stLearn_seurat[ind2], as.numeric)
                                    }
                                    stLearn_seurat_pval <- apply(stLearn_seurat[ind1], 1, my.min)
                                    stLearn_seurat_FDR <- apply(stLearn_seurat[ind2], 1, my.min)
                                    
                                    data.frame(genes = stLearn_seurat$gene,
                                               adj.p.value = stLearn_seurat_FDR,
                                               p.value = stLearn_seurat_pval,
                                               method = "StLearn_FindAllMarkers")
                                  })
  
  
  MERINGUE_files <- list.files(
    file, pattern = "^meringue_results *.*csv$", full.names = TRUE)
  MERINGUE_results <- lapply(MERINGUE_files, read.csv)
  
  nnSVG_files <- list.files(
    file, pattern = "^nnSVG_results *.*csv$", full.names = TRUE)
  nnSVG_results <- lapply(nnSVG_files, read.csv)
  # 
  SPARK_files <- list.files(
    file, pattern = "^spark_results *.*csv$", full.names = TRUE)
  SPARK_results <- lapply(SPARK_files, read.csv)
  # 
  SPARKX_files <- list.files(
    file, pattern = "^sparkx_results *.*csv$", full.names = TRUE)
  SPARKX_results <- lapply(SPARKX_files, read.csv)
  
  spatialDE_files <- list.files(
    file, pattern = "*_spatialDE_results *.*csv$", full.names = TRUE)
  spatialDE_results <- lapply(spatialDE_files, read.csv)
  
  spatialDE2_files <- list.files(
    file, pattern = "*_spatialDE2_results *.*csv$", full.names = TRUE)
  spatialDE2_results <- lapply(spatialDE2_files, read.csv)
  
  spaGCN_files <- list.files(
    file, pattern = "*_SpaGCN_results.csv$", full.names = TRUE)
  spaGCN_results <- lapply(spaGCN_files, read.csv)
  
  #message("Check sample names for each mathod: (SPARK for LIBD only contains the first sample 151507).")
  (names_edgeR_counts_BayesSpace <- str_sub(edgeR_counts_BayesSpace_files,name_head,name_tail))
  (names_edgeR_counts_stLearn <- str_sub(edgeR_counts_stLearn_files,name_head,name_tail))
  (names_MERINGUE <- str_sub(MERINGUE_files,name_head,name_tail))
  (names_SPARK <- str_sub(SPARK_files,name_head,name_tail))
  (names_SPARKX <- str_sub(SPARKX_files,name_head,name_tail))
  (names_nnSVG <- str_sub(nnSVG_files,name_head,name_tail))
  (names_seurat_BayesSpace <- str_sub(seurat_BayesSpace_files,-36,-31))
  (names_seurat_stLearn <- str_sub(seurat_stLearn_files,-33,-28))
  (names_scran_BayesSpace <- str_sub(scran_BayesSpace_files,-35,-30))
  (names_scran_stLearn <- str_sub(scran_stLearn_files,-32,-27))
  
  if(data == 'LIBD'){
    (names_spatialDE <- str_sub(spatialDE_files,-28,-23))
    (names_spatialDE2 <- str_sub(spatialDE2_files,-29,-24))
    (names_spaGCN <- str_sub(spaGCN_files,-25,-20))
  }else{
    (names_spatialDE <- str_sub(spatialDE_files,-31,-23))
    (names_spatialDE2 <- str_sub(spatialDE2_files,-32,-24))
    (names_spaGCN <- str_sub(spaGCN_files,-28,-20))
  }
  
  names(edgeR_counts_BayesSpace_results) <- names_edgeR_counts_BayesSpace
  names(edgeR_counts_stLearn_results) <- names_edgeR_counts_stLearn
  names(MERINGUE_results) <- names_MERINGUE
  names(SPARK_results) <- names_SPARK
  names(SPARKX_results) <- names_SPARKX
  names(nnSVG_results) <- names_nnSVG
  names(spatialDE_results) <- names_spatialDE
  names(spatialDE2_results) <- names_spatialDE2
  names(spaGCN_results) <- names_spaGCN
  names(seurat_BayesSpace_results) <- names_seurat_BayesSpace
  names(seurat_stLearn_results) <- names_seurat_stLearn
  names(scran_BayesSpace_results) <- names_scran_BayesSpace
  names(scran_stLearn_results) <- names_scran_stLearn
  
  colnames <- c("pval","genes","method")
  edgeR_counts_BayesSpace_sig <- lapply(edgeR_counts_BayesSpace_results, function(d) {
    d <- d[complete.cases(d),]
    d <- d[order(d$FDR, decreasing = FALSE), ] 
    d[d$FDR< threshold, ]$genes
  })
  
  l <- lapply(edgeR_counts_BayesSpace_results, "[",  c("PValue","genes"))
  l2 <-lapply(l, function(x) 
    cbind(x, method = "BayesSpace_edgeR"))
  l2 <- lapply(l2, setNames, colnames)
  l2 <- bind_rows(l2, .id = "sample")
  edgeR_counts_BayesSpace_results2 <- l2
  
  edgeR_counts_stLearn_sig <- lapply(edgeR_counts_stLearn_results, function(d) {
    d <- d[complete.cases(d),]
    d <- d[order(d$FDR, decreasing = FALSE), ] 
    d[d$FDR< threshold, ]$genes
  })
  
  l <- lapply(edgeR_counts_stLearn_results, "[",  c("PValue","genes"))
  l2 <-lapply(l, function(x) 
    cbind(x, method = "StL
          arn_edgeR"))
  l2 <- lapply(l2, setNames, colnames)
  l2 <- bind_rows(l2, .id = "sample")
  edgeR_counts_stLearn_results2 <- l2
  
  
  MERINGUE_sig <- lapply(MERINGUE_results, function(d) {
    d <- d[complete.cases(d),]
    d <- d[order(d$p.adj, decreasing = FALSE), ] 
    d[d$p.adj< threshold, ]$X
  })
  
  l <- lapply(MERINGUE_results, "[",  c("p.value","X"))
  l2 <-lapply(l, function(x) 
    cbind(x, method = "MERINGUE"))
  l2 <- lapply(l2, setNames, colnames)
  l2 <- bind_rows(l2, .id = "sample")
  MERINGUE_results2 <- l2
  
  
  SPARK_sig <- lapply(SPARK_results, function(d) {
    d <- d[complete.cases(d),]
    d <- d[order(d$adjusted_pvalue, decreasing = FALSE), ]
    d[d$adjusted_pvalue< threshold, ]$X
  })
  l <- lapply(SPARK_results, "[",  c("combined_pvalue","X"))
  l2 <-lapply(l, function(x) 
    cbind(x, method = "SPARK"))
  l2 <- lapply(l2, setNames, colnames)
  l2 <- bind_rows(l2, .id = "sample")
  SPARK_results2 <- l2
  
  SPARKX_sig <- lapply(SPARKX_results, function(d) {
    d <- d[complete.cases(d),]
    d <- d[order(d$adjustedPval, decreasing = FALSE), ] 
    d[d$adjustedPval< threshold, ]$X
  })
  l <- lapply(SPARKX_results, "[",  c("combinedPval","X"))
  l2 <-lapply(l, function(x) 
    cbind(x, method = "SPARKX"))
  l2 <- lapply(l2, setNames, colnames)
  l2 <- bind_rows(l2, .id = "sample")
  SPARKX_results2 <- l2
  
  nnSVG_sig <- lapply(nnSVG_results, function(d) {
    d <- d[complete.cases(d),]
    d <- d[order(d$padj, decreasing = FALSE), ]
    d[d$padj< threshold, ]$X
  })
  l <- lapply(nnSVG_results, "[",  c("pval","X"))
  l2 <-lapply(l, function(x)
    cbind(x, method = "nnSVG"))
  l2 <- lapply(l2, setNames, colnames)
  l2 <- bind_rows(l2, .id = "sample")
  nnSVG_results2 <- l2
  
  spatialDE_sig <- lapply(spatialDE_results, function(d) {
    d <- d[complete.cases(d),]
    d <- d[order(d$qval, decreasing = FALSE), ] 
    d[d$qval < threshold, ]$g
  })
  l <- lapply(spatialDE_results, "[",  c("pval","g"))
  l2 <-lapply(l, function(x) 
    cbind(x, method = "SpatialDE"))
  l2 <- lapply(l2, setNames, colnames)
  l2 <- bind_rows(l2, .id = "sample")
  SpatialDE_results2 <- l2
  
  spatialDE2_sig <- lapply(spatialDE2_results, function(d) {
    d <- d[complete.cases(d),]
    d <- d[order(d$padj, decreasing = FALSE), ] 
    d[d$padj < threshold, ]$gene
  })
  l <- lapply(spatialDE2_results, "[",  c("pval","gene"))
  l2 <-lapply(l, function(x) 
    cbind(x, method = "SpatialDE2"))
  l2 <- lapply(l2, setNames, colnames)
  l2 <- bind_rows(l2, .id = "sample")
  SpatialDE2_results2 <- l2
  
  
  spaGCN_sig <- lapply(spaGCN_results, function(d) {
    d <- d[complete.cases(d),]
    d <- d[order(d$min_pvals_adj, decreasing = FALSE), ] 
    d[d$min_pvals_adj < threshold, ]$genes
  })
  l <- lapply(spaGCN_results, "[",  c("min_pvals_adj","genes"))
  l2 <-lapply(l, function(x) 
    cbind(x, method = "SpaGCN"))
  l2 <- lapply(l2, setNames, colnames)
  l2 <- bind_rows(l2, .id = "sample")
  SpaGCN_results2 <- l2
  
  seurat_BayesSpace_sig <- lapply(seurat_BayesSpace_results, function(d) {
    d <- d[complete.cases(d),]
    d <- d[order(d$adj.p.value, decreasing = FALSE), ] 
    d[d$adj.p.value < threshold, ]$genes
  })
  l <- lapply(seurat_BayesSpace_results, "[",  c("p.value","genes"))
  l2 <-lapply(l, function(x) 
    cbind(x, method = "BayesSpace_FindAllMarkers"))
  l2 <- lapply(l2, setNames, colnames)
  l2 <- bind_rows(l2, .id = "sample")
  seurat_BayesSpace_results2 <- l2
  
  seurat_stLearn_sig <- lapply(seurat_stLearn_results, function(d) {
    d <- d[complete.cases(d),]
    d <- d[order(d$adj.p.value, decreasing = FALSE), ] 
    d[d$adj.p.value < threshold, ]$genes
  })
  l <- lapply(seurat_stLearn_results, "[",  c("p.value","genes"))
  l2 <-lapply(l, function(x) 
    cbind(x, method = "StLearn_FindAllMarkers"))
  l2 <- lapply(l2, setNames, colnames)
  l2 <- bind_rows(l2, .id = "sample")
  seurat_stLearn_results2 <- l2
  
  scran_BayesSpace_sig <- lapply(scran_BayesSpace_results, function(d) {
    d <- d[complete.cases(d),]
    d <- d[order(d$adj.p.value, decreasing = FALSE), ] 
    d[d$adj.p.value < threshold, ]$genes
  })
  l <- lapply(scran_BayesSpace_results, "[",  c("p.value","genes"))
  l2 <-lapply(l, function(x) 
    cbind(x, method = "BayesSpace_findMarkers"))
  l2 <- lapply(l2, setNames, colnames)
  l2 <- bind_rows(l2, .id = "sample")
  scran_BayesSpace_results2 <- l2
  
  scran_stLearn_sig <- lapply(scran_stLearn_results, function(d) {
    d <- d[complete.cases(d),]
    d <- d[order(d$adj.p.value, decreasing = FALSE), ] 
    d[d$adj.p.value < threshold, ]$genes
  })
  l <- lapply(scran_stLearn_results, "[",  c("p.value","genes"))
  l2 <-lapply(l, function(x) 
    cbind(x, method = "StLearn_findMarkers"))
  l2 <- lapply(l2, setNames, colnames)
  l2 <- bind_rows(l2, .id = "sample")
  scran_stLearn_results2 <- l2
  
  if(Manual == TRUE){
    edgeR_counts_Manual_files <- list.files(
      file, pattern = "Manual_edgeR_results *.*csv$", full.names = TRUE)
    edgeR_counts_Manual_results <- lapply(edgeR_counts_Manual_files, read.csv)
    (names_edgeR_counts_Manual <- str_sub(edgeR_counts_Manual_files,name_head,name_tail))
    names(edgeR_counts_Manual_results) <- names_edgeR_counts_Manual
    edgeR_counts_Manual_sig <- lapply(edgeR_counts_Manual_results, function(d) {
      d <- d[complete.cases(d),]
      d <- d[order(d$FDR, decreasing = FALSE), ] 
      d[d$FDR< threshold, ]$genes
    })
    
    l <- lapply(edgeR_counts_Manual_results, "[",  c("PValue","genes"))
    l2 <-lapply(l, function(x) 
      cbind(x, method = "Manual_edgeR"))
    l2 <- lapply(l2, setNames, colnames)
    l2 <- bind_rows(l2, .id = "sample")
    edgeR_counts_Manual_results2 <- l2
    
    seurat_Manual_files <- list.files(
      file, pattern =   "Manual_seurat_results.rda$", full.names = TRUE)
    seurat_Manual_results <- lapply(seurat_Manual_files, 
                                     function(xx) {
                                       load(xx, envir=globalenv())
                                       Manual_seurat_results <- as.data.frame( seuratResults)
                                       rm(seuratResults)
                                       Manual_seurat <- stats::reshape(Manual_seurat_results, idvar = "gene", 
                                                                        v.names = c("p_val", "p_val_adj"), 
                                                                        timevar = "cluster", direction = "wide")
                                       
                                       ind1 <- which(sub("\\..*", "", colnames(Manual_seurat)) == "p_val")
                                       ind2 <- which(sub("\\..*", "", colnames(Manual_seurat)) == "p_val_adj")
                                       ## Check whether types of p-val and FDR columns are numeric instead of character -> change the minimum value
                                       if(sum(sapply(Manual_seurat[ind1], class) != 'numeric') != 0){
                                         Manual_seurat[ind1] <- sapply(Manual_seurat[ind1], as.numeric)
                                         Manual_seurat[ind2] <- sapply(Manual_seurat[ind2], as.numeric)
                                       }
                                       Manual_seurat_pval <- apply(Manual_seurat[ind1], 1, my.min)
                                       Manual_seurat_FDR <- apply(Manual_seurat[ind2], 1, my.min)
                                       
                                       data.frame(genes = Manual_seurat$gene,
                                                  adj.p.value = Manual_seurat_FDR,
                                                  p.value = Manual_seurat_pval,
                                                  method = "Manual_FindAllMarkers")
                                     })
    
    scran_Manual_files <- list.files(
      file, pattern =   "Manual_scran_results.rda$", full.names = TRUE)
    scran_Manual_results <- lapply(scran_Manual_files, 
                                    function(xx) {
                                      load(xx, envir=globalenv())
                                      Manual_scran_results <- Reduce(
                                        function(x, y) merge(x, y, by = 'genes', all = T),
                                        lapply(scranResults, function(x) { x$genes <- rownames(x); x }))
                                      #do.call(cbind, scranResults)
                                      rm(scranResults)
                                      ind1 <- which(sub("\\..*", "", colnames(Manual_scran_results)) == "p")
                                      ind2 <- which(sub("\\..*", "", colnames(Manual_scran_results)) == "FDR")
                                      ## Check whether types of p-val and FDR columns are numeric instead of character -> change the minimum value
                                      if(sum(sapply(Manual_scran_results[ind1], class) != 'numeric') != 0){
                                        Manual_scran_results[ind1] <- sapply(Manual_scran_results[ind1], as.numeric)
                                        Manual_scran_results[ind2] <- sapply(Manual_scran_results[ind2], as.numeric)
                                      }
                                      Manual_scran_pval <- apply(Manual_scran_results[ind1], 1, my.min)
                                      Manual_scran_FDR <- apply(Manual_scran_results[ind2], 1, my.min)
                                      
                                      data.frame(genes = Manual_scran_results$genes,
                                                 adj.p.value = Manual_scran_FDR,
                                                 p.value = Manual_scran_pval,
                                                 method = "Manual_findMarkers")
                                    })
    
    (names_seurat_Manual <- str_sub(seurat_Manual_files,-32,-27))
    (names_scran_Manual <- str_sub(scran_Manual_files,-31,-26))
    names(seurat_Manual_results) <- names_seurat_Manual
    names(scran_Manual_results) <- names_scran_Manual
    
    seurat_Manual_sig <- lapply(seurat_Manual_results, function(d) {
      d <- d[complete.cases(d),]
      d <- d[order(d$adj.p.value, decreasing = FALSE), ] 
      d[d$adj.p.value < threshold, ]$genes
    })
    l <- lapply(seurat_Manual_results, "[",  c("p.value","genes"))
    l2 <-lapply(l, function(x) 
      cbind(x, method = "Manual_FindAllMarkers"))
    l2 <- lapply(l2, setNames, colnames)
    l2 <- bind_rows(l2, .id = "sample")
    seurat_Manual_results2 <- l2
    
    scran_Manual_sig <- lapply(scran_Manual_results, function(d) {
      d <- d[complete.cases(d),]
      d <- d[order(d$adj.p.value, decreasing = FALSE), ] 
      d[d$adj.p.value < threshold, ]$genes
    })
    l <- lapply(scran_Manual_results, "[",  c("p.value","genes"))
    l2 <-lapply(l, function(x) 
      cbind(x, method = "Manual_findMarkers"))
    l2 <- lapply(l2, setNames, colnames)
    l2 <- bind_rows(l2, .id = "sample")
    scran_Manual_results2 <- l2
    
    all_results <- rbind(edgeR_counts_Manual_results2,
                         edgeR_counts_BayesSpace_results2, edgeR_counts_stLearn_results2,
                         scran_Manual_results2,
                         scran_BayesSpace_results2, scran_stLearn_results2,
                         seurat_Manual_results2,
                         seurat_BayesSpace_results2, seurat_stLearn_results2,
                         MERINGUE_results2, 
                         SPARK_results2,nnSVG_results2, 
                         SPARKX_results2, SpatialDE2_results2, SpatialDE_results2, SpaGCN_results2)
    
    list_sig = list(edgeR_counts_Manual_sig,
                    edgeR_counts_BayesSpace_sig, edgeR_counts_stLearn_sig,
                    scran_Manual_sig,
                    scran_BayesSpace_sig, scran_stLearn_sig,
                    seurat_Manual_sig,
                    seurat_BayesSpace_sig, seurat_stLearn_sig,
                    MERINGUE_sig, 
                    SPARK_sig,nnSVG_sig, 
                    SPARKX_sig, spatialDE2_sig, spatialDE_sig, spaGCN_sig,
                    all_results)
  }else{
    all_results <- rbind(edgeR_counts_BayesSpace_results2, edgeR_counts_stLearn_results2,
                         scran_BayesSpace_results2, scran_stLearn_results2,
                         seurat_BayesSpace_results2, seurat_stLearn_results2,
                         MERINGUE_results2, 
                         SPARK_results2,nnSVG_results2, 
                         SPARKX_results2, SpatialDE2_results2, SpatialDE_results2, SpaGCN_results2)
    
    list_sig = list(edgeR_counts_BayesSpace_sig, edgeR_counts_stLearn_sig,
                    scran_BayesSpace_sig, scran_stLearn_sig,
                    seurat_BayesSpace_sig, seurat_stLearn_sig,
                    MERINGUE_sig, 
                    SPARK_sig,nnSVG_sig, 
                    SPARKX_sig, spatialDE2_sig, spatialDE_sig, spaGCN_sig,
                    all_results)
  }
  
  return(list_sig)
}


jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

# Get upper triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

jaccord_btw_rep <- function(
    list_test,
    path = '~',
    save_name = "_jaccard_index_corrplot.pdf",
    sample_names=sample_names,
    method,
    #Manual = F,
    top_genes = 200,
    ...
){
  list_test <- lapply(list_test,head,n=top_genes)
  # check
  stopifnot(length(list_test[1]) <= top_genes)
  result_list <- NULL
  for(i in c(1:length(sample_names))){
    row = mapply(jaccard,list_test,list_test[i])
    result_list[[i]] = as.data.frame(row)
  }
  mat_jaccord <- dplyr::bind_cols(l = result_list)
  
  colnames(mat_jaccord) <- sample_names
  rownames(mat_jaccord) <- colnames(mat_jaccord)
  lower_tri <- get_lower_tri(mat_jaccord)
  melted_cormat <- reshape2::melt(as.matrix(lower_tri), na.rm = TRUE)
  melted_cormat$Var1 <- as.character(melted_cormat$Var1)
  melted_cormat$Var2 <- as.character(melted_cormat$Var2)
  
  # Heatmap
  title = paste0("Method_", method)
  pdf(file=paste0(path,"/",title,save_name), onefile=FALSE) # or other device
  corrplot::corrplot(as.matrix(mat_jaccord),type="upper",tl.col="black",
                     mar=c(0,0,1,0), col.lim = c(0,1),tl.cex = 0.8)
  
  dev.off()
}

## compare to `ave_jaccard_btw_rep_melanoma`, 
# `ave_jaccard_btw_rep_LIBD` don't contain SPARK results but have 3 more results for Manual annotation (DESpace, findMarkers and FindAllMarkers)
ave_jaccard_btw_rep_LIBD <- function(
    list_test_all,
    path = '/',
    #save_name = "_jaccard_index_corrplot.pdf",
    sample_names=sample_names,
    top_genes = 200,
    #width = 12, height = 10,
    data = 'melanoma',
    col_method = colors_method,
    #Manual = F,
    ...
){
  for(j in c(1:10, 12:16)){
    list_test = list_test_all[[j]]
    list_test <- lapply(list_test,head,n=top_genes)
    # check
    stopifnot(length(list_test[1]) <= top_genes)
    result_list <- NULL
    for(i in c(1:length(sample_names))){
      row = mapply(jaccard,list_test,list_test[i])
      result_list[[i]] = as.data.frame(row)
    }
    mat_jaccord <- dplyr::bind_cols(l = result_list)
    
    colnames(mat_jaccord) <- sample_names
    rownames(mat_jaccord) <- colnames(mat_jaccord)
    lower_tri <- get_lower_tri(mat_jaccord)
    melted_cormat <- reshape2::melt(as.matrix(lower_tri), na.rm = TRUE)
    melted_cormat$Var1 <- as.character(melted_cormat$Var1)
    melted_cormat$Var2 <- as.character(melted_cormat$Var2)
    
    # Heatmap
    assign(paste0("method",j),as.matrix(mat_jaccord))
  }
  
  df <- as.data.frame(c(mean(method1[upper.tri(method1, diag = FALSE)]),
                        mean(method2[upper.tri(method2, diag = FALSE)]),
                        mean(method3[upper.tri(method3, diag = FALSE)]),
                        mean(method4[upper.tri(method4, diag = FALSE)]),
                        mean(method5[upper.tri(method5, diag = FALSE)]),
                        mean(method6[upper.tri(method6, diag = FALSE)]),
                        mean(method7[upper.tri(method7, diag = FALSE)]),
                        mean(method8[upper.tri(method8, diag = FALSE)]),
                        mean(method9[upper.tri(method9, diag = FALSE)]),
                        mean(method10[upper.tri(method10, diag = FALSE)]),
                        mean(method12[upper.tri(method12, diag = FALSE)]),
                        mean(method13[upper.tri(method13, diag = FALSE)]),
                        mean(method14[upper.tri(method14, diag = FALSE)]),
                        mean(method15[upper.tri(method15, diag = FALSE)]),
                        mean(method16[upper.tri(method16, diag = FALSE)])
                        
  ))
  
  colnames(df) <- c("Value")
  rownames(df) <- c("Manual_DESpace",
                    "BayesSpace_DESpace",
                    "StLearn_DESpace",
                    "Manual_findMarkers",
                    "BayesSpace_findMarkers",
                    "StLearn_findMarkers",
                    "Manual_FindAllMarkers",
                    "BayesSpace_FindAllMarkers",
                    "StLearn_FindAllMarkers",
                    "MERINGUE", #"SPARK",
                    "nnSVG",
                    "SPARK-X", "SpatialDE2",
                    "SpatialDE","SpaGCN")
  df <- reshape::melt(t(df))
  df$X2 <- as.factor(df$X2)
  #all_colours = c(rep("#009E73",1), rep("#CC79A7",1), rep("#0072B2",1), rep("#D55E00",1), rep("#CD853F",1))
  gg_data <- df[,2:3]
  colnames(gg_data) <- c( "Method", "value")
  breaks_plots = c(0.1,0.2,0.3,0.4,0.5,0.6)
  gg_data$Method <- factor(gg_data$Method,levels=c("Manual_DESpace",#"edgeR_CPMs_Manual",
                                                   "BayesSpace_DESpace",
                                                   "StLearn_DESpace",#"edgeR_CPMs_BayesSpace", 
                                                   "Manual_findMarkers",
                                                   "BayesSpace_findMarkers",
                                                   "StLearn_findMarkers",
                                                   "Manual_FindAllMarkers",
                                                   "BayesSpace_FindAllMarkers",
                                                   "StLearn_FindAllMarkers",
                                                   "MERINGUE", #"SPARK",
                                                   "nnSVG",
                                                   "SPARK-X", "SpatialDE2",
                                                   "SpatialDE","SpaGCN"))
  
  
  #pdf(file=paste0(path,save_name,data,".pdf"), onefile=FALSE) # or other device
  #gg_data %>%
  gg_Bar <- ggplot(gg_data,
                   #aes_string(y = "Time_Minutes", fill = "Method"),
                   aes(x = reorder(Method,-value), y = value, fill = Method)) + 
    
    #ggplot(aes_string(x = "Method", y = "value", fill = "Method")) +  
    geom_bar( stat = "identity",position = position_dodge(0.8), width = 0.7) +  theme_bw() +   xlab("") +  
    ylab("Average Jaccaed index across samples") + #facet_wrap(~Sample_id) +
    #scale_y_sqrt( breaks = c(0.1,0.15,0.25,0.4,0.5,0.6) ) + 
    #scale_y_continuous( breaks = c(0.05,0.1,0.15,0.2) ) + 
    scale_fill_manual("Method", values = col_method)+
    scale_y_continuous( breaks = breaks_plots ) +
    theme(axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1),   
          axis.text.y= element_text(size=rel(1.8),angle = 0, hjust = 1),  
          axis.title = element_text(size=15,face = "bold"),      
          panel.grid.minor = element_blank(),        
          panel.grid.major.x = element_blank(),       
          panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),        
          aspect.ratio = 1,        
          legend.position = "none",
          legend.text=element_text(size=15))
  return(gg_Bar)
}

ave_jaccard_btw_rep_melanoma <- function(
    list_test_all,
    path = '~',
    #save_name = "_jaccard_index_corrplot.pdf",
    sample_names=sample_names,
    top_genes = 200,
    #width = 12, height = 10,
    data = 'melanoma',
    col_method = colors_method,
    #Manual = F,
    ...
){
  for(j in c(1:13)){
    list_test = list_test_all[[j]]
    list_test <- lapply(list_test,head,n=top_genes)
    # check
    stopifnot(length(list_test[1]) <= top_genes)
    result_list <- NULL
    for(i in c(1:length(sample_names))){
      row = mapply(jaccard,list_test,list_test[i])
      result_list[[i]] = as.data.frame(row)
    }
    mat_jaccord <- dplyr::bind_cols(l = result_list)
    
    colnames(mat_jaccord) <- sample_names
    rownames(mat_jaccord) <- colnames(mat_jaccord)
    lower_tri <- get_lower_tri(mat_jaccord)
    melted_cormat <- reshape2::melt(as.matrix(lower_tri), na.rm = TRUE)
    melted_cormat$Var1 <- as.character(melted_cormat$Var1)
    melted_cormat$Var2 <- as.character(melted_cormat$Var2)
    
    # Heatmap
    assign(paste0("method",j),as.matrix(mat_jaccord))
  }
  
  df <- as.data.frame(c(mean(method1[upper.tri(method1, diag = FALSE)]),
                        mean(method2[upper.tri(method2, diag = FALSE)]),
                        mean(method3[upper.tri(method3, diag = FALSE)]),
                        mean(method4[upper.tri(method4, diag = FALSE)]),
                        mean(method5[upper.tri(method5, diag = FALSE)]),
                        mean(method6[upper.tri(method6, diag = FALSE)]),
                        mean(method7[upper.tri(method7, diag = FALSE)]),
                        mean(method8[upper.tri(method8, diag = FALSE)]),
                        mean(method9[upper.tri(method9, diag = FALSE)]),
                        mean(method10[upper.tri(method10, diag = FALSE)]),
                        mean(method11[upper.tri(method11, diag = FALSE)]),
                        mean(method12[upper.tri(method12, diag = FALSE)]),
                        mean(method13[upper.tri(method13, diag = FALSE)])
                        
  ))
  
  colnames(df) <- c("Value")
  rownames(df) <- c("BayesSpace_DESpace",
                    "StLearn_DESpace",
                    "BayesSpace_findMarkers",
                    "StLearn_findMarkers",
                    "BayesSpace_FindAllMarkers",
                    "StLearn_FindAllMarkers",
                    "MERINGUE", "SPARK",
                    "nnSVG",
                    "SPARK-X", "SpatialDE2",
                    "SpatialDE","SpaGCN")
  df <- reshape::melt(t(df))
  df$X2 <- as.factor(df$X2)
  #all_colours = c(rep("#009E73",1), rep("#CC79A7",1), rep("#0072B2",1), rep("#D55E00",1), rep("#CD853F",1))
  gg_data <- df[,2:3]
  colnames(gg_data) <- c( "Method", "value")
  breaks_plots = c(0.05,0.1,0.15,0.2)
  gg_data$Method <- factor(gg_data$Method,levels=c("BayesSpace_DESpace",
                                                   "StLearn_DESpace",#"edgeR_CPMs_BayesSpace", 
                                                   "BayesSpace_findMarkers",
                                                   "StLearn_findMarkers",
                                                   "BayesSpace_FindAllMarkers",
                                                   "StLearn_FindAllMarkers",
                                                   "MERINGUE", "SPARK",
                                                   "nnSVG",
                                                   "SPARK-X", "SpatialDE2",
                                                   "SpatialDE","SpaGCN"))
  
  
  #pdf(file=paste0(path,save_name,data,".pdf"), onefile=FALSE) # or other device
  #gg_data %>%
  gg_Bar <- ggplot(gg_data,
                   #aes_string(y = "Time_Minutes", fill = "Method"),
                   aes(x = reorder(Method,-value), y = value, fill = Method)) + 
    
    #ggplot(aes_string(x = "Method", y = "value", fill = "Method")) +  
    geom_bar( stat = "identity",position = position_dodge(0.8), width = 0.7) +  theme_bw() +   xlab("") +  
    ylab("Average Jaccaed index across samples") + #facet_wrap(~Sample_id) +
    #scale_y_sqrt( breaks = c(0.1,0.15,0.25,0.4,0.5,0.6) ) + 
    #scale_y_continuous( breaks = c(0.05,0.1,0.15,0.2) ) + 
    scale_fill_manual("Method", values = col_method)+
    scale_y_continuous( breaks = breaks_plots ) +
    theme(axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1),   
          axis.text.y= element_text(size=rel(1.8),angle = 0, hjust = 1),  
          axis.title = element_text(size=15,face = "bold"),      
          panel.grid.minor = element_blank(),        
          panel.grid.major.x = element_blank(),       
          panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),        
          aspect.ratio = 1,        
          legend.position = "none",
          legend.text=element_text(size=15))
  return(gg_Bar)
}


########################### main simulation #############################
merge_results <- function(j,sample_names = sample_names, 
                               pattern_names = patterns[i],
                               path = path,new_data_path = new_data_path,
                               spatial_probs = c(0.5,0.9),default = FALSE,rep_i = 1){
  if(!is.null(j)){
    dir = paste0(path,sample_names[j],'/',pattern_names,'/probs_',spatial_probs[1],'_',
                 spatial_probs[2],default,'/')
    new_data_dir = paste0(new_data_path,sample_names[j],'/',pattern_names,'/probs_',spatial_probs[1],'_',
                          spatial_probs[2],default,'/')
  }else{
    dir = paste0(path,"SlideSeq2/",'/',pattern_names,'/probs_',spatial_probs[1],'_',
                 spatial_probs[2],default,'/')
    new_data_dir <- paste0(new_data_path,'/',pattern_names,'/probs_',spatial_probs[1],'_',
                           spatial_probs[2],default,'/')
  }
    
    
    file = paste0(dir,'results_all.csv')
    result <- read.csv(file,header = TRUE,sep =  '\t')
    
    ## Load additional results: 
    
    load(paste0(new_data_dir, "BayesSpace_scran_results.rda"))
    BayesSpace_scran_results <- Reduce(
      function(x, y) merge(x, y, by = 'genes', all = T),
      lapply(scranResults, function(x) { x$genes <- rownames(x); x }))
    # BayesSpace_scran_results <- do.call(cbind, scranResults)
    # WRONG: orders of genes of each element in the list are different
    rm(scranResults)
    ind1 <- which(sub("\\..*", "", colnames(BayesSpace_scran_results)) == "p")
    ind2 <- which(sub("\\..*", "", colnames(BayesSpace_scran_results)) == "FDR")
    
    #ind1 <- which(sub(".*?\\.", "", colnames(BayesSpace_scran_results)) == "p.value")
    #ind2 <- which(sub(".*?\\.", "", colnames(BayesSpace_scran_results)) == "FDR")
    ## Check whether types of p-val and FDR columns are numeric instead of character -> change the minimum value
    if(sum(sapply(BayesSpace_scran_results[ind1], class) != 'numeric') != 0){
      BayesSpace_scran_results[ind1] <- sapply(BayesSpace_scran_results[ind1], as.numeric)
      BayesSpace_scran_results[ind2] <- sapply(BayesSpace_scran_results[ind2], as.numeric)
    }
    
    BayesSpace_scran_pval <- apply(BayesSpace_scran_results[ind1], 1, min)
    BayesSpace_scran_FDR <- apply(BayesSpace_scran_results[ind2], 1, min)
    
    BayesSpace_scran_results <- data.frame(genes = BayesSpace_scran_results$genes,
                                           adj.p.value = BayesSpace_scran_FDR,
                                           p.value = BayesSpace_scran_pval,
                                           method = "BayesSpace_findMarkers")
    
    
    load(paste0(new_data_dir, "stLearn_scran_results.rda"))
    stLearn_scran_results <- Reduce(
      function(x, y) merge(x, y, by = 'genes', all = T),
      lapply(scranResults, function(x) { x$genes <- rownames(x); x }))
    rm(scranResults)
    
    ind1 <- which(sub("\\..*", "", colnames(stLearn_scran_results)) == "p")
    ind2 <- which(sub("\\..*", "", colnames(stLearn_scran_results)) == "FDR")
    ## Check whether types of p-val and FDR columns are numeric instead of character -> change the minimum value
    if(sum(sapply(stLearn_scran_results[ind1], class) != 'numeric') != 0){
      stLearn_scran_results[ind1] <- sapply(stLearn_scran_results[ind1], as.numeric)
      stLearn_scran_results[ind2] <- sapply(stLearn_scran_results[ind2], as.numeric)
    }
    stLearn_scran_pval <- apply(stLearn_scran_results[ind1], 1, min)
    stLearn_scran_FDR <- apply(stLearn_scran_results[ind2], 1, min)
    
    stLearn_scran_results <- data.frame(genes = stLearn_scran_results$genes,
                                        adj.p.value = stLearn_scran_FDR,
                                        p.value = stLearn_scran_pval,
                                        method = "StLearn_findMarkers")
    ## Seurat
    load(paste0(new_data_dir, "BayesSpace_seurat_results_rawcounts.rda"))
    BayesSpace_seurat_results <- as.data.frame( seuratResults)
    
    BayesSpace_seurat <- reshape(BayesSpace_seurat_results, idvar = "gene", 
                                 v.names = c("p_val", "p_val_adj"), 
                                 timevar = "cluster", direction = "wide")
    
    ind1 <- which(sub("\\..*", "", colnames(BayesSpace_seurat)) == "p_val")
    ind2 <- which(sub("\\..*", "", colnames(BayesSpace_seurat)) == "p_val_adj")
    ## Check whether types of p-val and FDR columns are numeric instead of character -> change the minimum value
    if(sum(sapply(BayesSpace_seurat[ind1], class) != 'numeric') != 0){
      BayesSpace_seurat[ind1] <- sapply(BayesSpace_seurat[ind1], as.numeric)
      BayesSpace_seurat[ind2] <- sapply(BayesSpace_seurat[ind2], as.numeric)
    }
    BayesSpace_seurat_pval <- apply(BayesSpace_seurat[ind1], 1, min)
    BayesSpace_seurat_FDR <- apply(BayesSpace_seurat[ind2], 1, min)
    
    BayesSpace_seurat_results <- data.frame(genes = BayesSpace_seurat$gene,
                                            adj.p.value = BayesSpace_seurat_FDR,
                                            p.value = BayesSpace_seurat_pval,
                                            method = "BayesSpace_FindAllMarkers")
    
    
    load(paste0(new_data_dir, "stLearn_seurat_results_rawcounts.rda"))
    stLearn_seurat_results <- as.data.frame( seuratResults)
    
    stLearn_seurat <- reshape(stLearn_seurat_results, idvar = "gene", 
                              v.names = c("p_val", "p_val_adj"), 
                              timevar = "cluster", direction = "wide")
    
    ind1 <- which(sub("\\..*", "", colnames(stLearn_seurat)) == "p_val")
    ind2 <- which(sub("\\..*", "", colnames(stLearn_seurat)) == "p_val_adj")
    ## Check whether types of p-val and FDR columns are numeric instead of character -> change the minimum value
    if(sum(sapply(stLearn_seurat[ind1], class) != 'numeric') != 0){
      stLearn_seurat[ind1] <- sapply(stLearn_seurat[ind1], as.numeric)
      stLearn_seurat[ind2] <- sapply(stLearn_seurat[ind2], as.numeric)
    }
    stLearn_seurat_pval <- apply(stLearn_seurat[ind1], 1, min)
    stLearn_seurat_FDR <- apply(stLearn_seurat[ind2], 1, min)
    
    stLearn_seurat_results <- data.frame(genes = stLearn_seurat$gene,
                                         adj.p.value = stLearn_seurat_FDR,
                                         p.value = stLearn_seurat_pval,
                                         method = "StLearn_FindAllMarkers")
    
    all_results <- rbind(result[, c("genes", "adj.p.value", "p.value", "method")],
                         BayesSpace_scran_results,
                         stLearn_scran_results,
                         BayesSpace_seurat_results,
                         stLearn_seurat_results)
    
    file = paste0(dir,'results_all_revisions.csv')
    write.csv(all_results, file, row.names = FALSE)
    return(all_results)
                               
}
roc_fdr_plot <- function(j,
    sample_names = sample_names, 
    pattern_names = patterns[i],
    path = path,colours = colours,new_data_path = new_data_path,
    spatial_probs = c(0.5,0.9),default = FALSE,rep_i = 1,
    methods_all = methods_all,
    methods_order = methods_order,
    shape_border = shape_border,
    shape_fill = shape_fill){
  
    dir = paste0(path,sample_names[j],'/',pattern_names,'/probs_',spatial_probs[1],'_',
                 spatial_probs[2],default,'/')
    file = paste0(dir,'results_all_revisions.csv')
    result <- read.csv(file,header = TRUE,sep =  ',')
    genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
    all_genes <- as.data.frame((unique(result$genes)))
    
    all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
    genes_id <- apply(genes,1,function(x) substr(x, 10,18))
    status <- ifelse(all_genes_id %in% genes_id, 1, 0)
    truth <- cbind(all_genes, status)
    # check
    sum(truth$status) == dim(genes)[1]
    if(sum(all_genes_id[length(all_genes_id)] == "_count")>=1){
      truth <- truth[1:dim(truth)[1]-1,]
    }
    result <- result %>% dplyr::distinct()
    setDT(result)
    result$p.value <- as.numeric(result$p.value)
    result$adj.p.value <- as.numeric(result$adj.p.value)
    data <- dcast(result, method ~ genes,value.var = 'adj.p.value', fun.aggregate=mean)
    #data <- dcast(result, method ~ genes,value.var = 'adj.p.value')
    df2 = truth[match(colnames(data)[-1], truth$`(unique(result$genes))`),]
    df2 <- as.data.frame(df2[,2])
    rownames(df2) <- colnames(data)[-1]
    colnames(df2) <- 'status'
    data2 <- dcast(result, method ~ genes,value.var = 'p.value', fun.aggregate=mean)
    data <- data %>%
      mutate(method =  factor(method, levels = methods_order)) %>%
      arrange(method)
    data2 <- data2 %>%
      mutate(method =  factor(method, levels = methods_order)) %>%
      arrange(method)
    
    pval = data.frame(t(data2[, -1])
    )
    colnames(pval) <- methods_all
    padj = data.frame(t(data[, -1])
    )
    colnames(padj) <- methods_all
    
    
    DF_COBRA <- COBRAData(pval,
                          padj,
                          truth = data.frame(df2))
    perf <- calculate_performance(DF_COBRA, binary_truth = "status",
                                  thrs = c(0.01, 0.05, 0.1, 0.2))
    cobra_plot <- prepare_data_for_plot(perf, colorscheme = colours, incloverall = FALSE,
                                        facetted = TRUE,conditionalfill = FALSE,
                                        keepmethods = c("BayesSpace_DESpace","StLearn_DESpace","SPARK",
                                                        "SPARK-X","SpatialDE","SpatialDE2","MERINGUE",
                                                        "SpaGCN","nnSVG",#"nnSVG_clusters","SPARK-X_clusters",
                                                        "BayesSpace_findMarkers","StLearn_findMarkers",
                                                        "BayesSpace_FindAllMarkers","StLearn_FindAllMarkers"
                                                        #,"limma_logCPM","limma_voomCPM", "limma_logcounts","DESeq2"
                                        ))
    my.cols <- colours
    (gg_roc_with_pval <- plot_roc(cobra_plot,linewidth=2) + 
        scale_color_manual(values = my.cols,
                           name = "",
                           breaks=methods_all,
                           labels=methods_all) +
        #scale_fill_brewer(breaks=methods, type="qual", palette=3)+
        # scale_color_manual(values = cobra_plot@plotcolors, name = "", labels = methods_order,limits = force) +
        theme(axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(6)),
              axis.text.y=element_text(size=rel(6)),
              axis.title.y = element_text(size=rel(6)),
              axis.title.x = element_text(size=rel(6)),
              legend.title=element_text(size=rel(2)),
              legend.text=element_text(size=rel(6)),
              legend.key.width=unit(2, "cm"),
              legend.position="bottom",
              legend.box="vertical", legend.margin=margin())) 
    ####
    if(spatial_probs[1] == 0.6){gg_roc_with_pval <- gg_roc_with_pval + labs(y="TPR_strong", x = "TPR_weak")}
    # ggsave(paste0(path_save,dataset,'/AggregatedSamples_',spatial_probs[1],'',spatial_probs[2],roc_save_name),
    # width = 3,height = 4
    # )
    ## Manually match "method" and point "shape"
    sel_shape = c(2,13,11,8,10,9,4,5,6,7,3,14,12)
    col = my.cols
    names(col) = methods_all
    #col = col[c()]
    # plot FDR/TPR curve
    (gg_fdr_with_pval <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                                         pointsize = 0, linewidth = 2)+
      scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.8)) +
      scale_color_manual(values = my.cols,
                         name = "",
                         breaks=methods_all,
                         labels=methods_all)+
      guides(colour = guide_legend(ncol = 2, byrow = FALSE,
                                   override.aes = list(shape = shape_fill[-1],
                                                       fill = col[-length(col)]) ) ) +
      geom_point(size = 10, aes(fill = method, colour = method, shape = method
      ),
      shape = rep(shape_border[sel_shape],4), 
      stroke = 2, 
      alpha = 1) + # stroke = line width
      scale_fill_manual(values =rep(col,4), 
                        name = "",
                        breaks=methods_all,
                        labels=methods_all)+
      geom_point(size = 10, aes(fill = method, colour = method, shape = method
      ),
      shape = rep(shape_fill[sel_shape],4), 
      stroke = 2, alpha = 0.25)+
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(6)),
            axis.text.y=element_text(size=rel(6)),
            axis.title.y = element_text(size=rel(6)),
            axis.title.x = element_text(size=rel(6)),
            legend.title=element_text(size=rel(2)),
            legend.text=element_text(size=rel(6)),
            legend.key.width=unit(2, "cm"),
            aspect.ratio = 1, legend.position="bottom",
            legend.box="vertical", legend.margin=margin()) )
    print(paste0("Sample: ", sample_names[j], " Pattern: ", pattern_names))
    return(list(gg_roc_with_pval,gg_fdr_with_pval ))
}

simple_auc <- function(TPR, FPR){
  # inputs already sorted, best scores first 
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}
# https://www.r-bloggers.com/2016/11/calculating-auc-the-area-under-a-roc-curve/
  
overall_roc_fdr_plot <- function(
  sample_names = sample_names, 
  pattern_names = patterns[i],
  path = path,colours = colours,new_data_path = new_data_path,
  spatial_probs = c(0.5,0.9),default = FALSE,rep_i = 1,
  methods_all = methods_all,
  methods_order = methods_order,
  shape_border = shape_border,
  shape_fill = shape_fill,
  plot = TRUE,
  auc = FALSE){
 
    combined_data <- data.frame(0)
    combined_data2 <- data.frame(0)
    combined_df2 <- NULL
    combined_result <- NULL
    combined_genes <- NULL
    combined_all_genes <- NULL
    
    for(kk in c(1:length(sample_names))){
      dir = paste0(path,sample_names[kk],'/',pattern_names,'/probs_',spatial_probs[1],'_',
                   spatial_probs[2],default,'/')
      file = paste0(dir,'results_all_revisions.csv')
      result <- read.csv(file,header = TRUE,sep =  ',')
      genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
      all_genes <- as.data.frame((unique(result$genes)))
      
      all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
      genes_id <- apply(genes,1,function(x) substr(x, 10,18))
      status <- ifelse(all_genes_id %in% genes_id, 1, 0)
      truth <- cbind(all_genes, status)
      # check
      sum(truth$status) == dim(genes)[1]
      if(sum(all_genes_id[length(all_genes_id)] == "_count")>=1){
        truth <- truth[1:dim(truth)[1]-1,]
      }
      result <- result %>% dplyr::distinct()
      setDT(result)
      
      data <- dcast(result, method ~ genes,value.var = 'adj.p.value')
      #data <- dcast(result, method ~ genes,value.var = 'adj.p.value')
      df2 = truth[match(colnames(data)[-1], truth$`(unique(result$genes))`),]
      df2 <- as.data.frame(df2[,2])
      rownames(df2) <- colnames(data)[-1]
      colnames(df2) <- 'status'
      data2 <- dcast(result, method ~ genes,value.var = 'p.value')
      
      colnames(data) <- paste0(colnames(data), ".",kk)
      colnames(data2) <- paste0(colnames(data2),".", kk)
      rownames(df2) <- paste0(rownames(df2), ".",kk)
      
      colnames(data)[1] <- "method"
      colnames(data2)[1] <- "method"
      
      combined_data <- cbind(combined_data,data)
      combined_data2 <- cbind(combined_data2,data2)
      combined_df2 <- rbind(combined_df2,df2)
    }
    
    combined_data <- combined_data[,-1]
    combined_data2 <- combined_data2[,-1]
    
    methods <- combined_data[,1]
    combined_data <- combined_data %>% as.data.frame() %>% dplyr::select(-contains("method"))
    combined_data <- t(combined_data)
    colnames(combined_data) <- methods
    
    methods <- combined_data2[,1]
    combined_data2 <- combined_data2 %>% as.data.frame() %>% dplyr::select(-contains("method"))
    combined_data2 <- t(combined_data2)
    colnames(combined_data2) <- methods
    
    genes_id <- rownames(combined_data2)
    combined_data2 <- apply(combined_data2, 2,            # Specify own function within apply
                            function(x) as.numeric(as.character(x)))
    rownames(combined_data2) <- genes_id
    
    combined_data <- apply(combined_data, 2,            # Specify own function within apply
                           function(x) as.numeric(as.character(x)))
    rownames(combined_data) <- genes_id
    
    ind <- match(methods_order, colnames(combined_data))
    combined_data <-  combined_data %>% as.data.frame() %>% dplyr::select(all_of(ind))
    ind <- match(methods_order, colnames(combined_data2))
    combined_data2 <-  combined_data2 %>% as.data.frame() %>% dplyr::select(all_of(ind))
    
    pval = data.frame(combined_data2)
    colnames(pval) <- methods_all
    padj = data.frame(combined_data)
    colnames(padj) <- methods_all
    
    if(auc == TRUE){
      # DF_pval <- sapply(pval, function(x) {
      #   r = roc(combined_df2$status, x)
      #   auc(r)
      # })
      
      ## same as: 
      DF_COBRA <- COBRAData(pval,
                            padj,
                            truth = data.frame(combined_df2))
      perf <- calculate_performance(DF_COBRA, binary_truth = "status",aspects = "roc")
      roc_pval <- iCOBRA::roc(perf)
      DF_pval <- sapply(methods_all, function(x) {
          tmp_roc <- roc_pval[roc_pval$method == x,]
          with(tmp_roc, simple_auc(TPR, FPR))
      })

      # DF_padj <- sapply(padj, function(x) {
      #   r = roc(combined_df2$status, x)
      #   auc(r)
      # })
      return(list(DF_pval))
    }
    
    if(plot == TRUE){
      DF_COBRA <- COBRAData(pval,
                            padj,
                            truth = data.frame(combined_df2))
        perf <- calculate_performance(DF_COBRA, binary_truth = "status",
                                      thrs = c(0.01, 0.05, 0.1, 0.2))
        cobra_plot <- prepare_data_for_plot(perf, colorscheme = colours, incloverall = FALSE,
                                            facetted = TRUE,conditionalfill = FALSE,
                                            keepmethods = methods_all)
        my.cols <- colours
        (gg_roc_with_pval <- plot_roc(cobra_plot,linewidth=1.5) +
            scale_color_manual(values = my.cols,
                               name = "",
                               breaks=methods_all,
                               labels=methods_all) +
            #scale_fill_brewer(breaks=methods, type="qual", palette=3)+
            # scale_color_manual(values = cobra_plot@plotcolors, name = "", labels = methods_order,limits = force) +
            theme(axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(6)),
                  axis.text.y=element_text(size=rel(6)),
                  axis.title.y = element_text(size=rel(6)),
                  axis.title.x = element_text(size=rel(6)),
                  legend.title=element_text(size=rel(2)),
                  legend.text=element_text(size=rel(6)),
                  legend.key.width=unit(2, "cm"),
                  legend.position="bottom",
                  legend.box="vertical", legend.margin=margin())) 
        ####
        if(spatial_probs[1] == 0.6){gg_roc_with_pval <- gg_roc_with_pval + labs(y="TPR_strong", x = "TPR_weak")}
        # ggsave(paste0(path_save,dataset,'/AggregatedSamples_',spatial_probs[1],'',spatial_probs[2],roc_save_name),
        # width = 3,height = 4
        # )
        ## Manually match "method" and point "shape"
        sel_shape = c(1,12,10,7,9,8,3,4,5,6,2,13,11)+1
        col = my.cols[c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
        # plot FDR/TPR curve
        gg_fdr_with_pval <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                                             pointsize = 0, linewidth = 2)+
          scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.8)) +
          scale_color_manual(values = col,
                             name = "",
                             breaks=methods_all,
                             labels=methods_all) +
          guides(colour = guide_legend(ncol = 2, byrow = FALSE,
                                       override.aes = list(shape = shape_fill[-1],
                                                           fill = col) ) ) +
          geom_point(size = 10, aes(fill = method, colour = method, shape = method
          ),
          shape = rep(shape_border[sel_shape],4), 
          stroke = 2, 
          alpha = 1) + # stroke = line width
          geom_point(size = 10, aes(fill = method, colour = method, shape = method
          ),
          shape = rep(shape_fill[sel_shape],4), 
          stroke = 2, alpha = 0.25)+
          #scale_fill_manual(values = c(rep(c(my.cols[c(1,2,3,4,5,6)],rep("white",9)),2),"white","white","white"), guide = "none")+
          scale_fill_manual(values = rep(col,4),
                            name = "",
                            breaks=methods_all,
                            labels=methods_all)+
          theme(strip.background = element_blank(),
                strip.text = element_blank(),
                axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(6)),
                axis.text.y=element_text(size=rel(6)),
                axis.title.y = element_text(size=rel(6)),
                axis.title.x = element_text(size=rel(6)),
                legend.title=element_text(size=rel(2)),
                legend.text=element_text(size=rel(6)),
                legend.key.width=unit(2, "cm"),
                aspect.ratio = 1, legend.position="bottom",
                legend.box="vertical", legend.margin=margin())
        print(paste0("Sample: ", sample_names[kk], " Pattern: ", pattern_names))
        return(list(gg_roc_with_pval,gg_fdr_with_pval ))
    }
}

overall_roc_fdr_plot_cerebellum <- function(
    pattern_names = patterns[i],
    path = path,colours = colours,new_data_path = new_data_path,
    spatial_probs = c(0.5,0.9),default = FALSE,rep_i = 1,
    methods_all = methods_all,
    methods_order = methods_order,
    shape_border = shape_border,
    shape_fill = shape_fill){
  
  dir = paste0(path,"SlideSeq2/",'/',pattern_names,'/probs_',spatial_probs[1],'_',
               spatial_probs[2],default,'/')
  file = paste0(dir,'results_all_revisions.csv')
  result <- read.csv(file,header = TRUE,sep =  ',')
  genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
  all_genes <- as.data.frame((unique(result$genes)))
  
  all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
  if(pattern_names %in% c("mixture_patch", "mixture_reverse_patch")){
    status <- ifelse(all_genes$`(unique(result$genes))` %in% genes$V1, 1, 0)
  }else{
    status <- ifelse(all_genes$`(unique(result$genes))` %in% genes$gene_ids, 1, 0)
  }
  truth <- cbind(all_genes, status)
  # check
  sum(truth$status) == dim(genes)[1]
  if(sum(all_genes_id[length(all_genes_id)] == "_count")>=1){
    truth <- truth[1:dim(truth)[1]-1,]
  }
  result <- result %>% dplyr::distinct()
  setDT(result)
  result <- result[result$method %in% methods_order,]
  result$p.value <- as.numeric(result$p.value)
  result$adj.p.value <- as.numeric(result$adj.p.value)
  data <- dcast(result, method ~ genes,value.var = 'adj.p.value', fun.aggregate=mean)
  #data <- dcast(result, method ~ genes,value.var = 'adj.p.value')
  df2 = truth[match(colnames(data)[-1], truth$`(unique(result$genes))`),]
  df2 <- as.data.frame(df2[,2])
  rownames(df2) <- colnames(data)[-1]
  colnames(df2) <- 'status'
  data2 <- dcast(result, method ~ genes,value.var = 'p.value', fun.aggregate=mean)
  data <- data %>%
    mutate(method =  factor(method, levels = methods_order)) %>%
    arrange(method)
  data2 <- data2 %>%
    mutate(method =  factor(method, levels = methods_order)) %>%
    arrange(method)
  
  pval = data.frame(t(data2[, -1])
  )
  colnames(pval) <- methods_all
  padj = data.frame(t(data[, -1])
  )
  colnames(padj) <- methods_all
  ## SpaGCN do not provide p-values. Here SpaGCN column in pval is actually `1-in_out_group_ratio`.
  ## It contains negative values. So we replace SpaGCN column in pval by SpaGCN column in padj.
  if(!is.null(pval$SpaGCN)){pval$SpaGCN = padj$SpaGCN}
  DF_COBRA <- COBRAData(pval,
                        padj,
                        truth = data.frame(df2))
  perf <- calculate_performance(DF_COBRA, binary_truth = "status",
                                thrs = c(0.01, 0.05, 0.1, 0.2))
  cobra_plot <- prepare_data_for_plot(perf, colorscheme = colours, incloverall = FALSE,
                                      facetted = TRUE,conditionalfill = FALSE,
                                      keepmethods = methods_all)
  my.cols <- colours
  (gg_roc_with_pval <- plot_roc(cobra_plot,linewidth=2) + 
      scale_color_manual(values = my.cols,
                         name = "",
                         breaks=methods_all,
                         labels=methods_all) +
      #scale_fill_brewer(breaks=methods, type="qual", palette=3)+
      # scale_color_manual(values = cobra_plot@plotcolors, name = "", labels = methods_order,limits = force) +
      theme(axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(6)),
            axis.text.y=element_text(size=rel(6)),
            axis.title.y = element_text(size=rel(6)),
            axis.title.x = element_text(size=rel(6)),
            legend.title=element_text(size=rel(2)),
            legend.text=element_text(size=rel(6)),
            legend.key.width=unit(2, "cm"),
            legend.position="bottom",
            legend.box="vertical", legend.margin=margin())) 
  ####
  if(spatial_probs[1] == 0.6){gg_roc_with_pval <- gg_roc_with_pval + labs(y="TPR_strong", x = "TPR_weak")}
  # ggsave(paste0(path_save,dataset,'/AggregatedSamples_',spatial_probs[1],'',spatial_probs[2],roc_save_name),
  # width = 3,height = 4
  # )
  ## Manually match "method" and point "shape"
  if(pattern_names %in% c("mixture_patch", "mixture_reverse_patch")){
    sel_shape = c(2,13,11,8,10,4,5,6,7,3,14,12)
    col = my.cols[1:12]
    shape = shape_fill[-c(1,9)]
    #col_value = c(rep(c(my.cols[c(1,2,3,4,5,6)],rep("white",8)),2),"white","white","white")
  }else{
    sel_shape = c(2,13,11,8,10,9,4,5,6,7,3,14,12)
    col = my.cols[1:13]
    shape = shape_fill[-c(1)]
    #col_value = c(rep(c(my.cols[c(1,2,3,4,5,6)],rep("white",8)),2),"white","white","white")
  }
  names(col) = methods_all
  #col = col[c()]
  # plot FDR/TPR curve
  gg_fdr_with_pval <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                                       pointsize = 0, linewidth = 2)+
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.8)) +
    scale_color_manual(values = col,
                       name = "",
                       breaks=methods_all,
                       labels=methods_all) +
    guides(colour = guide_legend(ncol = 2, byrow = FALSE,
                                 override.aes = list(shape = shape,
                                                     fill = col) ) ) +
    geom_point(size = 10, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_border[sel_shape],4), 
    stroke = 2, 
    alpha = 1) + # stroke = line width
  geom_point(size = 10, aes(fill = method, colour = method, shape = method
  ),
  shape = rep(shape_fill[sel_shape],4), 
  stroke = 2, alpha = 0.25)+
    #scale_fill_manual(values =col_value, guide = "none")+
    scale_fill_manual(values = rep(col,4),
                      name = "",
                      breaks=methods_all,
                      labels=methods_all)+
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(6)),
          axis.text.y=element_text(size=rel(6)),
          axis.title.y = element_text(size=rel(6)),
          axis.title.x = element_text(size=rel(6)),
          legend.title=element_text(size=rel(2)),
          legend.text=element_text(size=rel(6)),
          legend.key.width=unit(2, "cm"),
          aspect.ratio = 1, legend.position="bottom",
          legend.box="vertical", legend.margin=margin())
  print(paste0(" Pattern: ", pattern_names))

  return(list(gg_roc_with_pval, gg_fdr_with_pval))
}


####################### Multiple samples simulation #######################
fct_spatialLIBD <- function(ll){
  ind1 <- which(sub("\\_.*", "", colnames(ll)) == "p")
  ind2 <- which(sub("\\_.*", "", colnames(ll)) == "fdr")
  ## Check whether types of p-val and FDR columns are numeric instead of character -> change the minimum value
  if(sum(sapply(ll[ind1], class) != 'numeric') != 0){
    ll[ind1] <- sapply(ll[ind1], as.numeric)
    ll[ind2] <- sapply(ll[ind2], as.numeric)
  }
  enrichment_pval <- apply(ll[ind1], 1, min)
  enrichment_FDR <- apply(ll[ind2], 1, min)
  
  data.frame(genes = ll$gene,
             FDR = enrichment_FDR,
             PValue = enrichment_pval,
             methods = "spatialLIBD")
}

overall_roc_fdr_plot_multi_samples <- function(
    path = "~",
    spatial_probs = c(0.5,0.8),
    # roc_save_name='_ROC.pdf',
    # fdr_save_name = '_FDR.pdf',
    all_colours = all_colours,
    shape_border = c(0, 1),
    shape_fill = c(15, 16),
    dataset = "LIBD",
    pattern_names#,
    #sample_names
){
  colours = c(all_colours[c(3,12,14)],"burlywood3", "darkgoldenrod1", "yellow", "white")
  # points, borders:
  shape_border = c(0,2,5,6,3,1)
  # points, fill:
  shape_fill = c(15,17,23,25,3,1)
  
  path_data = paste0(path, dataset, pattern_names, "/probs_", spatial_probs[1], "_", spatial_probs[2], "FALSE/")
  
  load(paste0(path_data,"result_SV_edgeR_counts.rda"))
  result1 = spatial_gene_results; rm(spatial_gene_results)
  load(paste0(path_data,"result_SV_scran_findMarkers.rda"))
  result2 = spatial_gene_results; rm(spatial_gene_results)
  load(paste0(path_data,"result_SV_seurat_FindAllMarkers.rda"))
  result3 = spatial_gene_results; rm(spatial_gene_results)
  load(paste0(path_data,"result_SV_spatialLIBD_registration.rda"))
  result4 = spatial_gene_results; rm(spatial_gene_results)
  
  ## scran:findMarkers -> return a list of results; names(list) = cluster names
  multi_scran_results <- Reduce(
    function(x, y) merge(x, y, by = 'genes', all = T),
    lapply(result2, function(x) { x$genes <- rownames(x); x }))
  
  ind1 <- which(sub("\\..*", "", colnames(multi_scran_results)) == "p")
  ind2 <- which(sub("\\..*", "", colnames(multi_scran_results)) == "FDR")
  ## Check whether types of p-val and FDR columns are numeric instead of character -> change the minimum value
  if(sum(sapply(multi_scran_results[ind1], class) != 'numeric') != 0){
    multi_scran_results[ind1] <- sapply(multi_scran_results[ind1], as.numeric)
    multi_scran_results[ind2] <- sapply(multi_scran_results[ind2], as.numeric)
  }
  multi_scran_pval <- apply(multi_scran_results[ind1], 1, min)
  multi_scran_FDR <- apply(multi_scran_results[ind2], 1, min)
  
  multi_scran_results <- data.frame(genes = multi_scran_results$genes,
                                    FDR = multi_scran_FDR,
                                    PValue = multi_scran_pval,
                                    methods = "findMarkers")
  
  ## seurat::FindAllMarkers -> return a dataframe; multiple p values for each gene (depends on the number of clusters) 
  multi_seurat <- reshape(result3[[1]], idvar = "gene", 
                          v.names = c("p_val", "p_val_adj"), 
                          timevar = "cluster", direction = "wide")
  
  ind1 <- which(sub("\\..*", "", colnames(multi_seurat)) == "p_val")
  ind2 <- which(sub("\\..*", "", colnames(multi_seurat)) == "p_val_adj")
  if(sum(sapply(multi_seurat[ind1], class) != 'numeric') != 0){
    multi_seurat[ind1] <- sapply(multi_seurat[ind1], as.numeric)
    multi_seurat[ind2] <- sapply(multi_seurat[ind2], as.numeric)
  }
  multi_seurat_pval <- apply(multi_seurat[ind1], 1, min)
  multi_seurat_FDR <- apply(multi_seurat[ind2], 1, min)
  
  multi_seurat_results <- data.frame(genes = multi_seurat$gene,
                                     FDR = multi_seurat_FDR,
                                     PValue = multi_seurat_pval,
                                     methods = "FindAllMarkers")
  
  ## spatialLIBD registration: return a list; the length of list = 4 or 6 (only mixture and inverted mixture patterns) 
  # results_anova, results_anova_nan,
  # results_enrichment, results_enrichment_nan,
  # results_pairwise, results_pairwise_nan
  results_registration <- lapply(result4, fct_spatialLIBD)
  
  if(length(result4) == 4){
    ind_enrich <- 1;ind_pair <- 3
  }else{
    ind_anova <- 1;ind_enrich <- 3;ind_pair <- 5
    results_anova <- results_registration[[ind_anova]]
    results_anova$methods <- "spatialLIBD_anova"
  }
  
  ## spatialLIBD::registration_stats_enrichment
  results_enrichment <- results_registration[[ind_enrich]]
  results_enrichment$methods <- "spatialLIBD_enrichment"
  
  ## spatialLIBD::registration_pairwise
  results_pairwise <- results_registration[[ind_pair]]
  results_pairwise$methods <- "spatialLIBD_pairwise"
  
  result <- rbind(multi_scran_results, multi_seurat_results,
                  results_enrichment, results_pairwise)
  methods_order = c(#'Single_Sample_DESpace',
    'DESpace',
    "FindAllMarkers","findMarkers",
    "spatialLIBD_enrichment","spatialLIBD_pairwise"
    #'stLearn_single_DESpace','stLearn_multi_DESpace'
  )
  if(pattern_names %in% c("MixClusters_pattern","MixClusters_reverse_pattern")){
    result <- rbind(result,
                    results_anova)
    methods_order <- c(methods_order,"spatialLIBD_anova")
    ind_shape <- c(1,2,3,6,4,5)
  }else{
    colours <- colours[-6]
    shape_border <- shape_border[-6]
    shape_fill <- shape_fill[-6]
    ind_shape <- c(1:5)
  }
  
  result$sample <- ""
  file = paste0(path_data, 'probs_',spatial_probs[2], '_selected_genes.txt')
  genes <- read.csv(file,header = TRUE,sep =  '\t')
  all_genes <- as.data.frame((unique(result$genes)))
  
  all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
  genes_id <- apply(genes,1,function(x) substr(x, 10,18))
  status <- ifelse(all_genes_id %in% genes_id, 1, 0)
  truth <- cbind(all_genes, status)
  # check
  sum(truth$status) == dim(genes)[1]
  
  result <- result %>% dplyr::distinct()
  # setDT(result)
  ################################
  single_Bayes = as.data.table(result1 %>% filter(method == "SV_single_edgeR" & cluster == "Original"))
  multi_Bayes = as.data.table(result1 %>% filter(method == "SV_multi_edgeR" & cluster == "Original"))
  single_Bayes[, methods := 'Single_Sample_DESpace']
  multi_Bayes[, methods := 'DESpace']
  
  single_Bayes = single_Bayes[,.(genes, FDR, PValue, sample, methods)]
  multi_Bayes = multi_Bayes[,.(genes, FDR, PValue, sample, methods)]
  
  samples = unique(single_Bayes$sample)
  combined_data <- data.frame(0)
  combined_data2 <- data.frame(0)
  combined_df2 <- NULL
  combined_result <- NULL
  combined_genes <- NULL
  combined_all_genes <- NULL
  
  for(ii in c(1:length(samples))){
    df1 <- rbind(
      
      #single_Bayes %>% filter(sample == samples[ii]),
      multi_Bayes
    )
    combine_results = rbind(setDT(df1), setDT(result), fill=TRUE)
    
    data <- dcast(combine_results, methods ~ genes,value.var = 'FDR')
    df2 = truth[match(colnames(data)[-1], truth$`(unique(result$genes))`),]
    df2 <- as.data.frame(df2[,2])
    rownames(df2) <- colnames(data)[-1]
    colnames(df2) <- 'status'
    data2 <- dcast(combine_results, methods ~ genes,value.var = 'PValue')
    
    colnames(data) <- paste0(colnames(data), ".",ii)
    colnames(data2) <- paste0(colnames(data2),".", ii)
    rownames(df2) <- paste0(rownames(df2), ".",ii)
    
    colnames(data)[1] <- "method"
    colnames(data2)[1] <- "method"
    
    data <- data[match(methods_order, as.vector(data[,1])$method),]
    data2 <- data2[match(methods_order, as.vector(data2[,1])$method),]
    #print(data[,1])
    data <- data[,-1]
    data2 <- data2[,-1]
    
    #data <- data[,1:(dim(data)[2]-1)]
    #data2 <- data2[,1:(dim(data2)[2]-1)]
    combined_data <- cbind(combined_data,data)
    combined_data2 <- cbind(combined_data2,data2)
    combined_df2 <- rbind(combined_df2,df2)
  }
  
  combined_data <- combined_data[,-1]
  combined_data2 <- combined_data2[,-1]
  
  combined_data <- t(combined_data)
  colnames(combined_data) <- methods_order
  combined_data2 <- t(combined_data2)
  colnames(combined_data2) <- methods_order
  
  genes_id <- rownames(combined_data2)
  combined_data2 <- apply(combined_data2, 2,            # Specify own function within apply
                          function(x) as.numeric(as.character(x)))
  rownames(combined_data2) <- genes_id
  
  combined_data <- apply(combined_data, 2,            # Specify own function within apply
                         function(x) as.numeric(as.character(x)))
  rownames(combined_data) <- genes_id
  combined_data <- combined_data[,match(methods_order, colnames(combined_data))]
  combined_data2 <- combined_data2[,match(methods_order, colnames(combined_data2))]
  
  DF_COBRA <- COBRAData(pval = data.frame(combined_data2
  ),
  padj = data.frame(combined_data
  ),
  truth = data.frame(combined_df2))
  # truth = 1 for SV genes and 0 for uniform genes.
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # ROC, FDR plots:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  perf <- calculate_performance(DF_COBRA, binary_truth = "status",
                                thrs = c(0.01, 0.05, 0.1, 0.2))
  cobra_plot <- prepare_data_for_plot(perf, colorscheme = colours,incloverall = FALSE,
                                      facetted = TRUE,conditionalfill = FALSE)
  cobra_plot@plotcolors <- cobra_plot@plotcolors[1:(4+length(result4)/2)]
  
  (gg_roc <- plot_roc(cobra_plot,linewidth=2)  + 
      #scale_fill_brewer(breaks=methods, type="qual", palette=3)+
      # scale_color_manual(values = cobra_plot@plotcolors, name = "", labels = methods_order,limits = force) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                       hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            legend.position = "top") )
  gg_roc <- gg_roc +
    guides(color=guide_legend(nrow=1,byrow=TRUE))
  
  #ggsave(paste0(path,roc_save_name),width = 7, height = 5)
  col = colours[-length(colours)]
  names(col) = methods_order
  # plot FDR/TPR curve
  gg_fdr <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                             pointsize = 0, linewidth = 2)+ 
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.8)) +
    scale_color_manual(values = col,
                       name = "",
                       breaks=methods_order,
                       labels=methods_order)+
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(3)),
          axis.text.y=element_text(size=rel(3)),
          axis.title.y = element_text(size=rel(3)),
          axis.title.x = element_text(size=rel(3)),
          legend.title=element_text(size=rel(2)),
          legend.text=element_text(size=rel(3)),
          legend.key.width=unit(2, "cm"),
          aspect.ratio = 1, legend.position="bottom",
          legend.box="vertical", legend.margin=margin())  +
    guides(colour = guide_legend(ncol = 2, byrow = FALSE,
                                 override.aes = list(shape = shape_fill,
                                                     fill = colours[-length(colours)]) ) ) +
    geom_point(size = 10, aes(fill = method, colour =method, shape = method), 
               shape = rep(shape_border[ind_shape],4), stroke = 2, alpha = 1) + # stroke = line width
    scale_fill_manual(values =col, guide = "none") +
    geom_point(size = 10, aes(fill = method, colour = method, shape = method), 
               shape = rep(shape_fill[ind_shape],4), stroke = 2, alpha = 0.25)
  return(list(gg_roc, gg_fdr))
  
}


roc_fdr_plot_multi_samples_DESpace <- function(
    i,
    sample_names,
    path = "~",
    spatial_probs = c(0.5,0.8),
    # roc_save_name='_ROC.pdf',
    # fdr_save_name = '_FDR.pdf',
    colours = colours,
    shape_border = c(0, 1),
    shape_fill = c(15, 16),
    dataset = "LIBD",
    pattern_names#,
    #sample_names
){
  path_data = paste0(path, dataset, "/", pattern_names, "/probs_", spatial_probs[1], "_", spatial_probs[2], "FALSE/")
  # for(q in seq_along(sample_names)) { 
  #  path_sub = paste0(path_data, "probs_", spatial_probs[1], "_", spatial_probs[2], "_final_object",sample_names[q], '.rda')
  #  load(path_sub)
  #  if(!exists("combined_object")){combined_object = final_object; remove(final_object)}
  #  assign(paste0("sce", q), combined_object)
  # }
  #  gene1 <- rownames(sce1);gene2 <- rownames(sce2);gene3 <- rownames(sce3)
  #  sum(gene1 == gene2)== 5000; sum(gene1 == gene3)== 5000
  
  load(paste0(path_data,"result_SV_edgeR_counts.rda"))
  result = result1 = spatial_gene_results
  file = paste0(path_data, 'probs_',spatial_probs[2], '_selected_genes.txt')
  genes <- read.csv(file,header = TRUE,sep =  '\t')
  all_genes <- as.data.frame((unique(result$genes)))
  
  all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
  genes_id <- apply(genes,1,function(x) substr(x, 10,18))
  status <- ifelse(all_genes_id %in% genes_id, 1, 0)
  truth <- cbind(all_genes, status)
  # check
  sum(truth$status) == dim(genes)[1]
  
  result <- result %>% dplyr::distinct()
  setDT(result)
  ################################
  methods_order = c('Single_Sample_DESpace','Multiple_Samples_DESpace'#,
                    #'stLearn_single_DESpace','stLearn_multi_DESpace'
  )
  single_Bayes = as.data.table(result1 %>% filter(method == "SV_single_edgeR" & cluster == "Original"))
  multi_Bayes = as.data.table(result1 %>% filter(method == "SV_multi_edgeR" & cluster == "Original"))
  single_Bayes[, method := 'Single_Sample_DESpace']
  multi_Bayes[, method := 'Multiple_Samples_DESpace']
  
  single_Bayes = single_Bayes[,.(genes, FDR, PValue, sample, method)]
  multi_Bayes = multi_Bayes[,.(genes, FDR, PValue, sample, method)]
  
  combine_results = rbind(
    single_Bayes %>% filter(sample == sample_names[i]),
    multi_Bayes
  )
  setDT(combine_results)
  data <- dcast(combine_results, method ~ genes,value.var = 'FDR')
  df2 = truth[match(colnames(data)[-1], truth$"(unique(result$genes))"),]
  df2 <- as.data.frame(df2[,2])
  rownames(df2) <- colnames(data)[-1]
  colnames(df2) <- 'status'
  data2 <- dcast(combine_results, method ~ genes,value.var = 'PValue')
  
  data <- data %>%
    mutate(method = factor(method, levels = methods_order)) %>%
    arrange(method)
  data2 <- data2 %>%
    mutate(method = factor(method, levels = methods_order)) %>%
    arrange(method)
  pval = data.frame(t(data2[, -1])
  )
  colnames(pval) <- methods_order
  padj = data.frame(t(data[, -1])
  )
  colnames(padj) <- methods_order
  
  
  DF_COBRA <- COBRAData(pval = data.frame(pval
  ),
  padj = data.frame(padj
  ),
  truth = data.frame(df2))
  # truth = 1 for SV genes and 0 for uniform genes.
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # ROC, FDR plots:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  perf <- calculate_performance(DF_COBRA, binary_truth = "status",
                                thrs = c(0.01, 0.05, 0.1, 0.2))
  cobra_plot <- prepare_data_for_plot(perf, colorscheme = colours,incloverall = FALSE,
                                      facetted = TRUE,conditionalfill = FALSE)
  cobra_plot@plotcolors <- cobra_plot@plotcolors[1:2]
  col = colours[-length(colours)]
  names(col) = methods_order
  (gg_roc <- plot_roc(cobra_plot,linewidth=2)  + 
      #scale_fill_brewer(breaks=methods, type="qual", palette=3)+
      # scale_color_manual(values = cobra_plot@plotcolors, name = "", labels = methods_order,limits = force) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                       hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            legend.position = "top") )
  gg_roc <- gg_roc +
    guides(color=guide_legend(nrow=1,byrow=TRUE))
  
  #ggsave(paste0(path,roc_save_name),width = 7, height = 5)
  col = colours[-length(colours)]
  names(col) = methods_order
  # plot FDR/TPR curve
  gg_fdr <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                             pointsize = 0, linewidth = 2)+ 
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.8)) +
    scale_color_manual(values = col,
                       name = "",
                       breaks=methods_order,
                       labels=methods_order)+
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(3)),
          axis.text.y=element_text(size=rel(3)),
          axis.title.y = element_text(size=rel(3)),
          axis.title.x = element_text(size=rel(3)),
          legend.title=element_text(size=rel(2)),
          legend.text=element_text(size=rel(3)),
          legend.key.width=unit(2, "cm"),
          aspect.ratio = 1, legend.position="bottom",
          legend.box="vertical", legend.margin=margin())  +
    guides(colour = guide_legend(ncol = 2, byrow = FALSE,
                                 override.aes = list(shape = shape_fill,
                                                     fill = colours[-length(colours)]) ) ) +
    geom_point(size = 10, aes(fill = method, colour =method, shape = method), 
               shape = rep(shape_border[c(2,1)],4), stroke = 2, alpha = 1) + # stroke = line width
    scale_fill_manual(values =col, guide = "none") +
    geom_point(size = 10, aes(fill = method, colour = method, shape = method), 
               shape = rep(shape_fill[c(2,1)],4), stroke = 2, alpha = 0.25)
  
  return(list(gg_roc, gg_fdr))
  
}
overall_roc_fdr_plot_multi_samples_DESpace <- function(
    path = "~",
    spatial_probs = c(0.5,0.8),
    # roc_save_name='_ROC.pdf',
    # fdr_save_name = '_FDR.pdf',
    colours = colours,
    shape_border = c(0, 1),
    shape_fill = c(15, 16),
    dataset = "LIBD",
    pattern_names#,
    #sample_names
){
  path_data = paste0(path, dataset, pattern_names, "/probs_", spatial_probs[1], "_", spatial_probs[2], "FALSE/")
  # for(q in seq_along(sample_names)) { 
  #  path_sub = paste0(path_data, "probs_", spatial_probs[1], "_", spatial_probs[2], "_final_object",sample_names[q], '.rda')
  #  load(path_sub)
  #  if(!exists("combined_object")){combined_object = final_object; remove(final_object)}
  #  assign(paste0("sce", q), combined_object)
  # }
  #  gene1 <- rownames(sce1);gene2 <- rownames(sce2);gene3 <- rownames(sce3)
  #  sum(gene1 == gene2)== 5000; sum(gene1 == gene3)== 5000
  
  load(paste0(path_data,"result_SV_edgeR_counts.rda"))
  result = result1 = spatial_gene_results
  file = paste0(path_data, 'probs_',spatial_probs[2], '_selected_genes.txt')
  genes <- read.csv(file,header = TRUE,sep =  '\t')
  all_genes <- as.data.frame((unique(result$genes)))
  
  all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
  genes_id <- apply(genes,1,function(x) substr(x, 10,18))
  status <- ifelse(all_genes_id %in% genes_id, 1, 0)
  truth <- cbind(all_genes, status)
  # check
  sum(truth$status) == dim(genes)[1]
  
  result <- result %>% dplyr::distinct()
  setDT(result)
  ################################
  methods_order = c('Single_Sample_DESpace','Multiple_Samples_DESpace'#,
                    #'stLearn_single_DESpace','stLearn_multi_DESpace'
  )
  single_Bayes = as.data.table(result1 %>% filter(method == "SV_single_edgeR" & cluster == "Original"))
  multi_Bayes = as.data.table(result1 %>% filter(method == "SV_multi_edgeR" & cluster == "Original"))
  single_Bayes[, methods := 'Single_Sample_DESpace']
  multi_Bayes[, methods := 'Multiple_Samples_DESpace']
  
  single_Bayes = single_Bayes[,.(genes, FDR, PValue, sample, methods)]
  multi_Bayes = multi_Bayes[,.(genes, FDR, PValue, sample, methods)]
  
  samples = unique(single_Bayes$sample)
  combined_data <- data.frame(0)
  combined_data2 <- data.frame(0)
  combined_df2 <- NULL
  combined_result <- NULL
  combined_genes <- NULL
  combined_all_genes <- NULL
  
  for(ii in c(1:length(samples))){
    combine_results = rbind(
      
      single_Bayes %>% filter(sample == samples[ii]),
      multi_Bayes
    )
    setDT(combine_results)
    data <- dcast(combine_results, methods ~ genes,value.var = 'FDR')
    df2 = truth[match(colnames(data)[-1], truth$`(unique(result$genes))`),]
    df2 <- as.data.frame(df2[,2])
    rownames(df2) <- colnames(data)[-1]
    colnames(df2) <- 'status'
    data2 <- dcast(combine_results, methods ~ genes,value.var = 'PValue')
    
    colnames(data) <- paste0(colnames(data), ".",ii)
    colnames(data2) <- paste0(colnames(data2),".", ii)
    rownames(df2) <- paste0(rownames(df2), ".",ii)
    
    colnames(data)[1] <- "method"
    colnames(data2)[1] <- "method"
    
    #data <- data[,1:(dim(data)[2]-1)]
    #data2 <- data2[,1:(dim(data2)[2]-1)]
    combined_data <- cbind(combined_data,data)
    combined_data2 <- cbind(combined_data2,data2)
    combined_df2 <- rbind(combined_df2,df2)
  }
  
  combined_data <- combined_data[,-1]
  combined_data2 <- combined_data2[,-1]
  
  combined_data <- t(combined_data)
  methods <- combined_data[1,]
  combined_data <- combined_data[-1,]
  colnames(combined_data) <- methods
  
  combined_data2 <- t(combined_data2)
  combined_data2 <- combined_data2[-1,]
  colnames(combined_data2) <- methods
  
  genes_id <- rownames(combined_data2)
  combined_data2 <- apply(combined_data2, 2,            # Specify own function within apply
                          function(x) as.numeric(as.character(x)))
  rownames(combined_data2) <- genes_id
  
  combined_data <- apply(combined_data, 2,            # Specify own function within apply
                         function(x) as.numeric(as.character(x)))
  rownames(combined_data) <- genes_id
  combined_data <- combined_data[,match(methods_order, colnames(combined_data))]
  combined_data2 <- combined_data2[,match(methods_order, colnames(combined_data2))]
  
  DF_COBRA <- COBRAData(pval = data.frame(combined_data2
  ),
  padj = data.frame(combined_data
  ),
  truth = data.frame(combined_df2))
  # truth = 1 for SV genes and 0 for uniform genes.
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # ROC, FDR plots:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  perf <- calculate_performance(DF_COBRA, binary_truth = "status",
                                thrs = c(0.01, 0.05, 0.1, 0.2))
  cobra_plot <- prepare_data_for_plot(perf, colorscheme = colours,incloverall = FALSE,
                                      facetted = TRUE,conditionalfill = FALSE)
  cobra_plot@plotcolors <- cobra_plot@plotcolors[1:2]
  
  (gg_roc <- plot_roc(cobra_plot,linewidth=2)  + 
      #scale_fill_brewer(breaks=methods, type="qual", palette=3)+
      # scale_color_manual(values = cobra_plot@plotcolors, name = "", labels = methods_order,limits = force) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                       hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            legend.position = "top") )
  gg_roc <- gg_roc +
    guides(color=guide_legend(nrow=1,byrow=TRUE))
  
  #ggsave(paste0(path,roc_save_name),width = 7, height = 5)
  col = colours[-length(colours)]
  names(col) = methods_order
  # plot FDR/TPR curve
  gg_fdr <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                             pointsize = 0, linewidth = 2)+ 
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.8)) +
    scale_color_manual(values = col,
                       name = "",
                       breaks=methods_order,
                       labels=methods_order)+
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(3)),
          axis.text.y=element_text(size=rel(3)),
          axis.title.y = element_text(size=rel(3)),
          axis.title.x = element_text(size=rel(3)),
          legend.title=element_text(size=rel(2)),
          legend.text=element_text(size=rel(3)),
          legend.key.width=unit(2, "cm"),
          aspect.ratio = 1, legend.position="bottom",
          legend.box="vertical", legend.margin=margin())  +
    guides(colour = guide_legend(ncol = 2, byrow = FALSE,
                                 override.aes = list(shape = shape_fill,
                                                     fill = colours[-length(colours)]) ) ) +
    geom_point(size = 10, aes(fill = method, colour =method, shape = method), 
               shape = rep(shape_border[c(2,1)],4), stroke = 2, alpha = 1) + # stroke = line width
    scale_fill_manual(values =col, guide = "none") +
    geom_point(size = 10, aes(fill = method, colour = method, shape = method), 
               shape = rep(shape_fill[c(2,1)],4), stroke = 2, alpha = 0.25)
  return(list(gg_roc, gg_fdr))
  
}