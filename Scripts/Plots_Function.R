# Settings:
# 5 Methods: edgeR_counts, edgeR_CPMs, SPARK, MERINGUE, SpatialDE;
# Filter: min 50 of non-zero spots per gene;
# FDR = 0.05.

#NULL (permuted) data:
# density of raw p-values;
#Real data:
# UpsetR plot;
# overall number of significant results;
# Jaccard index btw replicates: keep top xx results (~1,000 ? LIBD and ~200 ? melanoma);
# unique patterns: top xx from each methods that are not in top yy from other methods (xx and yy different on  LIBD and melanoma dataset);
######################################################################################################################################################################
######################################################################################################################################################################
########################################### Real data ###################################################################
################################################################################################################################
###########################################################################################################################
library(stringr)
library(ggplot2)
library(reshape2)
library(gridExtra) # for grid.arrange
library(grid) 
library(ggnewscale)
library(dplyr)
library(scales)
library(RColorBrewer)
library(SingleCellExperiment)
library(scater)
library(ggpubr)

library(BiocParallel)
#library(scater)
library(reshape)
library(data.table)
#library(RColorBrewer)
library(cowplot)
#library(ggplot2)
#library(dplyr)
library(iCOBRA)
library("wesanderson")


library(ComplexHeatmap)
library(forcats)

process_data <- function(
    path = '~/Desktop/master_thesis/Data_re_run/Real_data/',
    data = 'melanoma', # or 'LIBD'
    Manual = F, # if there are results from 'Manual clusters'
    threshold = 0.05,
    name_head = -13,
    name_tail = -5,
    ...
){
  file = paste0(path,data,"/")
  edgeR_counts_BayesSpace_files <- list.files(
    file, pattern = "BayesSpace_edgeR_results *.*csv$", full.names = TRUE)
  edgeR_counts_BayesSpace_results <- lapply(edgeR_counts_BayesSpace_files, read.csv)
  
  edgeR_counts_stLearn_files <- list.files(
    file, pattern = "^stLearn_edgeR_results *.*csv$", full.names = TRUE)
  edgeR_counts_stLearn_results <- lapply(edgeR_counts_stLearn_files, read.csv)
  
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
    
    all_results <- rbind(edgeR_counts_Manual_results2,
                         edgeR_counts_BayesSpace_results2, edgeR_counts_stLearn_results2,
                         MERINGUE_results2, 
                         SPARK_results2,nnSVG_results2, 
                         SPARKX_results2, SpatialDE2_results2, SpatialDE_results2, SpaGCN_results2)
    
    list_sig = list(edgeR_counts_Manual_sig,
                    edgeR_counts_BayesSpace_sig, edgeR_counts_stLearn_sig,
                    MERINGUE_sig, 
                    SPARK_sig,nnSVG_sig, 
                    SPARKX_sig, spatialDE2_sig, spatialDE_sig, spaGCN_sig,
                    all_results)
  }else{
    all_results <- rbind(edgeR_counts_BayesSpace_results2, edgeR_counts_stLearn_results2,
                         MERINGUE_results2, 
                         SPARK_results2,nnSVG_results2, 
                         SPARKX_results2, SpatialDE2_results2, SpatialDE_results2, SpaGCN_results2)
    
    list_sig = list(edgeR_counts_BayesSpace_sig, edgeR_counts_stLearn_sig,
                    MERINGUE_sig, 
                    SPARK_sig,nnSVG_sig, 
                    SPARKX_sig, spatialDE2_sig, spatialDE_sig, spaGCN_sig,
                    all_results)
  }
  
  return(list_sig)
}

library(UpSetR)
UpsetR_plot <- function(
    list_sig,
    path = '~/Desktop/master_thesis/Data_re_run/Real_data/plots',
    save_name = "_UpsetR_plot.pdf",
    #methods = c(edgeR_counts, edgeR_CPMs, SPARK, MERINGUE, SpatialDE),
    sample_names,
    Manual = F,
    ...
){
  for (i in 1:length(sample_names)) {
    fn <- paste0(path, "/UpsetR_plot")
    if(!file.exists(fn)) dir.create(fn, recursive = TRUE)
    fn <- paste0(fn,"/", sample_names[i], save_name)
    edgeR_counts_BayesSpace_sig <- list_sig[[1]]
    edgeR_counts_stLearn_sig <- list_sig[[2]]
    #edgeR_CPMs_BayesSpace_sig <- list_sig[[2]]
    MERINGUE_sig <- list_sig[[3]]
    nnSVG_sig <- list_sig[[4]]
    SPARK_sig <- list_sig[[5]]
    SPARKX_sig <- list_sig[[6]]
    spatialDE_sig <- list_sig[[7]]
    spatialDE2_sig <- list_sig[[8]]
    spaGCN_sig <- list_sig[[9]]
    
    pdf(fn, width = 8, height = 6)
    if(Manual == F){
      print(upset(
        fromList(list(BayesSpace_DESpace = as.character(edgeR_counts_BayesSpace_sig[[i]]),
                      StLearn_DESpace = as.character(edgeR_counts_stLearn_sig[[i]]),
                      MERINGUE = as.character(MERINGUE_sig[[i]]), 
                      nnSVG = as.character(nnSVG_sig[[i]]),
                      SPARK = as.character(SPARK_sig[[i]]),
                      SPARKX = as.character(SPARKX_sig[[i]]),
                      SpatialDE = as.character(spatialDE_sig[[i]]),
                      SpatialDE2 = as.character(spatialDE2_sig[[i]]),
                      SpaGCN = as.character(spaGCN_sig[[i]])
                      #edgeR_CPMs = as.character(edgeR_CPMs_BayesSpace_sig[[i]]),
                      #SV_edgeR_NoNorm = as.character(edgeR_NoNorm_sig[[i]]),
                      #SV_edgeR_CPMs_NoNorm = as.character(edgeR_CPMs_NoNorm_sig[[i]]),
                      
                      
        )),
        order.by = "freq", 
        decreasing = TRUE, 
        sets = rev(c("BayesSpace_DESpace","StLearn_DESpace",#"edgeR_CPMs",
                     "MERINGUE","nnSVG",
                     "SPARK","SPARKX",
                     "SpatialDE", "SpatialDE2","SpaGCN" )), 
        keep.order = TRUE, 
        mainbar.y.label = sample_names[i]#,
        #nintersects  = 15
      ))} else if(Manual == T){
        edgeR_counts_Manual_sig <- list_sig[[10]]
        #edgeR_CPMs_Manual_sig <- list_sig[[7]]
        
        print(upset(
          fromList(list(BayesSpace_DESpace = as.character(edgeR_counts_BayesSpace_sig[[i]]),
                        #edgeR_CPMs_BayesSpace = as.character(edgeR_CPMs_BayesSpace_sig[[i]]),
                        Manual_edgeR = as.character(edgeR_counts_Manual_sig[[i]]),
                        #edgeR_CPMs_Manual = as.character(edgeR_CPMs_Manual_sig[[i]]),
                        SpatialDE = as.character(spatialDE_sig[[i]]), 
                        MERINGUE = as.character(MERINGUE_sig[[i]]), 
                        SPARK = as.character(SPARK_sig[[i]])
                        
          )),
          order.by = "freq", 
          decreasing = TRUE, 
          sets = rev(c("BayesSpace_DESpace",#"edgeR_CPMs_BayesSpace",
                       "Manual_edgeR",#"edgeR_CPMs_Manual",
                       "SpatialDE", "MERINGUE","SPARK" )), 
          keep.order = TRUE, 
          mainbar.y.label = sample_names[i]#,
          #nintersects  = 15
        ))
      }
    dev.off()
  }
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

jaccard_btw_methods <- function(
    list_sig,
    path = '~/Desktop/master_thesis/Data_re_run/Real_data/plots',
    save_name = "_jaccard_index_corrplot.pdf",
    sample_names,
    Manual = F,
    top_genes = 200,
    ...
){
  for (i in (1:length(sample_names))){
    print(i)
    edgeR_counts_BayesSpace_sig <- list_sig[[1]]
    #edgeR_CPMs_BayesSpace_sig <- list_sig[[2]]
    MERINGUE_sig <- list_sig[[3]]
    SPARK_sig <- list_sig[[4]]
    spatialDE_sig <- list_sig[[5]]
    
    edgeR_counts_BayesSpace <- as.vector(unlist(edgeR_counts_BayesSpace_sig[i]))
    #edgeR_CPMs_BayesSpace <- as.vector(unlist(edgeR_CPMs_BayesSpace_sig [i]))
    SpatialDE <- as.vector(unlist(spatialDE_sig[i]))
    SPARK <- as.vector(unlist(SPARK_sig[i]))
    MERINGUE <- as.vector(unlist(MERINGUE_sig[i]))
    
    if(Manual == T){
      edgeR_counts_Manual_sig <- list_sig[[6]]
      #edgeR_CPMs_Manual_sig <- list_sig[[7]]
      edgeR_counts_Manual <- as.vector(unlist(edgeR_counts_Manual_sig[i]))
      #edgeR_CPMs_Manual <- as.vector(unlist(edgeR_CPMs_Manual_sig [i]))
      
      list_test <- list(edgeR_counts_BayesSpace[1:min(length(edgeR_counts_BayesSpace),top_genes)],
                        #edgeR_CPMs_BayesSpace[1:min(length(edgeR_CPMs_BayesSpace),top_genes)],
                        edgeR_counts_Manual[1:min(length(edgeR_counts_Manual),top_genes)],
                        #edgeR_CPMs_Manual[1:min(length(edgeR_CPMs_Manual),top_genes)],
                        SpatialDE[1:min(length(SpatialDE),top_genes)], 
                        MERINGUE[1:min(length(MERINGUE),top_genes)], 
                        SPARK[1:min(length(SPARK),top_genes)])
      
      row5 <- mapply(jaccard,list_test,list_test[5])
      #row7 <- mapply(jaccard,list_test,list_test[7])
    }
    
    
    if(Manual == F){list_test <- list(edgeR_counts_BayesSpace[1:min(length(edgeR_counts_BayesSpace),top_genes)],
                                      #edgeR_CPMs_BayesSpace[1:min(length(edgeR_CPMs_BayesSpace),top_genes)],
                                      SpatialDE[1:min(length(SpatialDE),top_genes)], 
                                      MERINGUE[1:min(length(MERINGUE),top_genes)], 
                                      SPARK[1:min(length(SPARK),top_genes)])}
    
    row1 <- mapply(jaccard,list_test,list_test[1])
    row2 <- mapply(jaccard,list_test,list_test[2])
    row3 <- mapply(jaccard,list_test,list_test[3])
    row4 <- mapply(jaccard,list_test,list_test[4])
    #row5 <- mapply(jaccard,list_test,list_test[5])
    
    if(Manual == T){
      mat_jaccord <- as.data.frame(cbind(row1,row2,row3,row4,row5))
      colnames(mat_jaccord) <- c("BayesSpace_DESpace","Manual_edgeR",
                                 "SpatialDE","MERINGUE","SPARK")}
    else if (Manual == F){
      mat_jaccord <- as.data.frame(cbind(row1,row2,row3,row4))
      colnames(mat_jaccord) <- c("BayesSpace_DESpace","SpatialDE","MERINGUE","SPARK")}
    
    rownames(mat_jaccord) <- colnames(mat_jaccord)
    lower_tri <- get_lower_tri(mat_jaccord)
    melted_cormat <- reshape2::melt(as.matrix(lower_tri), na.rm = TRUE)
    # Heatmap
    
    title = paste0("Sample ", sample_names[i])
    pdf(file=paste0(path,"/",title,save_name), onefile=FALSE) # or other device
    corrplot::corrplot(as.matrix(mat_jaccord),type="upper",tl.col="black",
                       mar=c(0,0,1,0), col.lim = c(0,1),tl.cex = 0.8)
    
    dev.off()
  }
}
#devtools::install_github("Laurae2/Laurae")

jaccord_btw_rep <- function(
    list_test,
    path = '~/Desktop/master_thesis/Data_re_run/Real_data/plots',
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

#hist-like plot with the average jaccard index (1 for LIBD and 1 for melanoma).

ave_jaccard_btw_rep <- function(
    list_test_all,
    path = '~/Desktop/master_thesis/Data_re_run/Real_data/plots',
    #save_name = "_jaccard_index_corrplot.pdf",
    sample_names=sample_names,
    top_genes = 200,
    #width = 12, height = 10,
    data = 'melanoma',
    colors_method = colors_method,
    #Manual = F,
    ...
){
  for(j in c(1:9)){
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
                          mean(method9[upper.tri(method9, diag = FALSE)])
                          
    ))
    
    colnames(df) <- c("Value")
    rownames(df) <- c("BayesSpace_DESpace","StLearn_DESpace",
                      "MERINGUE", 
                      "SPARK","nnSVG","SPARK-X",
                      "SpatialDE2", "SpatialDE",
                      "SpaGCN")
    colors_method <- colors_method
    df <- reshape::melt(t(df))
    df$X2 <- as.factor(df$X2)
    #all_colours = c(rep("#009E73",1), rep("#CC79A7",1), rep("#0072B2",1), rep("#D55E00",1), rep("#CD853F",1))
    gg_data <- df[,2:3]
    colnames(gg_data) <- c( "Method", "value")
    breaks_plots = c(0.05,0.1,0.15,0.2)
    gg_data$Method <- factor(gg_data$Method,levels=c("BayesSpace_DESpace","StLearn_DESpace",
                                                     "MERINGUE", "nnSVG",
                                                     "SPARK","SPARK-X",
                                                     "SpatialDE", "SpatialDE2",
                                                     "SpaGCN"))
  
  
  gg_Bar <- ggplot(gg_data,
         #aes_string(y = "Time_Minutes", fill = "Method"),
         aes(x = reorder(Method,-value), y = value, fill = Method)) + 
    
    #ggplot(aes_string(x = "Method", y = "value", fill = "Method")) +  
    geom_bar( stat = "identity",position = position_dodge(0.8), width = 0.7) +  theme_bw() +   xlab("") +  
  ylab("") +
     # ylab("Average Jaccaed index across samples") + #facet_wrap(~Sample_id) +
    #scale_y_sqrt( breaks = c(0.1,0.15,0.25,0.4,0.5,0.6) ) + 
    #scale_y_continuous( breaks = c(0.05,0.1,0.15,0.2) ) + 
    scale_fill_manual("Method", values = colors_method)+
    scale_y_continuous( breaks = breaks_plots ) +
    theme(axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1),   
          axis.text.y= element_text(size=rel(1.8),angle = 0, hjust = 1), 
          axis.title = element_blank(),#element_text(size=15,face = "bold"),      
          panel.grid.minor = element_blank(),        
          panel.grid.major.x = element_blank(),       
          panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),        
          aspect.ratio = 1,        
          legend.position = "none",
          legend.text=element_text(size=15))
  return(gg_Bar)
  #dev.off()
}
ave_jaccard_btw_rep_LIBD <- function(
    list_test_all,
    path = '~/Desktop/master_thesis/Data_re_run/Real_data/plots',
    #save_name = "_jaccard_index_corrplot.pdf",
    sample_names=sample_names,
    top_genes = 200,
    #width = 12, height = 10,
    data = 'melanoma',
    colors_method = colors_method,
    #Manual = F,
    ...
){
  for(j in c(1:4, 6:10)){
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
                        #mean(method5[upper.tri(method5, diag = FALSE)]),
                        mean(method6[upper.tri(method6, diag = FALSE)]),
                        mean(method7[upper.tri(method7, diag = FALSE)]),
                        mean(method8[upper.tri(method8, diag = FALSE)]),
                        mean(method9[upper.tri(method9, diag = FALSE)]),
                        mean(method10[upper.tri(method10, diag = FALSE)])
                        
  )
  )
  
  colnames(df) <- c("Value")
  rownames(df) <- c("Manual_DESpace",
                    "BayesSpace_DESpace",
                    "StLearn_DESpace",#"edgeR_CPMs_BayesSpace", 
                    "MERINGUE", #"SPARK",
                    "nnSVG",
                    "SPARK-X", "SpatialDE2",
                    "SpatialDE","SpaGCN")
  # colors_method <- c(
  #   "BayesSpace_edgeR" = "#A6CEE3",#"#FFDB6D", #"#C4961A",
  #   "StLearn_edgeR" = "#1F78B4",
  #   "Manual_edgeR" = "#08306B",
  #   "MERINGUE" = "#F4A582",#"#FC9272",#"#EF3B2C",#"#E7298A",#"#FFFF99",#"#F4EDCA",#"#D16103", 
  #   #"#9E9AC8",#"#6A51A3",#"#B15928",
  #   #"SPARK" = "#FFD92F",#"#E6AB02",#"#CAB2D6",#"#C3D7A4",#"#52854C",
  #   "SPARKX" = "#A6761D",#"#6A3D9A",
  #   "SpatialDE" = "#B2DF8A",#"#4E84C4"#"#293352"
  #   "SpatialDE2" = "#33A02C",
  #   "SpaGCN" = "#666666",
  #   "nnSVG" = "#C51B7D"
  # )
  df <- reshape::melt(t(df))
  df$X2 <- as.factor(df$X2)
  #all_colours = c(rep("#009E73",1), rep("#CC79A7",1), rep("#0072B2",1), rep("#D55E00",1), rep("#CD853F",1))
  gg_data <- df[,2:3]
  colnames(gg_data) <- c( "Method", "value")
  breaks_plots = c(0.1,0.2,0.3,0.4,0.5,0.6)
  gg_data$Method <- factor(gg_data$Method,levels=c("Manual_DESpace",#"edgeR_CPMs_Manual",
                                                   "BayesSpace_DESpace",
                                                   "StLearn_DESpace",#"edgeR_CPMs_BayesSpace", 
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
    scale_fill_manual("Method", values = colors_method)+
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
plot_exprs <- function(
    sce_one, 
    method = NULL, 
    gene_list, 
    cbp1 = colorRampPalette(rev(brewer.pal(n=9, name="BrBG")))(100)[55:100], 
    i=1, 
    sample_names=NULL,
    path = "~",
    coordinates = c("col","row"),
    width = 15, height = 8,
    point_size = 2.5,
    save_plot = TRUE,
    legend = TRUE,
    title = TRUE,
    ...
){
  sce_filter = subset(sce_one, rownames(rowData(sce_one)) %in% gene_list)
  xy_coords <- data.frame(
    x_coord = as.numeric(as.character(colData(sce_filter)[, coordinates[1]])), 
    y_coord = as.numeric(as.character(colData(sce_filter)[, coordinates[2]]))
  )
  # extract expression levels (UMI counts)
  exprs <- counts(sce_filter)
  dim(exprs)
  if(is.null(colData(sce_filter)$sizeFactor)){
    sizeFactors(sce_filter) <- colSums(assay(sce_filter))}
  d_plot <- cbind(
    barcode = colData(sce_filter)$sizeFactor, 
    xy_coords, 
    as.data.frame(as.matrix(t(exprs)))
  )
  if (dim(d_plot)[2] > 4){
    d_plot[,4:dim(d_plot)[2]] <- apply(d_plot[,4:dim(d_plot)[2]],2, function(x) if((max(x)-min(x)) != 0){(x-min(x))/(max(x)-min(x))})
  }else if (dim(d_plot)[2] == 4){
    x <- d_plot[,4]
    d_plot[,4] <- ((x-min(x))/(max(x)-min(x)))
  }
  
  d_plot <- reshape2::melt(
    d_plot, 
    id.vars = c("barcode", "x_coord", "y_coord"), 
    variable.name = "gene_id", 
    value.name = "Expression"
  )
  #max_UMI <- max(d_plot$Expression)
  # split into subplots since otherwise too many panels and code too slow
  # up to n_subplots depending on number of panels
  n_subplots <- 10
  
  n_col <- 10
  n_row <- 18
  #d_plot$subplot <- rep(1:n_subplots, each = n_col * n_row *    ncol(sce_low_padj))[1:(nrow(sce_low_padj) * ncol(sce_low_padj))]
  d_plot$subplot <- rep(1:n_subplots, each = n_col * n_row * ncol(sce_filter))[1:dim(d_plot)[1]]
  
  
  for (z in unique(d_plot$subplot)) {
    d_plot_sub <- d_plot[d_plot$subplot == z, ]
    
    plots <- 
      
      ggplot(d_plot_sub, aes(x = x_coord, y = y_coord, color = Expression)) + 
      facet_wrap(~ factor(gene_id), ncol = n_col) +
      #scale_x_discrete(limits = factor(order(levels(d_plot_sub$x_coord)))) +
      #scale_y_discrete(limits = factor((order(levels(d_plot_sub$y_coord))))) +
      geom_point(size = point_size, alpha = 1,shape = 46) + 
      scale_color_gradientn(colors = cbp1) +
      #values = rescale(c(0, max_UMI/2, max_UMI*3/4, max_UMI)), 
      #limits = c(0, max_UMI)) + 
      coord_fixed() + 
      theme_bw() + 
      theme(strip.text = element_text(#face="bold",
        size=10,lineheight=5.0), 
        #strip.placement = "outside",
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.ticks = element_blank()
      ) 
    #scale_y_reverse() +
    if(title == TRUE){
      if(is.null(sample_names)){
        plots <- plots + labs(title = "Expressions of top genes")
      }else{
        plots <- plots + labs(title = paste0("Expressions of top genes identified by ",method,": sample ", sample_names[i]))} 
    }
    
    if(legend == FALSE){
      #plots <- plots + theme_void() 
      plots <- plots + theme(legend.position = "none",
                             strip.background = element_blank(),
                             panel.border = element_blank()) 
      
    }
    
    
    if(save_plot == TRUE){
      filename <- paste0(path,"expression_plots_",method,"_", sample_names[i], ".png")
      print(plots)
      ggsave(filename, width = width, height = height)
    }
  }
  return(plots)
}

plot_exprs_pipe <- function(sce,
                            i,
                            top_genes = 100,
                            plot_genes = 50,
                            path = '~/Desktop/master_thesis/Data_re_run/Real_data/plots',
                            list_sig,
                            sample_names = sample_names,
                            coordinates  = c("col","row"),
                            Manual = F,
                            point_size = 0.05,
                            width = 15, height = 8,
                            ...
){
  sce <- logNormCounts(sce)
  assays(sce)$cpm <- calculateCPM(sce)
  
  colData(sce) = cbind(colData(sce), perCellQCMetrics(sce))
  rowData(sce) = cbind(rowData(sce), perFeatureQCMetrics(sce))
  sce_one = sce
  
  print(i)
  edgeR_counts_BayesSpace_sig <- list_sig[[1]]
  #edgeR_CPMs_BayesSpace_sig <- list_sig[[2]]
  MERINGUE_sig <- list_sig[[3]]
  SPARK_sig <- list_sig[[4]]
  spatialDE_sig <- list_sig[[5]]
  
  if(Manual == T){
    edgeR_counts_Manual_sig <- list_sig[[6]]
    #edgeR_CPMs_Manual_sig <- list_sig[[7]]
    edgeR_counts_Manual <- as.vector(unlist(edgeR_counts_Manual_sig[i]))
    #edgeR_CPMs_Manual <- as.vector(unlist(edgeR_CPMs_Manual_sig [i]))
  }
  
  edgeR_counts_BayesSpace <- as.vector(unlist(edgeR_counts_BayesSpace_sig[i]))
  #edgeR_CPMs_BayesSpace <- as.vector(unlist(edgeR_CPMs_BayesSpace_sig [i]))
  SpatialDE <- as.vector(unlist(spatialDE_sig[i]))
  SPARK <- as.vector(unlist(SPARK_sig[i]))
  MERINGUE <- as.vector(unlist(MERINGUE_sig[i]))
  
  `%notin%` <- Negate(`%in%`)
  
  if(Manual == F){
    # unique for spatialDE
    gene1 <- SpatialDE[SpatialDE %notin% c(edgeR_counts_BayesSpace[1:top_genes],
                                           SPARK[1:top_genes],
                                           MERINGUE[1:top_genes])]
    # unique for SV
    gene2 <- edgeR_counts_BayesSpace[edgeR_counts_BayesSpace %notin% c(SpatialDE[1:top_genes],
                                                                       SPARK[1:top_genes],
                                                                       MERINGUE[1:top_genes])]  
    
    #gene2_2 <- edgeR_CPMs_BayesSpace[edgeR_CPMs_BayesSpace %notin% c(SpatialDE[1:top_genes],
    #                                                    SPARK[1:top_genes],
    #                                                   MERINGUE[1:top_genes])]  
    # unique for meringue
    gene3 <- MERINGUE[MERINGUE %notin% c(edgeR_counts_BayesSpace[1:top_genes],
                                         SPARK[1:top_genes],
                                         SpatialDE[1:top_genes])]
    # unique for spark
    gene4 <- SPARK[SPARK %notin% c(edgeR_counts_BayesSpace[1:top_genes],
                                   SpatialDE[1:top_genes],
                                   MERINGUE[1:top_genes])]
    
    length(gene1); length(gene2); #length(gene2_2); 
    length(gene3); length(gene4)
    
    gene1 <- gene1[1:plot_genes]
    gene2 <- gene2[1:plot_genes]
    #gene2_2 <- gene2_2[1:plot_genes]
    gene3 <- gene3[1:plot_genes]
    gene4 <- gene4[1:plot_genes]
    
    if(!is.na(gene1[1])){plot_exprs(sce_one=sce, method = " SpatialDE",gene_list = gene1, 
                                    i=i, sample_names=sample_names,path=path,coordinates = coordinates,
                                    point_size = point_size,width = width, height = height)}
    
    if(!is.na(gene2[1])){plot_exprs(sce_one=sce, method = " BayesSpace_DESpace",gene_list = gene2, 
                                    i=i, sample_names=sample_names,path=path,coordinates = coordinates,
                                    point_size = point_size,width = width, height = height)}
    
    #if(!is.na(gene2_2[1])){plot_exprs(sce_one=sce, method = " edgeR_CPMs_BayesSpace",gene_list = gene2_2, 
    #  i=i, sample_names=sample_names,path=path,coordinates = coordinates)}
    
    if(!is.na(gene3[1])){plot_exprs(sce_one=sce, method = " MERINGUE",gene_list = gene3, 
                                    i=i, sample_names=sample_names,path=path,coordinates = coordinates,
                                    point_size = point_size,width = width, height = height)}
    
    if(!is.na(gene4[1])){plot_exprs(sce_one=sce, method = " SPARK",gene_list = gene4, 
                                    i=i, sample_names=sample_names,path=path,coordinates = coordinates,
                                    point_size = point_size,width = width, height = height)}
  }
  if(Manual == T){
    # unique for spatialDE
    gene1 <- SpatialDE[SpatialDE %notin% c(edgeR_counts_Manual[1:top_genes],
                                           SPARK[1:top_genes],
                                           MERINGUE[1:top_genes])]
    # unique for SV
    gene2 <- edgeR_counts_Manual[edgeR_counts_Manual %notin% c(SpatialDE[1:top_genes],
                                                               SPARK[1:top_genes],
                                                               MERINGUE[1:top_genes])]  
    
    # gene2_2 <- edgeR_CPMs_Manual[edgeR_CPMs_Manual %notin% c(SpatialDE[1:top_genes],
    #                                                                SPARK[1:top_genes],
    #                                                                MERINGUE[1:top_genes])] 
    
    gene2_3 <- edgeR_counts_BayesSpace[edgeR_counts_BayesSpace %notin% c(SpatialDE[1:top_genes],
                                                                         SPARK[1:top_genes],
                                                                         MERINGUE[1:top_genes])]  
    
    # gene2_4 <- edgeR_CPMs_BayesSpace[edgeR_CPMs_BayesSpace %notin% c(SpatialDE[1:top_genes],
    #                                                        SPARK[1:top_genes],
    #                                                       MERINGUE[1:top_genes])]  
    # unique for meringue
    gene3 <- MERINGUE[MERINGUE %notin% c(edgeR_counts_Manual[1:top_genes],
                                         SPARK[1:top_genes],
                                         SpatialDE[1:top_genes])]
    # unique for spark
    gene4 <- SPARK[SPARK %notin% c(edgeR_counts_Manual[1:top_genes],
                                   SpatialDE[1:top_genes],
                                   MERINGUE[1:top_genes])]
    
    length(gene1); length(gene2); length(gene2_3); length(gene3); length(gene4)
    
    gene1 <- gene1[1:plot_genes]
    gene2 <- gene2[1:plot_genes]
    # gene2_2 <- gene2_2[1:plot_genes]
    gene2_3 <- gene2_3[1:plot_genes]
    #gene2_4 <- gene2_4[1:plot_genes]
    gene3 <- gene3[1:plot_genes]
    gene4 <- gene4[1:plot_genes]
    
    if(!is.na(gene1[1])){plot_exprs(sce_one=sce, method = " SpatialDE",gene_list = gene1, 
                                    i=i, sample_names=sample_names,path=path,coordinates = coordinates,
                                    point_size = point_size,width = width, height = height)}
    
    if(!is.na(gene2[1])){plot_exprs(sce_one=sce, method = " Manual_edgeR",gene_list = gene2, 
                                    i=i, sample_names=sample_names,path=path,coordinates = coordinates,
                                    point_size = point_size,width = width, height = height)}
    
    #    if(!is.na(gene2_2[1])){plot_exprs(sce_one=sce, method = " edgeR_CPMs_Manual",gene_list = gene2_2, 
    #                                    i=i, sample_names=sample_names,path=path,coordinates = coordinates)}
    
    if(!is.na(gene2_3[1])){plot_exprs(sce_one=sce, method = " BayesSpace_DESpace",gene_list = gene2_3, 
                                      i=i, sample_names=sample_names,path=path,coordinates = coordinates,
                                      point_size = point_size,width = width, height = height)}
    
    #  if(!is.na(gene2_4[1])){plot_exprs(sce_one=sce, method = " edgeR_CPMs_BayesSpace",gene_list = gene2_4, 
    #                                  i=i, sample_names=sample_names,path=path,coordinates = coordinates)}
    #
    if(!is.na(gene3[1])){plot_exprs(sce_one=sce, method = " MERINGUE",gene_list = gene3, 
                                    i=i, sample_names=sample_names,path=path,coordinates = coordinates,
                                    point_size = point_size,width = width, height = height)}
    
    if(!is.na(gene4[1])){plot_exprs(sce_one=sce, method = " SPARK",gene_list = gene4, 
                                    i=i, sample_names=sample_names,path=path,coordinates = coordinates,
                                    point_size = point_size,width = width, height = height)}
  }
}

######################################################################################################################################################################
######################################################################################################################################################################
########################################### Simulation ###################################################################
################################################################################################################################
###########################################################################################################################
#Figures Simulated data:
# ROC and FDR (-> and ROC/FDR aggregated across samples);

#library(spatialLIBD)
# library(BiocParallel)
# library(scater)
# library(reshape)
# library(data.table)
# library(RColorBrewer)
# library(cowplot)
# library(ggplot2)
# library(dplyr)
# library(iCOBRA)
# library("wesanderson")
roc_fdr_plot <- function(
    i, 
    sample_names, 
    pattern_names = "BayesSpace_clusters_patch",
    path = '~/Desktop/master_thesis/Data_re_run/simulation/',
    spatial_probs = c(0.5,0.9),
    default = FALSE,
    methods = c("edgeR_counts","edgeR_CPMs","SPARK","MERINGUE","SpatialDE"),
    methods_order = c("SV_edgeR_counts","SV_edgeR_CPMs","spark","meringue","spatialDE"),
    rep_i = 1,
    # roc_save_name='_ROC.pdf',
    # fdr_save_name = '_FDR.pdf',
    pointsize = 2.5,
    dataset = 'melanoma',
    ...
){
  dir = paste0(path,"results/",sample_names[i],'/',pattern_names,'/probs_',spatial_probs[1],'_',
               spatial_probs[2],default,'/')
  #file = paste0(dir,rep_i,'_probs_',spatial_probs[1],'_',
  # spatial_probs[2],'_results.txt')
  file = paste0(dir,'results_all.csv')
  result <- read.csv(file,header = TRUE,sep =  '\t')
  #genes <- read.csv(paste0(dir,'probs_',spatial_probs[2],'_selected_genes.txt'),header = TRUE,sep =  '\t')
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
  data <- dcast(result, method ~ genes,value.var = 'adj.p.value', fun.aggregate=mean)
  #data <- dcast(result, method ~ genes,value.var = 'adj.p.value')
  df2 = truth[match(colnames(data)[-1], truth$`(unique(result$genes))`),]
  df2 <- as.data.frame(df2[,2])
  rownames(df2) <- colnames(data)[-1]
  colnames(df2) <- 'status'
  data2 <- dcast(result, method ~ genes,value.var = 'p.value', fun.aggregate=mean)
  #data2 <- dcast(result, method ~ genes,value.var = 'p.value')
  
  #data <- data %>%
  # slice(match(methods_order, method))
  #data2 <- data2 %>%
  # slice(match(methods_order, method))
  data <- data %>%
    mutate(method =  factor(method, levels = methods_order)) %>%
    arrange(method)
  data2 <- data2 %>%
    mutate(method =  factor(method, levels = methods_order)) %>%
    arrange(method)
  
  pval = data.frame(BayesSpace_DESpace = t(data2[1,-1]),
                    #edgeR_CPMs = t(data2[2,-1]),
                    StLearn_DESpace = t(data2[2,-1]),
                    SPARK = t(data2[3,-1]),
                    SPARK_X = t(data2[4,-1]),
                    # SPARK_X_wrt_clusters = t(data2[5,-1]),
                    MERINGUE = t(data2[5,-1]),
                    SpaGCN = t(data[6,-1]),
                    nnSVG = t(data2[7,-1]),
                    # nnSVG_wrt_clusters = t(data2[9,-1]),
                    SpatialDE = t(data2[8,-1]),
                    SpatialDE2 = t(data2[9,-1])
  )
  colnames(pval) <- c("BayesSpace_DESpace",#"edgeR_CPMs",
                      "StLearn_DESpace",
                      "SPARK",
                      "SPARK-X",#"SPARK-X_wrt_clusters",
                      "MERINGUE", "SpaGCN","nnSVG",#"nnSVG_wrt_clusters",
                      "SpatialDE","SpatialDE2")
  padj = data.frame(BayesSpace_DESpace = t(data[1,-1]),
                    #edgeR_CPMs = t(data[2,-1]),
                    StLearn_DESpace = t(data[2,-1]),
                    SPARK = t(data[3,-1]),
                    SPARK_X = t(data[4,-1]),
                    #SPARK_X_wrt_clusters = t(data[5,-1]),
                    MERINGUE = t(data[5,-1]),
                    SpaGCN = t(data[6,-1]),
                    nnSVG = t(data[7,-1]),
                    #nnSVG_wrt_clusters = t(data[9,-1]),
                    SpatialDE = t(data[8,-1]),
                    SpatialDE2 = t(data[9,-1])
  )
  colnames(padj) <- c("BayesSpace_DESpace",#"edgeR_CPMs",
                      "StLearn_DESpace","SPARK",
                      "SPARK-X",#"SPARK-X_wrt_clusters",
                      "MERINGUE", "SpaGCN","nnSVG",#"nnSVG_wrt_clusters",
                      "SpatialDE","SpatialDE2")
  DF_COBRA <- COBRAData(#pval,
    padj,
    truth = data.frame(df2))
  perf <- calculate_performance(DF_COBRA, binary_truth = "status")
  cobra_plot <- prepare_data_for_plot(perf, colorscheme = "Paired",incloverall = FALSE,
                                      facetted = TRUE,conditionalfill = FALSE,
                                      keepmethods = methods
  )
  # "Accent"
  # plot ROC curve
  #library(RColorBrewer)
  #brewer.pal(8, "Set3")
  #cobra_plot@plotcolors <- brewer.pal(8, "Set3")
  #cobra_plot@plotcolors <- c("#7FC97F", "#BEAED4", "#FDC086", "#386CB0", "#F0027F", "#BF5B17", "#666666")
  #cobra_plot@plotcolors <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#B3DE69","#FCCDE5","#FDB462")
  cobra_plot@plotcolors <- c("#A6CEE3", "#1F78B4", "#FFD92F", "#A6761D",
                                      "#F4A582", "#666666", "#C51B7D", 
                                      "#B2DF8A", "#33A02C")
                                      
  names(cobra_plot@plotcolors) <- c("BayesSpace_DESpace","StLearn_DESpace",
                                    "SPARK","SPARK-X","MERINGUE","SpaGCN","nnSVG",
                                    "SpatialDE","SpatialDE2")
  
  
  folder_path = paste0(path,"plots/",dataset,"/",'probs_',spatial_probs[1],'_',
                       spatial_probs[2],default, "/", pattern_names,"/")
  if(!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
  print(folder_path)
  p <- plot_roc(cobra_plot)
  
  
  # (p <- p+ scale_colour_brewer(name = "Methods", palette = "Set3",
  #                              breaks=methods,
  #                              labels=methods) +
  #    theme(legend.position = "bottom"))
  #     
  # (p <- p+scale_color_dutchmasters(palette = "pearl_earring", discrete = TRUE) +
  #          theme(legend.position = "bottom"))
  #my.cols <- piratepal(palette = "pony")
  # my.cols <- c("#EB5291FF", "#FBBB68FF", "#F5BACFFF", "#9DDAF5FF", "#6351A0FF", "#525052FF", 
  #              "#FEF79EFF", "#1794CEFF", 
  #              "#16A08CFF" )
  # my.cols = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
  #             "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
  #             "#CAB2D6", "#6A3D9A","#FFFF99")
  my.cols = c("#A6CEE3", "#1F78B4", "#FFD92F", "#A6761D",
                       "#F4A582", "#666666", "#C51B7D", 
                       "#B2DF8A", "#33A02C")
                       (p <- p + scale_color_manual(values = my.cols,
                                                    name = "Methods", 
                                                    breaks=methods,
                                                    labels=methods) +
                           theme(legend.position = "bottom"))
  #          
  if(spatial_probs[1] ==  0.6){print(p + labs(y="TPR_strong", x = "TPR_weak"))}  
  
  ggsave(paste0(folder_path,sample_names[i],'_','methods',roc_save_name), width = 10,height = 7
  )
  
  # plot FDR/TPR curve
  plot_fdrtprcurve(cobra_plot,pointsize = pointsize)+
    scale_x_sqrt(breaks = c(0.00, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) +
    scale_color_manual(values = my.cols,name = "Methods",
                       breaks=methods,
                       labels=methods)+
    theme(legend.position = "bottom")
  # scale_colour_brewer(name = "Methods",palette = "Accent",
  #                     breaks=methods,
  #                     labels=methods)
  #scale_color_manual(values = plotcolors(cobra_plot), labels = methods_order,name = "", limits = force) 
  
  #scale_y_continuous(trans = pow_trans(2))  
  ggsave(paste0(folder_path,sample_names[i],'_','methods',fdr_save_name), width = 10,height = 7
  )
  #saveRDS(list(cobradata = cobra_plot, pattern = pattern_names),
  #       paste0(folder_path,
  #              'methods',".rds"))
  
  DF_COBRA <- COBRAData(pval,
                        padj,
                        truth = data.frame(df2))
  perf <- calculate_performance(DF_COBRA, binary_truth = "status")
  cobra_plot <- prepare_data_for_plot(perf, colorscheme = "Accent",incloverall = FALSE,
                                      facetted = TRUE,conditionalfill = FALSE,
                                      keepmethods = methods
  )
  
  plot_fdrtprcurve(cobra_plot,pointsize = pointsize)+
    scale_x_sqrt(breaks = c(0.00, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) +
    scale_color_manual(values = my.cols,
                       name = "Methods",
                       breaks=methods,
                       labels=methods)+
    theme(legend.position = "bottom")
  # scale_colour_brewer(name = "Methods",palette = "Accent",
  #                     breaks=methods,
  #                     labels=methods)
  
  ggsave(paste0(folder_path,sample_names[i],'_','methods_with_pval',fdr_save_name),width = 10,height = 7
  )
  
  # plot ROC curve
}

overall_roc_fdr_plot <- function(
    sample_names, 
    pattern_names = "BayesSpace_clusters_patch",
    path = '~/Desktop/master_thesis/Data_re_run/simulation/',
    spatial_probs = c(0.5,0.9),
    methods_all = c("BayesSpace_DESpace","StLearn_DESpace",
                    "SPARK","SPARK-X","MERINGUE",
                    "SpaGCN","nnSVG","SpatialDE","SpatialDE2"
    ),
    methods_order = c("SV_edgeR_counts","StLearn_DESpace",
                      "spark", "spark_x","meringue","SpaGCN","nnSVG",
                      "spatialDE","spatialDE2"
    ),
    default = FALSE,
    rep_i = 1,
    # roc_save_name='_ROC.pdf',
    # fdr_save_name = '_FDR.pdf',
    pointsize = 2.5,
    dataset = 'melanoma',
    colours = colours,#path_save = path_save,
    shape_border = c(0, 1, 2, 5, 6, 3, 4, 8, 11),
    shape_fill = c(15, 16, 17, 23, 25, 3, 4, 8, 11),
    ...
){
  combined_data <- data.frame(0)
  combined_data2 <- data.frame(0)
  combined_df2 <- NULL
  combined_result <- NULL
  combined_genes <- NULL
  combined_all_genes <- NULL
  
  for(i in c(1:length(sample_names))){
    dir = paste0(path,sample_names[i],'/',pattern_names,'/probs_',spatial_probs[1],'_',
                 spatial_probs[2],default,'/')
    #file = paste0(dir,rep_i,'_probs_',spatial_probs[1],'_',
    #              spatial_probs[2],'_results.txt')
    file = paste0(dir,'results_all.csv')
    result <- read.csv(file,header = TRUE,sep =  '\t')
    #result <- result %>% filter(method != "spark")
    genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
    all_genes <- as.data.frame((unique(result$genes)))
    
    all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
    genes_id <- apply(genes,1,function(x) substr(x, 10,18))
    status <- ifelse(all_genes_id %in% genes_id, 1, 0)
    truth <- cbind(all_genes, status)
    # check
    sum(truth$status) == dim(genes)[1]
    if(all_genes_id[length(all_genes_id)] == "_count"){
      truth <- truth[1:dim(truth)[1]-1,]
    }
    
    setDT(result)
    data <- dcast(result, method ~ genes,value.var = 'adj.p.value', fun.aggregate=mean)
    df2 = truth[match(colnames(data)[-1], truth$`(unique(result$genes))`),]
    df2 <- df2[1:(dim(df2)[1]-1),]
    df2 <- as.data.frame(df2[,2])
    rownames(df2) <- colnames(data)[-c(1,length(colnames(data)))]
    
    colnames(df2) <- 'status'
    data2 <- dcast(result, method ~ genes,value.var = 'p.value', fun.aggregate=mean)
    
    colnames(data) <- paste0(colnames(data), ".",i)
    colnames(data2) <- paste0(colnames(data2),".", i)
    rownames(df2) <- paste0(rownames(df2), ".",i)
    
    colnames(data)[1] <- "method"
    colnames(data2)[1] <- "method"
    
    data <- data[,1:(dim(data)[2]-1)]
    data2 <- data2[,1:(dim(data2)[2]-1)]
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
  
  # have to check whether there are 2 rows showing method's name
  #combined_data <- combined_data[-c(1001,2002),]
  # combined_data <- combined_data[-c(5000,10000),] ### ONLY FOR LIBD
  
  combined_data2 <- t(combined_data2)
  combined_data2 <- combined_data2[-1,]
  colnames(combined_data2) <- methods
  #combined_data2 <- combined_data2[-c(1001,2002),]
  #  combined_data2 <- combined_data2[-c(5000,10000),] ### ONLY FOR LIBD
  
  genes_id <- rownames(combined_data2)
  combined_data2 <- apply(combined_data2, 2,            # Specify own function within apply
                          function(x) as.numeric(as.character(x)))
  rownames(combined_data2) <- genes_id
  
  combined_data <- apply(combined_data, 2,            # Specify own function within apply
                         function(x) as.numeric(as.character(x)))
  rownames(combined_data) <- genes_id
  combined_data <- combined_data[,match(methods_order, colnames(combined_data))]
  combined_data2 <- combined_data2[,match(methods_order, colnames(combined_data2))]
  colnames(combined_data) <- methods_all
  colnames(combined_data2) <- methods_all
  
  combined_data <- as.data.frame(combined_data[,c(1:length(methods_all))])
  combined_data2 <- as.data.frame(combined_data2[,c(1:length(methods_all))])
  
  DF_COBRA <- COBRAData(pval = data.frame(combined_data2
  ),
  padj = data.frame(combined_data
  ),
  truth = data.frame(combined_df2))
  colnames(DF_COBRA@pval)[4] <- 'SPARK-X'
  ## otherwise it would be 'SPARK.X'
  colnames(DF_COBRA@padj)[4] <- 'SPARK-X'
  # truth = 1 for SV genes and 0 for uniform genes.
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # ROC, FDR plots:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  perf <- calculate_performance(DF_COBRA, binary_truth = "status",
                                thrs = c(0.01, 0.05, 0.1, 0.2))
  cobra_plot <- prepare_data_for_plot(perf, colorscheme = colours, incloverall = FALSE,
                                      facetted = TRUE,conditionalfill = FALSE,
                                      keepmethods = c("BayesSpace_DESpace","StLearn_DESpace","SPARK","SPARK-X","SpatialDE","SpatialDE2","MERINGUE","SpaGCN","nnSVG"))
  # plot ROC curve
  # folder_path = paste0(path_save,dataset,"/",'AggregatedSamples_',spatial_probs[1],'_',
  #                      spatial_probs[2],default, "_", pattern_names,"/")
  # if(!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
  
  (gg_roc <- plot_roc(cobra_plot))
  
  #### New color scheme
  my.cols <-  colours
                         
   (gg_roc <- gg_roc + scale_color_manual(values = my.cols,
                               name = "", 
                               breaks=methods_all,
                               labels=methods_all) +
      theme(legend.position = "bottom")) 
  
  ####
  if(spatial_probs[1] == 0.6){print(gg_roc + labs(y="TPR_strong", x = "TPR_weak"))}
  # ggsave(paste0(path_save,dataset,'/AggregatedSamples_',spatial_probs[1],'_',spatial_probs[2],roc_save_name),
  #        width =  3,height = 4
  # )
  ## Manually match "method" and point "shape"
  sel_shape = c(1,7,9,8,3,4,5,6,2)
  col = my.cols[-length(my.cols)]
  names(col) = methods_all
  #col = col[c(1,9,5,6,7,8,2,4,3)]
  # plot FDR/TPR curve
  gg_fdr <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                   pointsize = 0, linewidth = 2)+
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.8)) +
    scale_color_manual(values = my.cols,
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
          legend.box="vertical", legend.margin=margin())  +
    guides(colour = guide_legend(ncol = 2, byrow = FALSE,
                                 override.aes = list(shape = shape_fill,
                                                     fill = colours[-length(colours)]) ) ) +
     geom_point(size = 6, aes(fill = method, colour = method, shape = method), 
                shape = rep(shape_border[sel_shape],4), stroke = 2, alpha = 1) + # stroke = line width
    scale_fill_manual(values =col, guide = "none") +
    geom_point(size = 6, aes(fill = method, colour = method, shape = method), 
               shape = rep(shape_fill[sel_shape],4), stroke = 2, alpha = 0.25)

  # plot_data_points <- fdrtpr(cobraplot)                 
return(list(gg_roc, gg_fdr))
  # ggsave(paste0(folder_path,'/AggregatedSamples_','probs_',spatial_probs[1],'_',spatial_probs[2],fdr_save_name),
  #        width = 3,height = 4
  # )
  # saveRDS(list(cobradata = cobra_plot, pattern = pattern_names),
  #         paste0(folder_path,
  #                "/AggregatedSamples_probs_",spatial_probs[1],'_',spatial_probs[2],".rds"))
  
}

roc_fdr_plot_cerebellum <- function(
    #i, 
  #sample_names, 
  pattern_names = "BayesSpace_clusters_patch",
  path = '~/Desktop/master_thesis/Data_re_run/simulation/',
  spatial_probs = c(0.5,0.9),
  default = FALSE,
  methods_all = c("edgeR_counts","edgeR_CPMs","SPARK","MERINGUE","SpatialDE"),
  methods_order = c("SV_edgeR_counts","SV_edgeR_CPMs","spark","meringue","spatialDE"),
  rep_i = 1,
  # roc_save_name='_ROC.pdf',
  # fdr_save_name = '_FDR.pdf',
  pointsize = 2.5,
  colours = colours,#path_save = path_save,
  shape_border = c(0, 1, 2, 5, 6, 3, 4, 8, 11),
  shape_fill = c(15, 16, 17, 23, 25, 3, 4, 8, 11),
  # dataset = 'melanoma',
  ...
){
  dir = paste0(path,"SlideSeq2/",pattern_names,'/probs_',spatial_probs[1],'_',
               spatial_probs[2],default,'/')
  #file = paste0(dir,rep_i,'_probs_',spatial_probs[1],'_',
  # spatial_probs[2],'_results.txt')
  file = paste0(dir,'results_all.csv')
  result <- read.csv(file,header = TRUE,sep =  '\t')
  #genes <- read.csv(paste0(dir,'probs_',spatial_probs[2],'_selected_genes.txt'),header = TRUE,sep =  '\t')
  genes <- read.csv(paste0(dir,'SV_genes.txt'),header = TRUE,sep =  '\t')
  all_genes <- as.data.frame((unique(result$genes)))
  
  #all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
  #genes_id <- apply(genes,1,function(x) substr(x, 10,18))
  # genes$V1 for mix and mix inverted
  if(pattern_names %in% c("mixture_patch", "mixture_reverse_patch")){
    status <- ifelse(all_genes$`(unique(result$genes))` %in% genes$V1, 1, 0)
  }else{
    status <- ifelse(all_genes$`(unique(result$genes))` %in% genes$gene_ids, 1, 0)
  }
  truth <- cbind(all_genes, status)
  # check
  sum(truth$status) == dim(genes)[1]
  if(sum(all_genes[length(all_genes)] == "_count")>=1){
    truth <- truth[1:dim(truth)[1]-1,]
  }
  result <- result %>% dplyr::distinct()
  setDT(result)
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
  
  pval = data.frame(t(data2[,-1]))
  identical(colnames(data2)[-1], rownames(pval))
  print(length(methods_all))
  print(length(methods_order))
  print(dim(pval))
  colnames(pval) <- methods_all
  padj = data.frame(t(data[,-1]))
  colnames(padj) <- methods_all
  identical(colnames(data)[-1], rownames(padj))
  ## SpaGCN do not provide p-values. Here SpaGCN column in pval is actually `1-in_out_group_ratio`.
  ## It contains negative values. So we replace SpaGCN column in pval by SpaGCN column in padj.
  if(!is.null(pval$SpaGCN)){pval$SpaGCN = padj$SpaGCN}
  DF_COBRA <- iCOBRA::COBRAData(pval,
                                padj,
                                truth = data.frame(df2))
  perf <- iCOBRA::calculate_performance(DF_COBRA, binary_truth = "status",
                                        thrs = c(0.01, 0.05, 0.1, 0.2))
  cobra_plot <- iCOBRA::prepare_data_for_plot(perf, colorscheme = colours,incloverall = FALSE,
                                              facetted = TRUE,conditionalfill = FALSE,
                                              keepmethods = methods_all)
  (gg_roc <- iCOBRA::plot_roc(cobra_plot))
  
  #### New color scheme
  my.cols <-  colours
  (gg_roc <- gg_roc + scale_color_manual(values = my.cols,
                                         name = "", 
                                         breaks=methods_all,
                                         labels=methods_all) +
      theme(legend.position = "bottom")) 
  #          
  if(spatial_probs[1] ==  0.6){print(p + labs(y="TPR_strong", x = "TPR_weak"))}  
  
  ## Manually match "method" and point "shape"
  if(length(methods_all) == 9) {
    sel_shape = c(1,7,9,8,3,4,5,6,2)
  }else{
    sel_shape = c(1,7,8,3,4,5,6,2)
  }
  col = my.cols[-length(my.cols)]
  names(col) = methods_all
  #col = col[c(1,9,5,6,7,8,2,4,3)]
  # plot FDR/TPR curve
  gg_fdr <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                             pointsize = 0, linewidth = 2)+
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.8)) +
    scale_color_manual(values = my.cols,
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
          legend.box="vertical", legend.margin=margin())  +
    guides(colour = guide_legend(ncol = 2, byrow = FALSE,
                                 override.aes = list(shape = shape_fill,
                                                     fill = colours[-length(colours)]) ) ) +
    geom_point(size = 6, aes(fill = method, colour = method, shape = method), 
               shape = rep(shape_border[sel_shape],4), stroke = 2, alpha = 1) + # stroke = line width
    scale_fill_manual(values =col, guide = "none") +
    geom_point(size = 6, aes(fill = method, colour = method, shape = method), 
               shape = rep(shape_fill[sel_shape],4), stroke = 2, alpha = 0.25)
  
  # plot_data_points <- fdrtpr(cobraplot)                 
  return(list(gg_roc, gg_fdr))
}

roc_fdr_plot_multi_samples <- function(
    path = "~/Desktop/master_thesis/Paper/Data/Simulated/MultiSample/",
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
  file = paste0(path_data, 'probs_0.8_selected_genes.txt')
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
    geom_point(size = 6, aes(fill = method, colour =method, shape = method), 
               shape = rep(shape_border[c(2,1)],4), stroke = 2, alpha = 1) + # stroke = line width
    scale_fill_manual(values =col, guide = "none") +
    geom_point(size = 6, aes(fill = method, colour = method, shape = method), 
               shape = rep(shape_fill[c(2,1)],4), stroke = 2, alpha = 0.25)
  return(list(gg_roc, gg_fdr))
  
}
#Sensitivity/Robustness analysis - consider alternatives to BayesSpace in a sub-set of the simulated data (the LIBD and melanoma manual/BayesSpace SV patterns):
#BayesSpace with n_clusters = 10;
#stLearn (set parameters approx. as in BayesSpace?);
#Internal-only comparison; 
# in each plot: counts_BayesSpace, counts_BayesSpace_10_clust, counts_stLearn_pca_kmeans, counts_stLearn_louvain, 
# CPMs_BayesSpace, CPMs_BayesSpace_10_clust, CPMs_stLearn_pca_kmeans, CPMs_stLearn_louvain, 

Robustness_roc_fdr_plot <- function(
    i, 
    sample_names, 
    pattern_names = "Manual_clusters_patch",
    path = '~/Desktop/master_thesis/Data_re_run/simulation/',
    spatial_probs = c(0.5,0.9),
    methods = c('counts_BayesSpace', 'counts_BayesSpace_10_clust', 'counts_stLearn_pca_kmeans', 'counts_stLearn_louvain'#, 
                #'CPMs_BayesSpace', 'CPMs_BayesSpace_10_clust', 'CPMs_stLearn_pca_kmeans', 'CPMs_stLearn_louvain'
    ),
    methods_order = c('counts_BayesSpace', 'counts_BayesSpace_10_clust', 'counts_stLearn_pca_kmeans', 'counts_stLearn_louvain'#, 
                      # 'CPMs_BayesSpace', 'CPMs_BayesSpace_10_clust', 'CPMs_stLearn_pca_kmeans', 'CPMs_stLearn_louvain'
    ),
    layer_names = 'layer_guess_reordered', #'spatial.cluster', 
    rep_i = 1,
    roc_save_name='_ROC.pdf',
    fdr_save_name = '_FDR.pdf',
    pointsize = 2.5,
    linewidth = 0.5,
    dataset = 'melanoma',
    default = FALSE,
    ...
){
  dir = paste0(path,"results/",sample_names[i],'/',pattern_names,'/probs_',spatial_probs[1],'_',
               spatial_probs[2])
  BayesSpace_file <- paste0(path,"results/",sample_names[i],'/',pattern_names,'/1_probs_',spatial_probs[1],'_',
                            spatial_probs[2],layer_names,'_results_edgeR_robustness.txt')
  pca_file <- paste0(path,"results/",sample_names[i],'/',pattern_names,'/1_probs_',spatial_probs[1],'_',
                     spatial_probs[2],'stLearn_pca_kmeans_results.txt')
  #pca_file <- paste0(path,"results/",sample_names[i],'/',pattern_names,'/1_probs_',spatial_probs[1],'_',
  #                 spatial_probs[2],'stLearn_louvain_results.txt')
  
  BayesSpace_results <- read.csv(BayesSpace_file,header = TRUE,sep =  '\t')
  pca_results <- read.csv(pca_file,header = TRUE,sep =  '\t')
  
  genes <- read.csv(paste0(dir,FALSE,'/probs_',
                           spatial_probs[2],'_selected_genes.txt'),header = TRUE,sep =  '\t')
  
  BayesSpace_results = as.data.table(BayesSpace_results)
  BayesSpace_results[method=='SV_edgeR_counts', method := 'BayesSpace_DESpace']
  BayesSpace_results[method=='counts_BayesSpace_10_clust', method := 'BayesSpace_10_clust_BayesSpace']
  
  pca_results = as.data.table(pca_results)
  pca_results[method=='SV_edgeR_counts', method := 'StLearn_DESpace']
  #pca_results[method=='SV_edgeR_CPMs', method := 'CPMs_stLearn_pca_kmeans']
  
  result <- rbind(BayesSpace_results,
                  pca_results)
  all_genes <- as.data.frame((unique(result$genes)))
  
  all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
  genes_id <- apply(genes,1,function(x) substr(x, 10,18))
  status <- ifelse(all_genes_id %in% genes_id, 1, 0)
  
  truth <- cbind(all_genes, status)
  
  # check
  sum(truth$status) == dim(genes)[1]
  if(all_genes_id[length(all_genes_id)] == "_count"){
    truth <- truth[1:dim(truth)[1]-1,]
  }
  
  results <- result %>% filter(method %in% methods)
  results <- results %>% dplyr::distinct()
  setDT(results)
  data <- dcast(results, method ~ genes,value.var = 'adj.p.value')
  
  
  df2 = truth[match(colnames(data)[-1], truth$`(unique(result$genes))`),]
  
  df2 <- as.data.frame(df2[,2])
  rownames(df2) <- colnames(data)[-1]
  colnames(df2) <- 'status'
  data2 <- dcast(results, method ~ genes,value.var = 'p.value')
  
  #data <- data %>%
  #slice(match(methods_order, method))
  #data2 <- data2 %>%
  #slice(match(methods_order, method))
  data <- data %>%
    mutate(method =  factor(method, levels = methods_order)) %>%
    arrange(method)
  data2 <- data2 %>%
    mutate(method =  factor(method, levels = methods_order)) %>%
    arrange(method)
  
  pval = data.frame(BayesSpace_DESpace = t(data2[1,-1]), 
                    BayesSpace_10_clust_DESpace = t(data2[2,-1]), 
                    StLearn_DESpace = t(data2[3,-1])#, 
                    #counts_stLearn_louvain = t(data2[4,-1])
                    #CPMs_BayesSpace = t(data2[5,-1]), 
                    #CPMs_BayesSpace_10_clust = t(data2[6,-1]),
                    #CPMs_stLearn_pca_kmeans = t(data2[7,-1]),
                    #CPMs_stLearn_louvain = t(data2[8,-1])
  )
  colnames(pval) <- c('BayesSpace_DESpace', 'BayesSpace_10_clust_DESpace', 'StLearn_DESpace'#, 'counts_stLearn_louvain'#, 
                      #'CPMs_BayesSpace', 'CPMs_BayesSpace_10_clust', 'CPMs_stLearn_pca_kmeans', 'CPMs_stLearn_louvain'
  )
  padj = data.frame(BayesSpace_DESpace = t(data[1,-1]), 
                    BayesSpace_10_clust_DESpace = t(data[2,-1]), 
                    StLearn_DESpace = t(data[3,-1])#, 
                    #counts_stLearn_louvain = t(data[4,-1])
                    #CPMs_BayesSpace = t(data[5,-1]), 
                    #CPMs_BayesSpace_10_clust = t(data[6,-1]),
                    #CPMs_stLearn_pca_kmeans = t(data[7,-1]),
                    #CPMs_stLearn_louvain = t(data[8,-1])
  )
  colnames(padj) <- c('BayesSpace_DESpace', 'BayesSpace_10_clust_DESpace', 'StLearn_DESpace'#, 'counts_stLearn_louvain'#, 
                      #'CPMs_BayesSpace', 'CPMs_BayesSpace_10_clust', 'CPMs_stLearn_pca_kmeans', 'CPMs_stLearn_louvain'
  )
  DF_COBRA <- COBRAData(pval,
                        padj,
                        truth = data.frame(df2))
  
  
  perf <- calculate_performance(DF_COBRA, binary_truth = "status")
  cobra_plot <- prepare_data_for_plot(perf, colorscheme = "Set2",incloverall = FALSE,
                                      facetted = TRUE,conditionalfill = FALSE,
                                      keepmethods = methods)
  #plotcolors <- cobra_plot@plotcolors
  # plot ROC curve
  folder_path = paste0(path,"plots/",dataset,"/",pattern_names,"/",sample_names[i],'/probs_',spatial_probs[1],'_',
                       spatial_probs[2],default)
  
  if(!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
  
  plot_roc(cobra_plot,linewidth=linewidth)  + 
    scale_fill_brewer(breaks=methods, type="qual", palette=3)+
    scale_color_manual(values = plotcolors(cobra_plot), name = "", limits = force) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.position = "top") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  ggsave(paste0(folder_path,'/Robustness',roc_save_name),width = 7, height = 5)
  
  # plot FDR/TPR curve
  plot_fdrtprcurve(cobra_plot,pointsize = pointsize,linewidth=linewidth) + 
    scale_x_sqrt(breaks = c(0.00, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) +
    # scale_y_continuous(trans = pow_trans(2))  +
    scale_color_manual(values = plotcolors(cobra_plot), name = "", limits = force) +
    scale_fill_manual(values = plotcolors(cobra_plot), name = "", 
                      limits = force, breaks = c('BayesSpace_DESpace', 'BayesSpace_10_clust_DESpace', 'StLearn_DESpace')) +
    
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.position = "top") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) #+
  # scale_fill_discrete(name = "", labels = c('BayesSpace_DESpace', 'BayesSpace_10_clust_edgeR', 'StLearn_DESpace'))
  
  #scale_x_sqrt(breaks = c(0.00, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) #+
  # scale_y_continuous(trans = pow_trans(2))  
  ggsave(paste0(folder_path,'/Robustness',fdr_save_name),width = 7, height = 5)
  
  saveRDS(list(cobradata = cobra_plot, pattern = pattern_names),paste0(folder_path,"/cobra_Robustness.rds"))
  
  
}

overall_Robustness_roc_fdr_plot <- function(
    sample_names, 
    pattern_names = "Manual_clusters_patch",
    path = '~/Desktop/master_thesis/Data_re_run/simulation/',
    spatial_probs = c(0.5,0.9),
    methods_all = c('BayesSpace_DESpace', 'BayesSpace_10_clust_DESpace', 'StLearn_DESpace'#, 'counts_stLearn_louvain'#, 
                    #'CPMs_BayesSpace', 'CPMs_BayesSpace_10_clust', 'CPMs_stLearn_pca_kmeans', 'CPMs_stLearn_louvain'
    ),
    methods_order = methods,
    layer_names = 'layer_guess_reordered', #'spatial.cluster', 
    rep_i = 1,
    roc_save_name='_ROC.pdf',
    fdr_save_name = '_FDR.pdf',
    pointsize = 2.5,
    linewidth = 0.5,
    dataset = 'melanoma',
    ...
){  
  combined_data <- data.frame(0)
  combined_data2 <- data.frame(0)
  combined_df2 <- NULL
  combined_result <- NULL
  combined_genes <- NULL
  combined_all_genes <- NULL
  
  
  for(i in c(1:length(sample_names))){
    dir = paste0(path,"results/",sample_names[i],'/',pattern_names,'/probs_',spatial_probs[1],'_',
                 spatial_probs[2])
    BayesSpace_file <- paste0(path,"results/",sample_names[i],'/',pattern_names,'/1_probs_',spatial_probs[1],'_',
                              spatial_probs[2],layer_names,'_results_edgeR_robustness.txt')
    pca_file <- paste0(path,"results/",sample_names[i],'/',pattern_names,'/1_probs_',spatial_probs[1],'_',
                       spatial_probs[2],'stLearn_pca_kmeans_results.txt')
    #pca_file <- paste0(path,"results/",sample_names[i],'/',pattern_names,'/1_probs_',spatial_probs[1],'_',
    #                  spatial_probs[2],'stLearn_louvain_results.txt')
    
    BayesSpace_results <- read.csv(BayesSpace_file,header = TRUE,sep =  '\t')
    #BayesSpace_10_results <- read.csv(BayesSpace_10_file,header = TRUE,sep =  '\t')
    #louvain_results <- read.csv(louvain_file,header = TRUE,sep =  '\t')
    pca_results <- read.csv(pca_file,header = TRUE,sep =  '\t')
    
    genes <- read.csv(paste0(dir,FALSE,'/probs_',
                             spatial_probs[2],'_selected_genes.txt'),header = TRUE,sep =  '\t')
    
    BayesSpace_results = as.data.table(BayesSpace_results)
    BayesSpace_results[method=='SV_edgeR_counts', method := 'BayesSpace_DESpace']
    BayesSpace_results[method=='counts_BayesSpace_10_clust', method := 'BayesSpace_10_clust_DESpace']
    
    
    pca_results = as.data.table(pca_results)
    pca_results[method=='SV_edgeR_counts', method := 'StLearn_DESpace']
    # pca_results[method=='SV_edgeR_CPMs', method := 'CPMs_stLearn_pca_kmeans']
    
    result <- rbind(BayesSpace_results,
                    pca_results)
    all_genes <- as.data.frame((unique(result$genes)))
    
    all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
    genes_id <- apply(genes,1,function(x) substr(x, 10,18))
    status <- ifelse(all_genes_id %in% genes_id, 1, 0)
    truth <- cbind(all_genes, status)
    # check
    sum(truth$status) == dim(genes)[1]
    if(all_genes_id[length(all_genes_id)] == "_count"){
      truth <- truth[1:dim(truth)[1]-1,]
    }
    
    results <- result %>% filter(method %in% methods_all)
    setDT(results)
    data <- dcast(results, method ~ genes,value.var = 'adj.p.value')
    df2 = truth[match(colnames(data)[-1], truth$`(unique(result$genes))`),]
    df2 <- as.data.frame(df2[,2])
    rownames(df2) <- colnames(data)[-1]
    colnames(df2) <- 'status'
    data2 <- dcast(results, method ~ genes,value.var = 'p.value')
    
    colnames(data) <- paste0(colnames(data), ".",i)
    colnames(data2) <- paste0(colnames(data2),".", i)
    rownames(df2) <- paste0(rownames(df2), ".",i)
    
    colnames(data)[1] <- "method"
    colnames(data2)[1] <- "method"
    
    data <- data[,1:(dim(data)[2]-1)]
    data2 <- data2[,1:(dim(data2)[2]-1)]
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
  colnames(combined_data) <- methods_all
  colnames(combined_data2) <- methods_all
  
  
  DF_COBRA <- COBRAData(pval = data.frame(combined_data2
  ),
  padj = data.frame(combined_data
  ),
  truth = data.frame(combined_df2))
  # truth = 1 for SV genes and 0 for uniform genes.
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # ROC, FDR plots:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  perf <- calculate_performance(DF_COBRA, binary_truth = "status")
  cobra_plot <- prepare_data_for_plot(perf, colorscheme = "Set2",incloverall = FALSE,
                                      facetted = TRUE,conditionalfill = FALSE,
                                      keepmethods = methods_all)
  
  library(RColorBrewer)
  brewer.pal(8, "Set2")
  cobra_plot@plotcolors <- brewer.pal(8, "Set2")
  names(cobra_plot@plotcolors) <- methods_all
  
  #folder_path = paste0(path,"plots/",sample_names[i],'/',pattern_names,'/probs_',spatial_probs[1],'_',
  #                     spatial_probs[2],FALSE)
  folder_path = paste0(path,"plots/",dataset,"/",pattern_names)
  
  if(!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
  
  plot_roc(cobra_plot,linewidth=linewidth)  + 
    scale_fill_brewer(breaks=methods, type="qual", palette=3)+
    scale_color_manual(values = plotcolors(cobra_plot), name = "", limits = force) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.position = "top") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  
  ggsave(paste0(folder_path,'/overall_Robustness',roc_save_name),width = 7, height = 5)
  
  # plot FDR/TPR curve
  plot_fdrtprcurve(cobra_plot,pointsize = pointsize,linewidth=linewidth)+ 
    scale_fill_brewer(breaks=methods_all, type="qual", palette=3)+
    scale_color_manual(values = plotcolors(cobra_plot), name = "", limits = force) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.position = "top") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))+ 
    #scale_x_continuous(breaks = c(0.0001,0.01,0.05,0.1,1),limits = c(0.0001,1),trans = "log") 
    scale_x_sqrt(breaks = c(0.00, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) #+
  # scale_y_continuous(trans = pow_trans(2))  
  ggsave(paste0(folder_path,'/overall_Robustness',fdr_save_name),width = 7, height = 5)
  
  saveRDS(list(cobradata = cobra_plot, pattern = pattern_names),paste0(folder_path,"/cobra_overall_robustness.rds"))
  
  
}


plot_theme <- function(stripsize, titlecol) {
  theme_grey() +
    theme(legend.position = "right",
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text = element_text(size = stripsize),
          strip.background = element_rect(fill = NA, colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.title = element_text(colour = titlecol))
}
plot_fpr_tpr <- function(cobraplot, title, stripsize, titlecol, pointsize,
                         xaxisrange, aspc) {
  if (aspc == "FPR")
    plot_data <- fpr(cobraplot)
  else if (aspc == "TPR")
    plot_data <- tpr(cobraplot)
  else if (aspc %in% c("SPEARMAN", "PEARSON")) {
    plot_data <- corr(cobraplot)
    plot_data$thr <- as.character(plot_data$thr)
  }
  
  if (!(isTRUE(facetted(cobraplot)))) {
    plot_data$method <- plot_data$fullmethod
  }
  
  ## Number of thresholds
  nthr <- length(unique(plot_data$thr))
  
  pp <- 
    ggplot(plot_data, aes_string(x = aspc, y = "method", group = "method")) +
    geom_point(size = pointsize + 1,
               aes_string(colour = "method", shape = "thr")) +
    scale_shape_manual(values = rep(19, nthr), guide = "none") + 
    scale_color_manual(values = plotcolors(cobraplot), name = "", limits = force) +
    xlim(xaxisrange[1], xaxisrange[2]) +
    plot_theme(stripsize = stripsize, titlecol = titlecol) +
    ggtitle(title)
  if (isTRUE(facetted(cobraplot))) {
    if (!is.finite(maxsplit(cobraplot)))
      msp <- length(unique(plot_data$splitval))
    else
      msp <- maxsplit(cobraplot)
    pp + facet_wrap(~ splitval, nrow = ceiling((msp + 1)/3))
  } else {
    pp
  }
}

#' Plot FPR
#'
#' Plot observed false positive rate (FPR) for given adjusted p-value
#' thresholds.
#'
#' @param cobraplot A \code{COBRAPlot} object.
#' @param title A character string giving the title of the plot.
#' @param stripsize A numeric value giving the size of the strip text, when the
#'   results are stratified by an annotation.
#' @param titlecol A character string giving the color of the title.
#' @param pointsize A numeric value giving the size of the plot characters.
#' @param xaxisrange A numeric vector with two elements, giving the lower and
#'   upper boundary of the x-axis, respectively.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @export
#' @author Charlotte Soneson
#' @examples
#' data(cobradata_example)
#' cobraperf <- calculate_performance(cobradata_example,
#'                                    binary_truth = "status", aspects = "fpr")
#' cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2",
#'                                    incltruth = TRUE)
#' plot_fpr(cobraplot, xaxisrange = c(0, 0.25))
plot_fpr <- function(cobraplot, title = "", stripsize = 15, titlecol = "black",
                     pointsize = 5, xaxisrange = c(0, 1)) {
  plot_fpr_tpr(cobraplot = cobraplot, title = title, stripsize = stripsize,
               titlecol = titlecol, pointsize = pointsize,
               xaxisrange = xaxisrange, aspc = "FPR")
}

#' Plot TPR
#'
#' Plot observed true positive rate (TPR) for given adjusted p-value thresholds.
#'
#' @param cobraplot A \code{COBRAPlot} object.
#' @param title A character string giving the title of the plot.
#' @param stripsize A numeric value giving the size of the strip text, when the
#'   results are stratified by an annotation.
#' @param titlecol A character string giving the color of the title.
#' @param pointsize A numeric value giving the size of the plot characters.
#' @param xaxisrange A numeric vector with two elements, giving the lower and
#'   upper boundary of the x-axis, respectively.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @export
#' @author Charlotte Soneson
#' @examples
#' data(cobradata_example)
#' cobraperf <- calculate_performance(cobradata_example,
#'                                    binary_truth = "status", aspects = "tpr")
#' cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2",
#'                                    incltruth = TRUE)
#' plot_tpr(cobraplot)
plot_tpr <- function(cobraplot, title = "", stripsize = 15, titlecol = "black",
                     pointsize = 5, xaxisrange = c(0, 1)) {
  plot_fpr_tpr(cobraplot = cobraplot, title = title, stripsize = stripsize,
               titlecol = titlecol, pointsize = pointsize,
               xaxisrange = xaxisrange, aspc = "TPR")
}

#' Plot correlations
#'
#' Plot correlations between observations and a continuous truth value.
#'
#' @param cobraplot A \code{COBRAPlot} object.
#' @param title A character string giving the title of the plot.
#' @param stripsize A numeric value giving the size of the strip text, when the
#'   results are stratified by an annotation.
#' @param titlecol A character string giving the color of the title.
#' @param pointsize A numeric value giving the size of the plot characters.
#' @param xaxisrange A numeric vector with two elements, giving the lower and
#'   upper boundary of the x-axis, respectively.
#' @param corrtype A character string giving the type of correlation to show.
#'   Either "pearson" or "spearman".
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @export
#' @author Charlotte Soneson
#' @examples
#' data(cobradata_example)
#' cobraperf <- calculate_performance(cobradata_example, cont_truth = "logFC",
#'                                    aspects = "corr")
#' cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2",
#'                                    incltruth = TRUE)
#' plot_corr(cobraplot, corrtype = "spearman")
plot_corr <- function(cobraplot, title = "", stripsize = 15, titlecol = "black",
                      pointsize = 5, xaxisrange = c(-1, 1),
                      corrtype = "pearson") {
  plot_fpr_tpr(cobraplot = cobraplot, title = title, stripsize = stripsize,
               titlecol = titlecol, pointsize = pointsize,
               xaxisrange = xaxisrange, aspc = toupper(corrtype))
}

## ---------------------- ROC or FPC --------------------------------- ##
#' @import ggplot2
plot_roc_fpc <- function(cobraplot, title, stripsize, titlecol, xaxisrange,
                         yaxisrange, maxnfdc, aspc, linewidth) {
  if (aspc == "roc")
    plot_data <- iCOBRA::roc(cobraplot)
  else if (aspc == "fpc")
    plot_data <- iCOBRA::fpc(cobraplot)
  if (!(isTRUE(facetted(cobraplot)))) {
    plot_data$method <- plot_data$fullmethod
  }
  
  ## Number of colors/linetypes
  nlevs <- length(unique(plot_data$method))
  
  pp <- ggplot(plot_data, aes_string(x = ifelse(aspc == "roc", "FPR", "topN"),
                                     y = ifelse(aspc == "roc", "TPR", "FP"),
                                     group = "method", colour = "method")) +
    geom_path(size = linewidth, aes_string(linetype = "method")) +
    scale_linetype_manual(values = rep("solid", nlevs), guide = "none") + 
    scale_color_manual(values = plotcolors(cobraplot), name = "", limits = force) +
    plot_theme(stripsize = stripsize, titlecol = titlecol) +
    xlim(ifelse(aspc == "roc", xaxisrange[1], 0),
         ifelse(aspc == "roc", xaxisrange[2], maxnfdc)) +
    ylim(ifelse(aspc == "roc", yaxisrange[1], 0),
         ifelse(aspc == "roc", yaxisrange[2],
                max(plot_data$FP[plot_data$topN <= maxnfdc]))) +
    ggtitle(title)
  if (isTRUE(facetted(cobraplot))) {
    if (!is.finite(maxsplit(cobraplot)))
      msp <- length(unique(plot_data$splitval))
    else
      msp <- maxsplit(cobraplot)
    pp + facet_wrap(~ splitval, nrow = ceiling((msp + 1)/3))
  } else {
    pp
  }
}

#' Plot ROC curves
#'
#' Plot receiver operating characteristics (ROC) curves.
#'
#' @param cobraplot A \code{COBRAPlot} object.
#' @param title A character string giving the title of the plot.
#' @param stripsize A numeric value giving the size of the strip text, when the
#'   results are stratified by an annotation.
#' @param titlecol A character string giving the color of the title.
#' @param xaxisrange A numeric vector with two elements, giving the lower and
#'   upper boundary of the x-axis, respectively.
#' @param yaxisrange A numeric vector with two elements, giving the lower and
#'   upper boundary of the y-axis, respectively.
#' @param linewidth The line width used for plotting
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @export
#' @author Charlotte Soneson
#' @examples
#' data(cobradata_example)
#' cobraperf <- calculate_performance(cobradata_example,
#'                                    binary_truth = "status", aspects = "roc")
#' cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2",
#'                                    incltruth = TRUE)
#' plot_roc(cobraplot)
plot_roc <- function(cobraplot, title = "", stripsize = 15, titlecol = "black",
                     xaxisrange = c(0, 1), yaxisrange = c(0, 1), 
                     linewidth = 1) {
  plot_roc_fpc(cobraplot = cobraplot, title = title, stripsize = stripsize,
               titlecol = titlecol, xaxisrange = xaxisrange,
               yaxisrange = yaxisrange, maxnfdc = NULL, aspc = "roc",
               linewidth = linewidth)
}

#' Plot FP curves
#'
#' Plot false positive curves, indicating the number of false positives among
#' the top-ranked N variables, for varying values of N.
#'
#' @param cobraplot A \code{COBRAPlot} object.
#' @param title A character string giving the title of the plot.
#' @param stripsize A numeric value giving the size of the strip text, when the
#'   results are stratified by an annotation.
#' @param titlecol A character string giving the color of the title.
#' @param maxnfdc A numeric value giving the largest N to consider.
#' @param linewidth The line width used for plotting
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @export
#' @author Charlotte Soneson
#' @examples
#' data(cobradata_example)
#' cobraperf <- calculate_performance(cobradata_example,
#'                                    binary_truth = "status", aspects = "fpc")
#' cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2",
#'                                    incltruth = TRUE)
#' plot_fpc(cobraplot, maxnfdc = 750)
plot_fpc <- function(cobraplot, title = "", stripsize = 15, titlecol = "black",
                     maxnfdc = 500, linewidth = 1) {
  plot_roc_fpc(cobraplot = cobraplot, title = title, stripsize = stripsize,
               titlecol = titlecol, xaxisrange = NULL, yaxisrange = NULL,
               maxnfdc = maxnfdc, aspc = "fpc", linewidth = linewidth)
}

## ------------------------- SCATTER --------------------------------- ##
#' Plot scatter plots
#'
#' Plot scatter plots, indicating the relationship between observed values and a
#' continuous truth.
#'
#' @param cobraplot A \code{COBRAPlot} object.
#' @param title A character string giving the title of the plot.
#' @param stripsize A numeric value giving the size of the strip text, when the
#'   results are stratified by an annotation.
#' @param titlecol A character string giving the color of the title.
#' @param pointsize A numeric value giving the size of the plot characters.
#' @param doflip A logical indicating whether to flip the axes when results are
#'   stratified by an annotation. By default (\code{doflip = FALSE}),
#'   stratification levels are shown as columns and methods as rows in the plot.
#' @param dolog A logical indicating whether to log10-transform values before
#'   plotting.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @export
#' @author Charlotte Soneson
#' @examples
#' data(cobradata_example)
#' cobraperf <- calculate_performance(cobradata_example, cont_truth = "logFC",
#'                                    aspects = "scatter")
#' cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2",
#'                                    incltruth = TRUE)
#' plot_scatter(cobraplot)
plot_scatter <- function(cobraplot, title = "", stripsize = 10,
                         titlecol = "black", pointsize = 3, doflip = FALSE,
                         dolog = FALSE) {
  plot_data <- scatter(cobraplot)
  
  if (isTRUE(dolog)) {
    plot_data$OBSERVATION[which(plot_data$OBSERVATION < 1e-305)] <- 0
    plot_data$TRUTH[which(plot_data$TRUTH < 1e-305)] <- 0
  }
  
  if (isTRUE(facetted(cobraplot))) {
    plot_data$fullmethod <- plot_data$method
  }
  pp <- ggplot(plot_data, aes_string(x = "OBSERVATION", y = "TRUTH",
                                     colour = "fullmethod")) +
    geom_point(size = pointsize) +
    scale_color_manual(values = plotcolors(cobraplot), name = "", limits = force) +
    plot_theme(stripsize = stripsize, titlecol = titlecol) +
    ggtitle(title)
  if (isTRUE(facetted(cobraplot))) {
    if (isTRUE(doflip))
      pp <- pp + facet_grid(splitval ~ method)
    else
      pp <- pp + facet_grid(method ~ splitval)
  } else {
    pp <- pp + facet_wrap(~ method)
  }
  if (isTRUE(dolog))
    pp <- pp + scale_x_log10() + scale_y_log10()
  pp
}


## ------------------- FDRTPR, FDRNBR or FSRNBR ---------------------- ##
#' @import ggplot2
plot_fdrcurve <- function(cobraplot, title, stripsize, titlecol, pointsize,
                          xaxisrange, yaxisrange, plottype, aspc, linewidth) {
  if (aspc == "TPR") {
    plot_data_lines <- fdrtprcurve(cobraplot)
    plot_data_points <- fdrtpr(cobraplot)
    xasp <- "FDR"
    yasp <- aspc
  } else if (aspc == "NBR") {
    plot_data_lines <- fdrnbrcurve(cobraplot)
    plot_data_points <- fdrnbr(cobraplot)
    xasp <- "FDR"
    yasp <- aspc
  } else if (aspc == "FSR") {
    plot_data_lines <- fsrnbrcurve(cobraplot)
    plot_data_points <- fsrnbr(cobraplot)
    xasp <- "FSR"
    yasp <- "NBR"
  }
  
  if ("curve" %in% plottype && length(plot_data_lines) == 0) {
    message(paste0("The provided 'plottype' argument includes 'curve' but ", 
                   "the required values are not calculated. Please include ", 
                   "the appropriate aspect in calculate_performance()."))
    plottype <- setdiff(plottype, "curve")
  }
  if ("points" %in% plottype && length(plot_data_points) == 0) {
    message(paste0("The provided 'plottype' argument includes 'points' but ", 
                   "the required values are not calculated. Please include ", 
                   "the appropriate aspect in calculate_performance()."))
    plottype <- setdiff(plottype, "points")
  }
  pp <- ggplot()
  
  thresholds <- sort(unique(as.numeric(gsub("thr", "", plot_data_points$thr))))
  plot_data_points$method2.satis <- paste0(plot_data_points$method,
                                           plot_data_points$satis)
  
  if (!(isTRUE(facetted(cobraplot)))) {
    plot_data_points$method <- plot_data_points$fullmethod
    plot_data_lines$method <- plot_data_lines$fullmethod
    plot_data_points$method2.satis <- paste0(plot_data_points$fullmethod,
                                             plot_data_points$satis)
  }
  
  ## Number of colors/linetypes
  nlevs <- length(unique(plot_data_lines$method))
  ## Number of thresholds
  nthr <- length(unique(plot_data_points$thr))
  
  if ("curve" %in% plottype && "points" %in% plottype) {
    pp <- ggplot(plot_data_lines, aes_string(x = xasp, y = yasp,
                                             group = "method",
                                             colour = "method")) +
      geom_vline(xintercept = seq(0, xaxisrange[2], 0.1),
                 colour = "lightgrey", linetype = "dashed") +
      geom_vline(xintercept = thresholds, linetype = "dashed") +
      geom_path(size = linewidth, aes_string(linetype = "method")) +
      scale_linetype_manual(values = rep("solid", nlevs), guide = "none") + 
      geom_point(data = plot_data_points, size = pointsize,
                 aes_string(fill = "method2.satis", colour = "method", 
                            shape = "thr"),
                 stroke = 1) +
      scale_shape_manual(values = rep(21, nthr), guide = "none") + 
      scale_fill_manual(values = plotcolors(cobraplot), guide = "none",
                        name = "") +
      scale_color_manual(values = plotcolors(cobraplot), name = "", limits = force) +
      ylim(ifelse(yasp == "TPR", yaxisrange[1], 0),
           ifelse(yasp == "TPR", yaxisrange[2],
                  max(c(0, plot_data_lines[[yasp]][plot_data_lines[[xasp]] <=
                                                     xaxisrange[2]])))) +
      scale_x_continuous(breaks = c(thresholds,
                                    seq(0, xaxisrange[2], 0.1)),
                         labels = c(thresholds, "", 
                                    seq(0, xaxisrange[2], 0.1)[-1]),
                         limits = c(xaxisrange[1], xaxisrange[2])) +
      plot_theme(stripsize = stripsize, titlecol = titlecol) +
      ggtitle(title)
    
  } else if ("curve" %in% plottype) {
    pp <- ggplot(plot_data_lines,
                 aes_string(x = xasp, y = yasp,
                            group = "method", colour = "method")) +
      geom_path(size = linewidth, aes_string(linetype = "method")) +
      scale_linetype_manual(values = rep("solid", nlevs), guide = "none") + 
      xlim(xaxisrange[1], xaxisrange[2]) +
      ylim(ifelse(yasp == "TPR", yaxisrange[1], 0),
           ifelse(yasp == "TPR", yaxisrange[2],
                  max(c(0, plot_data_lines[[yasp]][plot_data_lines[[xasp]] <=
                                                     xaxisrange[2]])))) +
      scale_color_manual(values = plotcolors(cobraplot), name = "", limits = force) +
      plot_theme(stripsize = stripsize, titlecol = titlecol) +
      ggtitle(title)
  } else if ("points" %in% plottype) {
    pp <- ggplot(plot_data_points, aes_string(x = xasp, y = yasp,
                                              group = "method")) +
      geom_vline(xintercept = seq(0, xaxisrange[2], 0.1),
                 colour = "lightgrey", linetype = "dashed") +
      geom_vline(xintercept = thresholds, linetype = "dashed") +
      geom_path(size = linewidth, 
                aes_string(colour = "method", linetype = "method")) +
      scale_linetype_manual(values = rep("solid", nlevs), guide = "none") + 
      geom_point(size = pointsize,
                 aes_string(fill = "method2.satis", colour = "method",
                            shape = "thr"),
                 stroke = 1) +
      scale_shape_manual(values = rep(21, nthr), guide = "none") + 
      scale_fill_manual(values = plotcolors(cobraplot), guide = "none",
                        name = "") +
      scale_color_manual(values = plotcolors(cobraplot), name = "", limits = force) +
      ylim(ifelse(yasp == "TPR", yaxisrange[1], 0),
           ifelse(yasp == "TPR", yaxisrange[2],
                  max(c(0, plot_data_lines[[yasp]][plot_data_lines[[xasp]] <=
                                                     xaxisrange[2]])))) +
      scale_x_continuous(breaks = c(thresholds,
                                    seq(0, xaxisrange[2], 0.1)),
                         labels = c(thresholds, "", 
                                    seq(0, xaxisrange[2], 0.1)[-1]),
                         limits = c(xaxisrange[1], xaxisrange[2])) +
      plot_theme(stripsize = stripsize, titlecol = titlecol) +
      ggtitle(title)
  }
  if (isTRUE(facetted(cobraplot))) {
    if (!is.finite(maxsplit(cobraplot))) {
      if (length(plot_data_lines) != 0)
        msp <- length(unique(plot_data_lines$splitval))
      else
        msp <- length(unique(plot_data_points$splitval))
    } else {
      msp <- maxsplit(cobraplot)
    }
    pp + facet_wrap(~ splitval, nrow = ceiling((msp + 1)/3))
  } else {
    pp
  }
}

#' Plot TPR vs FDR
#'
#' Plot observed true positive rate (TPR) vs observed false discovery rate
#' (FDR), for given adjusted p-value thresholds and/or as curves traced out by
#' considering all threshold values.
#'
#' @param cobraplot A \code{COBRAPlot} object.
#' @param title A character string giving the title of the plot.
#' @param stripsize A numeric value giving the size of the strip text, when the
#'   results are stratified by an annotation.
#' @param titlecol A character string giving the color of the title.
#' @param pointsize A numeric value giving the size of the plot characters.
#' @param xaxisrange A numeric vector with two elements, giving the lower and
#'   upper boundary of the x-axis, respectively.
#' @param yaxisrange A numeric vector with two elements, giving the lower and
#'   upper boundary of the y-axis, respectively.
#' @param plottype A character vector giving the type of plot to construct. Can
#'   be any combination of the two elements "curve" and "points".
#' @param linewidth The line width used for plotting
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @export
#' @author Charlotte Soneson
#' @examples
#' data(cobradata_example)
#' cobraperf <- calculate_performance(cobradata_example,
#'                                    binary_truth = "status",
#'                                    aspects = c("fdrtpr", "fdrtprcurve"))
#' cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2",
#'                                    incltruth = TRUE)
#' plot_fdrtprcurve(cobraplot, plottype = c("curve", "points"))
plot_fdrtprcurve <- function(cobraplot, title = "", stripsize = 15,
                             titlecol = "black", pointsize = 5,
                             xaxisrange = c(0, 1), yaxisrange = c(0, 1),
                             plottype = c("curve", "points"),
                             linewidth = 1) {
  plot_fdrcurve(cobraplot = cobraplot, title = title, stripsize = stripsize,
                titlecol = titlecol, pointsize = pointsize,
                xaxisrange = xaxisrange, yaxisrange = yaxisrange,
                plottype = plottype, aspc = "TPR", linewidth = linewidth)
}

#' Plot number of significant features vs FDR
#'
#' Plot the number of features considered significant vs observed false
#' discovery rate (FDR), for given adjusted p-value thresholds and/or as curves
#' traced out by considering all threshold values.
#'
#' @param cobraplot A \code{COBRAPlot} object.
#' @param title A character string giving the title of the plot.
#' @param stripsize A numeric value giving the size of the strip text, when the
#'   results are stratified by an annotation.
#' @param titlecol A character string giving the color of the title.
#' @param pointsize A numeric value giving the size of the plot characters.
#' @param xaxisrange A numeric vector with two elements, giving the lower and
#'   upper boundary of the x-axis, respectively.
#' @param plottype A character vector giving the type of plot to construct. Can
#'   be any combination of the two elements "curve" and "points".
#' @param linewidth The line width used for plotting
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @export
#' @author Charlotte Soneson
#' @examples
#' data(cobradata_example)
#' cobraperf <- calculate_performance(cobradata_example,
#'                                    binary_truth = "status",
#'                                    aspects = c("fdrnbr", "fdrnbrcurve"))
#' cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2",
#'                                    incltruth = TRUE)
#' plot_fdrnbrcurve(cobraplot, plottype = c("curve", "points"))
plot_fdrnbrcurve <- function(cobraplot, title = "", stripsize = 15,
                             titlecol = "black", pointsize = 5,
                             xaxisrange = c(0, 1),
                             plottype = c("curve", "points"),
                             linewidth = 1) {
  plot_fdrcurve(cobraplot = cobraplot, title = title, stripsize = stripsize,
                titlecol = titlecol, pointsize = pointsize,
                xaxisrange = xaxisrange, yaxisrange = NULL,
                plottype = plottype, aspc = "NBR",
                linewidth = linewidth)
}

#' Plot number of features with s-value below threshold vs FSR
#'
#' Plot the number of features with an s-value below a threshold vs the observed
#' false sign rate (FSR), for given adjusted p-value thresholds and/or as curves
#' traced out by considering all threshold values.
#'
#' @param cobraplot A \code{COBRAPlot} object.
#' @param title A character string giving the title of the plot.
#' @param stripsize A numeric value giving the size of the strip text, when the
#'   results are stratified by an annotation.
#' @param titlecol A character string giving the color of the title.
#' @param pointsize A numeric value giving the size of the plot characters.
#' @param xaxisrange A numeric vector with two elements, giving the lower and
#'   upper boundary of the x-axis, respectively.
#' @param plottype A character vector giving the type of plot to construct. Can
#'   be any combination of the two elements "curve" and "points".
#' @param linewidth The line width used for plotting
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @export
#' @author Charlotte Soneson
#' @examples
#' data(cobradata_example_sval)
#' cobraperf <- calculate_performance(cobradata_example_sval,
#'                                    cont_truth = "logFC",
#'                                    aspects = c("fsrnbr", "fsrnbrcurve"))
#' cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2",
#'                                    incltruth = TRUE)
#' plot_fsrnbrcurve(cobraplot, plottype = c("curve", "points"))
plot_fsrnbrcurve <- function(cobraplot, title = "", stripsize = 15,
                             titlecol = "black", pointsize = 5,
                             xaxisrange = c(0, 1),
                             plottype = c("curve", "points"),
                             linewidth = 1) {
  plot_fdrcurve(cobraplot = cobraplot, title = title, stripsize = stripsize,
                titlecol = titlecol, pointsize = pointsize,
                xaxisrange = xaxisrange, yaxisrange = NULL,
                plottype = plottype, aspc = "FSR",
                linewidth = linewidth)
}

## ------------------------ OVERLAP ---------------------------------- ##
#' Plot Venn diagram
#'
#' Plot a Venn diagram showing the overlaps among sets of significant feature
#' for a given adjusted p-value threshold. Optionally, the truth can be included
#' as a "perfect" method. Note that maximally five methods (including the truth,
#' if applicable) can be compared.
#'
#' @param cobraplot A \code{COBRAPlot} object.
#' @param ... Additional arguments to \code{limma::vennDiagram}.
#'
#' @return Nothing, displays a graph
#'
#' @export
#' @import limma
#' @author Charlotte Soneson
#' @examples
#' data(cobradata_example)
#' cobraperf <- calculate_performance(cobradata_example,
#'                                    binary_truth = "status",
#'                                    aspects = "overlap")
#' cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2",
#'                                    incltruth = TRUE)
#' plot_overlap(cobraplot)
plot_overlap <- function(cobraplot, ...) {
  overlap_table <- overlap(cobraplot)
  circle.col <- plotcolors(cobraplot)
  if (length(overlap_table) == 0)
    return(NULL)
  
  if (!is(overlap_table, "list")) {
    if (ncol(overlap_table) < 6) {
      if (ncol(overlap_table) == 5)
        cols <- rep(circle.col[colnames(overlap_table)], 2)[2:6]
      else
        cols <- circle.col[colnames(overlap_table)]
      vennDiagram(overlap_table, circle.col = cols, ...)
    } else {
      NULL
    }
  } else {
    ncl <- ceiling(sqrt(length(overlap_table)))
    nrw <- ceiling(length(overlap_table)/ncl)
    graphics::par(mfrow = c(nrw, ncl), mar = c(4, 1, 1, 1))
    for (i in seq_len(length(overlap_table))) {
      if (ncol(overlap_table[[i]]) < 6) {
        if (ncol(overlap_table[[i]]) == 5)
          cols <- rep(circle.col[colnames(overlap_table[[i]])], 2)[2:6]
        else
          cols <- circle.col[colnames(overlap_table[[i]])]
        vennDiagram(overlap_table[[i]], circle.col = cols,
                    main = paste0(splv(cobraplot), ":",
                                  names(overlap_table)[i]), ...)
      } else {
        NULL
      }
    }
    graphics::par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
  }
}

#' Create UpSet plots
#'
#' Generate UpSet plots showing the overlaps among sets of significant feature
#' for a given adjusted p-value threshold. Optionally, the truth can be included
#' as a "perfect" method. Note that if the results are stratified, only one 
#' category at a time can be displayed.
#'
#' @param cobraplot A \code{COBRAPlot} object.
#' @param stratum If results are stratified, the category to plot results for.
#'   Can be numeric or categorical (the name of the category).
#' @param nsets The number of methods to include. By default, it is determined
#'   automatically from the \code{cobraplot} object.
#' @param nintersects The number of set intersections to display. By default, it
#'   is determined automatically from the \code{cobraplot} object.
#' @param sets.bar.color The colors to use for the bars in the UpSet plot. By
#'   default, they are extracted from the \code{plotcolors} slot of the
#'   \code{cobraplot} object.
#' @param ... Additional arguments to \code{UpSetR::upset}.
#'
#' @return Nothing, displays a graph
#' @references 
#' Lex and Gehlenborg (2014): Points of view: Sets and intersections. Nature
#' Methods 11, 779.
#' 
#' Lex et al (2014): UpSet: Visualization of intersecting sets. IEEE
#' Transactions on Visualization and Computer Graphics 20(12), 1983-1992.
#' @export
#' @import UpSetR
#' @author Charlotte Soneson
#' @examples
#' data(cobradata_example)
#' cobraperf <- calculate_performance(cobradata_example,
#'                                    binary_truth = "status",
#'                                    aspects = "overlap")
#' cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2",
#'                                    incltruth = TRUE)
#' plot_upset(cobraplot)
#' plot_upset(cobraplot, order.by = "freq", decreasing = TRUE)
#' 
#' cobraperf <- calculate_performance(cobradata_example, 
#'                                    binary_truth = "status", 
#'                                    aspects = "overlap",
#'                                    splv = "expr_cat")
#' cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2", 
#'                                    incltruth = TRUE)
#' plot_upset(cobraplot, stratum = "[2.85e+00,1.45e+01)")
plot_upset <- function(cobraplot, stratum = NULL, nsets = NULL, 
                       nintersects = NULL, sets.bar.color = NULL, ...) {
  overlap_table <- overlap(cobraplot)
  if (length(overlap_table) == 0)
    return(NULL)
  
  if (!is(overlap_table, "list")) {
    plotorder <- 
      colnames(overlap_table)[order(colSums(overlap_table), 
                                    seq(seq_len(ncol(overlap_table))),
                                    decreasing = "true")]
    if (all(colSums(overlap_table) == 0)) return(NULL)
    if (is.null(nsets)) nsets <- ncol(overlap_table)
    if (is.null(nintersects)) nintersects <- 2^(ncol(overlap_table)) - 1
    if (is.null(sets.bar.color)) 
      sets.bar.color <- plotcolors(cobraplot)[plotorder]
    upset(overlap_table, nsets = nsets, nintersects = nintersects, 
          sets.bar.color = sets.bar.color, ...)
  } else {
    if (is.null(stratum)) stop("You must provide a stratum")
    plotorder <- 
      colnames(overlap_table[[stratum]])[
        order(colSums(overlap_table[[stratum]]), 
              seq(seq_len(ncol(overlap_table[[stratum]]))),
              decreasing = "true")]
    if (all(colSums(overlap_table[[stratum]]) == 0)) return(NULL)
    if (is.null(nsets)) nsets <- ncol(overlap_table[[stratum]])
    if (is.null(nintersects)) 
      nintersects <- 2^(ncol(overlap_table[[stratum]])) - 1
    if (is.null(sets.bar.color)) 
      sets.bar.color <- plotcolors(cobraplot)[plotorder]
    upset(overlap_table[[stratum]], nsets = nsets, nintersects = nintersects,
          sets.bar.color = sets.bar.color, ...)
  }
}

## -------------------------- Deviation ------------------------------ ##
#' Plot deviations
#'
#' Plot the deviations between observed scores and the continuous truth
#' variable.
#'
#' @param cobraplot A \code{COBRAPlot} object.
#' @param title A character string giving the title of the plot.
#' @param stripsize A numeric value giving the size of the strip text, when the
#'   results are stratified by an annotation.
#' @param titlecol A character string giving the color of the title.
#' @param xaxisrange A numeric vector with two elements, giving the lower and
#'   upper boundary of the x-axis, respectively.
#' @param plottype Either "boxplot" or "violinplot", indicating what type of
#'   plot to make.
#' @param dojitter A logical indicating whether to include jittered data points
#'   or not.
#' @param transf A character indicating the transformation to apply to the
#'   deviations before plotting. Must be one of "raw", "absolute" or "squared"
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @export
#' @author Charlotte Soneson
#' @examples
#' data(cobradata_example)
#' cobraperf <- calculate_performance(cobradata_example, cont_truth = "logFC",
#'                                    aspects = "deviation")
#' cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2",
#'                                    incltruth = TRUE)
#' plot_deviation(cobraplot)
plot_deviation <- function(cobraplot, title = "", stripsize = 15,
                           titlecol = "black", xaxisrange = NULL,
                           plottype = "boxplot",
                           dojitter = TRUE, transf = "raw") {
  stopifnot(transf %in% c("raw", "absolute", "squared"))
  plot_data <- deviation(cobraplot)
  if (transf == "absolute")
    plot_data$absDEVIATION <- abs(plot_data$DEVIATION)
  else if (transf == "squared")
    plot_data$sqDEVIATION <- plot_data$DEVIATION^2
  
  if (!(isTRUE(facetted(cobraplot)))) {
    plot_data$method <- plot_data$fullmethod
  }
  
  pp <- ggplot(
    plot_data, aes_string(x = "method",
                          y = ifelse(transf == "raw", "DEVIATION",
                                     ifelse(transf == "absolute",
                                            "absDEVIATION", "sqDEVIATION")),
                          group = "method", colour = "method")) +
    coord_flip() +
    scale_color_manual(values = plotcolors(cobraplot), name = "", limits = force) +
    plot_theme(stripsize = stripsize, titlecol = titlecol) +
    ggtitle(title)
  if (plottype == "boxplot") {
    if (isTRUE(dojitter))
      pp <- pp + geom_boxplot(outlier.size = 0)
    else
      pp <- pp + geom_boxplot()
  }
  else if (plottype == "violinplot")
    pp <- pp + geom_violin()
  if (isTRUE(dojitter))
    pp <- pp + geom_jitter(position = position_jitter(width = 0.1, height = 0),
                           size = 1.5)
  if (isTRUE(facetted(cobraplot))) {
    if (!is.finite(maxsplit(cobraplot)))
      msp <- length(unique(plot_data$splitval))
    else
      msp <- maxsplit(cobraplot)
    pp <- pp + facet_wrap(~ splitval,
                          nrow = ceiling((msp + 1)/3))
  }
  if (!is.null(xaxisrange))
    pp <- pp + ylim(xaxisrange[1], xaxisrange[2])
  pp
}

#' @inheritParams sce_image_clus
#' @param d A data.frame with the sample-level information. This is typically
#' obtained using `as.data.frame(colData(sce))`.
#' @param title The title for the plot.
#'
#' @return A [ggplot2][ggplot2::ggplot] object.
#' @export
#' @importFrom S4Vectors metadata
#' @family Spatial cluster visualization functions
#'
#' @examples
#'
#' ## Obtain the necessary data
#' if (!exists('ori_sce')) ori_sce <- fetch_data('sce')
#' sce_sub <- ori_sce[, ori_sce$sample_name == '151673']
#'
#' ## Use the manual color palette by Lukas M Weber
#' ## Don't plot the histology information
#' sce_image_clus_p(
#'     sce = sce_sub,
#'     d = as.data.frame(colData(sce_sub)),
#'     clustervar = 'layer_guess_reordered',
#'     sampleid = '151673',
#'     colors = libd_layer_colors,
#'     title = '151673 LIBD Layers',
#'     spatial = FALSE
#' )
#'

sce_image <-
  function(sce,
           #d,
           clustervar,
           coordinates = c("col","row"),
           path = '~',
           save = T,
           width = 6, 
           height = 4,
           #sampleid,
           #colors,
           #spatial,
           method = "",
           title = NULL) {
    
    ## Some variables
    key <- NULL
    d <- as.data.frame(colData(sce))
    n_clus <- length(unique(factor(d[,clustervar])))
    if (n_clus > 12)
    {colors <-Polychrome::palette36.colors(n_clus)}
    else
    {colors <-c(
      "#A6CEE3","#B2DF8A",
               "#FB9A99","#FDBF6F",
               "#CAB2D6","#FFFF99",
               "#1F78B4","#33A02C","#E31A1C","#FF7F00","#6A3D9A",
               "#B15928"
    )}
    
    p <- ggplot(d,
                aes(
                  x = as.numeric(as.character(d[, coordinates[1]])),
                  y = as.numeric(as.character(d[, coordinates[2]])),
                  fill = factor(!!sym(clustervar)),
                  key =  key
                ))
    
    p <- p +
      geom_point(shape = 21,
                 size = 1,
                 stroke = 0.01) +
      coord_cartesian(expand = TRUE) +
      scale_fill_manual(values = colors) +
      #xlim(0, max(sce$width)) +
      #ylim(max(sce$height), 0) +
      xlab("") + ylab("") +
      labs(fill = NULL) +
      guides(fill = guide_legend(override.aes = list(size = 3))) +
      #ggtitle(title) +
      #theme_set(theme_bw(base_size = 10)) +
      #theme(
      #    panel.grid.major = element_blank(),
      #    panel.grid.minor = element_blank(),
      #    panel.background = element_blank(),
      #    axis.line = element_blank(),
      #    axis.text = element_blank(),
      #    axis.ticks = element_blank()
      #)
      theme_bw() + 
      theme(strip.text = element_text(face="bold", size=6,lineheight=5.0), 
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(), 
            axis.text.y = element_blank(), 
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()
      )
    
    if (save == T){
      filename <- paste0(path,"expression_plots_",method, ".png")
      p
      ggsave(filename, width = width, height = height)
    }
    return(p)
  }
