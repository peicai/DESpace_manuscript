rm(list = ls())

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# aggregate all melanoma results:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
path = "./Real_data/results/"
path_save = "./Figures/Figures/Supplementary/"
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Functions.R")
sample_names = c("mel1_rep1", "mel1_rep2",
                 "mel2_rep1", "mel2_rep2",
                 "mel3_rep1", "mel3_rep2",
                 "mel4_rep1", "mel4_rep2")
list_result_melanoma <- list()
list_result_melanoma <- process_data(
   path = path,
   data = 'melanoma', Manual = FALSE,
   threshold = 2, # indicates that we do not set threshold
   name_head = -13,name_tail = -5,
   sig_genes = FALSE
)
names(list_result_melanoma) <- c("BayesSpace_DESpace", "StLearn_DESpace",
                              "BayesSpace_scran", "StLearn_scran",
                              "BayesSpace_seurat", "StLearn_seurat",
                              "MERINGUE", "SPARK","nnSVG", "SPARK-X", "SpatialDE2","SpatialDE","SpaGCN",
                              "All_results")
save(list_result_melanoma, file = "./Real_data/results/melanoma/All_methods_results.rda")

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# load results:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
load('./Real_data/results/melanoma/All_melanoma_results.rda')
length(data)

methods = sort(colnames(data[[1]]))

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# load associated genes
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
BUZZ_WORDS = list.files("melanoma")
BUZZ_WORDS = substring(BUZZ_WORDS, 1,  nchar(BUZZ_WORDS) - 5)

library(rjson)
library(doRNG)
GENES <- foreach(i =  1:length(BUZZ_WORDS), .combine = "c") %do% {
  a = unlist(fromJSON(file = paste0("melanoma/", BUZZ_WORDS[i],".json")))
  sel = names(a) == "Ensembl"
  list(a[sel])
}
str(GENES)
names(GENES) <- BUZZ_WORDS

sapply(GENES, length)

# keep HLA genes only:
a = unlist(fromJSON(file = paste0("melanoma/HLA.json")))
sel = names(a) == "Gene"
a = a[sel]
sel_HLA = grep("HLA", a)
a[sel_HLA]
rm(a); rm(sel)

GENES$HLA = GENES$HLA[1:19]

GENES$overall = unique(unlist(GENES))

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# loop over 12 samples:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
FOUND = list()
for(i in 1:length(data)){
  results = data[[i]]
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # top results:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  results[is.na(results)] = 1
  
  n_top = 500
  
  res = sapply(methods, function(method){
    sel = which(colnames(results) == method )
    ordering = order(results[,sel])
    genes = unique(rownames(results)[ ordering ])
    genes[1:n_top]
  })
  head(res)
  
  found = apply(res[,1:length(methods)],2, function(xx){
    sapply(GENES, function(gg){
      sum(xx %in% gg)
    })
  })
  found
  
  gene_names = unique(rownames(results))
  n_genes = length(gene_names)
  expected = sapply(GENES, function(gene){
    sum(gene_names %in% gene) * n_top/n_genes
  })
  #found = cbind(found, expected)
  
  FOUND[[i]] = found
}

str(FOUND)

res = do.call(rbind, FOUND)

rows = unique(rownames(res))

a = lapply(rows, function(row){
  sel_res = rownames(res) == row
  colSums(res[sel_res,])
  # sort(colSums(res[sel_res,]), decreasing = TRUE)
})
names(a) = rows
a

A = do.call(rbind, a)
A = A[,order(A[4,], decreasing = TRUE)]
A

library(xtable)

xtable(t(A), digits = 0)
\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrrrrrrrrr}
\hline
& BayesSpace\_edgeR & StLearn\_edgeR & BayesSpace\_seurat & StLearn\_seurat & SPARKX & StLearn\_scran & nnSVG & SPARK & SpatialDE2 & MERINGUE & SpatialDE & BayesSpace\_scran & SpaGCN \\ 
\hline
HLA & 65 & 58 & 40 & 50 & 48 & 57 & 64 & 63 & 50 & 64 & 56 & 50 & 36 \\ 
melanoma marker & 8 & 7 & 7 & 5 & 7 & 6 & 6 & 7 & 5 & 8 & 7 & 6 & 5 \\ 
melanoma & 1073 & 1062 & 1069 & 1049 & 1042 & 1035 & 1021 & 1011 & 1006 & 995 & 984 & 961 & 913 \\ 
overall & 1107 & 1094 & 1092 & 1078 & 1068 & 1064 & 1055 & 1044 & 1030 & 1029 & 1014 & 986 & 936 \\ 
\hline
\end{tabular}
\end{table}
