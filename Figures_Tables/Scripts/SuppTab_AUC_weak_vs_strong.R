## SuppTab 3: AUC - individual samples weak vs. strong
rm(list = ls())
setwd("~/DESpace_manuscript")
path = "./Simulation/Output/"
path_save = "./Figures/Figures/Supplementary/"
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Functions.R")
# remove the first "Manual_DESpace" color and method
colours = c(all_colours[2:14], "white")
methods_all = methods_all[2:14]
## For LIBD
pval_auc_LIBD = padj_auc_LIBD = list()
sample_names = c("151507", "151669", "151673")
patterns = c("bottom_patch","circle_patch",
             "Manual_clusters_patch",
             "mixture_patch","mixture_reverse_patch")

gg <- list()
setwd("./DESpace_data")
new_data_path <- "./LIBD/"
k=1

for(i in seq_along(patterns)){
  sample_names = c("151507", "151669", "151673")
  if(i %in% c(4,5)) spatial_probs = c(0.6,0.6)
  if(i %in% c(1,2,3)) spatial_probs = c(0.6,0.9)
  if(i == 1) sample_names <- sample_names[-3] # remove the sample 151673

  gg = overall_roc_fdr_plot(
    sample_names = sample_names, 
    pattern_names = patterns[i],
    path = path,colours = colours,new_data_path = new_data_path,
    spatial_probs = spatial_probs,default = FALSE,rep_i = 1,
    methods_all = methods_all,
    methods_order = methods_order,
    shape_border = shape_border,
    shape_fill = shape_fill,
    auc = TRUE)
  
  pval_auc_LIBD[[i]] = gg[[1]]
}

table <- do.call(rbind, pval_auc_LIBD)
tab <- t(table) %>% as.data.frame()
tab$Average <- rowMeans(tab)
tab <- tab[order(tab$Average, decreasing = TRUE),]
tab <- round(tab,2)
colnames(tab) <- c("Bottom/Right", "Circular", "Annotations", "Mixture", "Inverted mixture", "Average")
tab$Mean <- round(rowMeans(tab[,1:5]),2)
write.table(tab, sep = ",", file = paste0('LIBD_pval_simple_auc.csv'),
            row.names = TRUE)

## For 151673; bottom pattern: without 'SPARK'
# i <- 3; pattern_names <- c("bottom_patch")
# spatial_probs <- c(0.6, 0.9); default <- FALSE
# methods_all <- methods_all[-3] # remove 'SPARK'
# colours <- colours[-3]

## For melanoma
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Functions.R")
# remove the first "Manual_DESpace" color and method
colours = c(all_colours[2:14], "white")
methods_all = methods_all[2:14]

pval_auc_melanoma = padj_auc_melanoma = list()
sample_names = c("mel1_rep1", "mel2_rep1", "mel3_rep1", "mel4_rep1")
patterns = c("right_patch","circle_patch",
             "BayesSpace_clusters_patch",
             "mixture_patch","mixture_reverse_patch")
setwd("./DESpace_data")
new_data_path <- "./melanoma/"
for(i in seq_along(patterns)){
  if(i %in% c(4,5)) spatial_probs = c(0.6,0.6)
  if(i %in% c(1,2,3)) spatial_probs = c(0.6,0.9)

  gg = overall_roc_fdr_plot(
    sample_names = sample_names, 
    pattern_names = patterns[i],
    path = path,colours = colours,new_data_path = new_data_path,
    spatial_probs = spatial_probs,default = FALSE,rep_i = 1,
    methods_all = methods_all,
    methods_order = methods_order,
    shape_border = shape_border,
    shape_fill = shape_fill,
    auc = TRUE)
  pval_auc_melanoma[[i]] = gg[[1]]
}

table <- do.call(rbind, pval_auc_melanoma)
tab <- t(table) %>% as.data.frame()
tab$Average <- rowMeans(tab)
tab <- tab[order(tab$Average, decreasing = TRUE),]
tab <- round(tab,2)
colnames(tab) <- c("Bottom/Right", "Circular", "Annotations", "Mixture", "Inverted mixture", "Average")
tab$Mean <- round(rowMeans(tab[,1:5]),2)

write.table(tab, sep = ",", file = paste0('melanoma_pval_simple_auc.csv'),
            row.names = TRUE)
