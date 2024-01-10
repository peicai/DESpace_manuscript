rm(list = ls())
library(data.table)
library(iCOBRA)
library(S4Vectors)
library(stringr)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(pROC)

# Load spatially variable genes (ground truth)
SVG_list <- read.csv("./Simulation/Output/151507/split_mixture_patch/probs_0.5_0.9FALSE/probs_0.9_SVGs.txt",
                     header = TRUE,sep =  '\t')

################################ DESpace based on stLearn ########################

ll <- lapply(c(2,4,6,8,10,12), function(x){
  path <- "./Simulation/Output/151507/split_mixture_patch_"
  results_dir <- paste0(path,x,"/probs_0.5_0.9FALSE/")
  sce_dir <- paste0(results_dir, "sce_edgeR.rda")
  SVG_dir <- paste0(results_dir, "SV_genes.txt")
  DESpace_dir <- paste0(results_dir, "result_SV_edgeR_stLearn_louvain.rda")
  load(DESpace_dir)
  result <- spatial_gene_results
  genes <- read.csv(SVG_dir,header = TRUE,sep =  '\t')
  if(!all(genes$V1 %in% SVG_list$V1)){
    message("Please check the list of spatial variable genes; 
              it should be matched to the ground truth.")
  }
  l2 <-cbind(result, num_clusters = x)
  l2 <- l2[,c("genes","PValue", "FDR", "method", "num_clusters")]
  return(l2)
}
)
result <- do.call(rbind, ll)
all_genes <- as.data.frame((unique(result$genes)))

all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
genes_id <- apply(SVG_list,1,function(x) substr(x, 10,18))
status <- ifelse(all_genes_id %in% genes_id, 1, 0)
truth <- cbind(all_genes, status)
# check
sum(truth$status) == dim(SVG_list)[1]
if(sum(all_genes_id[length(all_genes_id)] == "_count")>=1){
  truth <- truth[1:dim(truth)[1]-1,]
}
result <- result %>% dplyr::distinct()
setDT(result)
result$PValue <- as.numeric(result$PValue)
result$FDR <- as.numeric(result$FDR)
data <- dcast(result, num_clusters ~ genes,value.var = 'FDR', fun.aggregate=mean)
#data <- dcast(result, method ~ genes,value.var = 'adj.p.value')
df2 = truth[match(colnames(data)[-1], truth$`(unique(result$genes))`),]
df2 <- as.data.frame(df2[,2])
rownames(df2) <- colnames(data)[-1]
colnames(df2) <- 'status'
data2 <- dcast(result, num_clusters ~ genes,value.var = 'PValue', fun.aggregate=mean)
methods_all <- data2$num_clusters
pval = data.frame(t(data2[, -1]))
colnames(pval) <- methods_all
methods_all <- data$num_clusters
padj = data.frame(t(data[, -1]))
colnames(padj) <- methods_all


DF_COBRA <- COBRAData(pval,
                      padj,
                      truth = data.frame(df2))
perf <- calculate_performance(DF_COBRA, binary_truth = "status",
                              thrs = c(0.01, 0.05, 0.1, 0.2))
cbbPalette <- c("#000000", "#FFD92F", "#56B4E9", "#33A02C", "#F4A582", "#C51B7D")

cobra_plot <- prepare_data_for_plot(perf, colorscheme = cbbPalette, incloverall = FALSE,
                                    facetted = TRUE,conditionalfill = FALSE,
                                    keepmethods = c("2","4","6","8","10","12")
)

shape_border = c(0, 1, 2, 5, 6, 3)
shape_fill = c(15, 16, 17, 23, 25, 3)
## Manually match "method" and point "shape"
sel_shape = c(5,6,1,2,3,4)
sel_fill = sel_shape
col <- c(cbbPalette)
cobra_plot@plotcolors <- cobra_plot@plotcolors[1:6]
(gg_stLearn <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                                pointsize = 0, linewidth = 2)+
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.8)) +
    scale_color_manual(values = rep(cbbPalette[c(5,6,1,2,3,4)],4),
                       name = "",
                       breaks=c('2','4','6','8','10','12'))+
    guides(colour = guide_legend(ncol = 6, byrow = FALSE,
                                 override.aes = list(shape = shape_fill,
                                                     fill = cbbPalette[c(5,6,1,2,3,4)]) ) ) +
    geom_point(size = 10, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_border[sel_shape],4), 
    stroke = 2, 
    alpha = 1) + # stroke = line width
    scale_fill_manual(values =rep(cbbPalette,4), 
                      name = "",
                      breaks = c('2','4','6','8','10','12'))+
    geom_point(size = 10, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_fill[sel_fill],4), 
    stroke = 2, alpha = 0.25)+
    labs(title = "StLearn") + 
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          title = element_text(size = rel(3)),
          axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(3)),
          axis.text.y=element_text(size=rel(3)),
          axis.title.y = element_text(size=rel(1)),
          axis.title.x = element_text(size=rel(1)),
          legend.title=element_text(size=rel(1)),
          legend.text=element_text(size=rel(3)),
          legend.key.width=unit(2.5, "cm"),
          aspect.ratio = 1, legend.position="bottom",
          legend.box="vertical", legend.margin=margin()) )

################################ DESpace based on BayesSpace ########################

ll <- lapply(c(2,4,6,8,10,12), function(x){
  path <- "~/Desktop/master_thesis/Data_re_run/simulation/results/151507/split_mixture_patch_"
  results_dir <- paste0(path,x,"/probs_0.5_0.9FALSE/")
  sce_dir <- paste0(results_dir, "sce_edgeR.rda")
  SVG_dir <- paste0(results_dir, "SV_genes.txt")
  DESpace_dir <- paste0(results_dir, "result_SV_edgeR_counts.rda")
  load(DESpace_dir)
  result <- spatial_gene_results
  genes <- read.csv(SVG_dir,header = TRUE,sep =  '\t')
  if(!all(genes$V1 %in% SVG_list$V1)){
    message("Please check the list of spatial variable genes; 
              it should be matched to the ground truth.")
  }
  l2 <-cbind(result, num_clusters = x)
  return(l2)
}
)
result <- do.call(rbind, ll)

all_genes <- as.data.frame((unique(result$genes)))

all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
genes_id <- apply(SVG_list,1,function(x) substr(x, 10,18))
status <- ifelse(all_genes_id %in% genes_id, 1, 0)
truth <- cbind(all_genes, status)
# check
sum(truth$status) == dim(SVG_list)[1]
if(sum(all_genes_id[length(all_genes_id)] == "_count")>=1){
  truth <- truth[1:dim(truth)[1]-1,]
}
result <- result %>% dplyr::distinct()
setDT(result)
result$p.value <- as.numeric(result$p.value)
result$adj.p.value <- as.numeric(result$adj.p.value)
data <- dcast(result, num_clusters ~ genes,value.var = 'adj.p.value', fun.aggregate=mean)
#data <- dcast(result, method ~ genes,value.var = 'adj.p.value')
df2 = truth[match(colnames(data)[-1], truth$`(unique(result$genes))`),]
df2 <- as.data.frame(df2[,2])
rownames(df2) <- colnames(data)[-1]
colnames(df2) <- 'status'
data2 <- dcast(result, num_clusters ~ genes,value.var = 'p.value', fun.aggregate=mean)
methods_all <- data2$num_clusters
pval = data.frame(t(data2[, -1]))
colnames(pval) <- methods_all
methods_all <- data$num_clusters
padj = data.frame(t(data[, -1]))
colnames(padj) <- methods_all


DF_COBRA <- COBRAData(pval,
                      padj,
                      truth = data.frame(df2))
perf <- calculate_performance(DF_COBRA, binary_truth = "status",
                              thrs = c(0.01, 0.05, 0.1, 0.2))
cbbPalette <- c("#000000", "#FFD92F", "#56B4E9", "#33A02C", "#F4A582", "#C51B7D")

cobra_plot <- prepare_data_for_plot(perf, colorscheme = cbbPalette, incloverall = FALSE,
                                    facetted = TRUE,conditionalfill = FALSE,
                                    keepmethods = c("2","4","6","8","10","12")
)

shape_border = c(0, 1, 2, 5, 6, 3)
shape_fill = c(15, 16, 17, 23, 25, 3)
## Manually match "method" and point "shape"
sel_shape = c(5,6,1,2,3,4)
sel_fill = sel_shape
col <- c(cbbPalette)
cobra_plot@plotcolors <- cobra_plot@plotcolors[1:6]
(gg_BayesSpace <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                                   pointsize = 0, linewidth = 2)+
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.8)) +
    scale_color_manual(values = rep(cbbPalette[c(5,6,1,2,3,4)],4),
                       name = "",
                       breaks=c('2','4','6','8','10','12'))+
    guides(colour = guide_legend(ncol = 6, byrow = FALSE,
                                 override.aes = list(shape = shape_fill,
                                                     fill = cbbPalette[c(5,6,1,2,3,4)]) ) ) +
    geom_point(size = 10, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_border[sel_shape],4), 
    stroke = 2, 
    alpha = 1) + 
    scale_fill_manual(values =rep(cbbPalette,4), 
                      name = "",
                      breaks = c('2','4','6','8','10','12'))+
    geom_point(size = 10, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_fill[sel_fill],4), 
    stroke = 2, alpha = 0.25)+
    labs(title = "BayesSpace") + 
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          title = element_text(size = rel(3)),
          axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(3)),
          axis.text.y=element_text(size=rel(3)),
          axis.title.y = element_text(size=rel(1)),
          axis.title.x = element_text(size=rel(1)),
          legend.title=element_text(size=rel(1)),
          legend.text=element_text(size=rel(3)),
          legend.key.width=unit(2.5, "cm"),
          aspect.ratio = 1, legend.position="bottom",
          legend.box="vertical", legend.margin=margin()) )

################################ DESpace based on real clusters ########################

ll <- lapply(c(2,4,6,8,10,12), function(x){
  print(x)
  path <- "~/Desktop/master_thesis/Data_re_run/simulation/results/151507/split_mixture_patch_"
  results_dir <- paste0(path,x,"/probs_0.5_0.9FALSE/")
  sce_dir <- paste0(results_dir, "sce_edgeR.rda")
  SVG_dir <- paste0(results_dir, "SV_genes.txt")
  DESpace_dir <- paste0(results_dir, "result_SV_edgeR_original_clusters.rda")
  load(DESpace_dir)
  result <- spatial_gene_results
  genes <- read.csv(SVG_dir,header = TRUE,sep =  '\t')
  if(!all(genes$V1 %in% SVG_list$V1)){
    message("Please check the list of spatial variable genes; 
              it should be matched to the ground truth.")
  }
  l2 <-cbind(result, num_clusters = x)
  l2 <- l2[,c("genes","PValue", "FDR", "method", "num_clusters")]
  return(l2)
}
)
result <- do.call(rbind, ll)

all_genes <- as.data.frame((unique(result$genes)))

all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
genes_id <- apply(SVG_list,1,function(x) substr(x, 10,18))
status <- ifelse(all_genes_id %in% genes_id, 1, 0)
truth <- cbind(all_genes, status)
# check
sum(truth$status) == dim(SVG_list)[1]
if(sum(all_genes_id[length(all_genes_id)] == "_count")>=1){
  truth <- truth[1:dim(truth)[1]-1,]
}
result <- result %>% dplyr::distinct()
setDT(result)
result$PValue <- as.numeric(result$PValue)
result$FDR <- as.numeric(result$FDR)
data <- dcast(result, num_clusters ~ genes,value.var = 'FDR', fun.aggregate=mean)
#data <- dcast(result, method ~ genes,value.var = 'adj.p.value')
df2 = truth[match(colnames(data)[-1], truth$`(unique(result$genes))`),]
df2 <- as.data.frame(df2[,2])
rownames(df2) <- colnames(data)[-1]
colnames(df2) <- 'status'
data2 <- dcast(result, num_clusters ~ genes,value.var = 'PValue', fun.aggregate=mean)
methods_all <- data2$num_clusters
pval = data.frame(t(data2[, -1]))
colnames(pval) <- methods_all
methods_all <- data$num_clusters
padj = data.frame(t(data[, -1]))
colnames(padj) <- methods_all


DF_COBRA <- COBRAData(pval,
                      padj,
                      truth = data.frame(df2))
perf <- calculate_performance(DF_COBRA, binary_truth = "status",
                              thrs = c(0.01, 0.05, 0.1, 0.2))
cbbPalette <- c("#000000", "#FFD92F", "#56B4E9", "#33A02C", "#F4A582", "#C51B7D")

cobra_plot <- prepare_data_for_plot(perf, colorscheme = cbbPalette, incloverall = FALSE,
                                    facetted = TRUE,conditionalfill = FALSE,
                                    keepmethods = c("2","4","6","8","10","12")
)

shape_border = c(0, 1, 2, 5, 6, 3)
shape_fill = c(15, 16, 17, 23, 25, 3)
## Manually match "method" and point "shape"
sel_shape = c(5,6,1,2,3,4)
sel_fill = sel_shape
col <- c(cbbPalette)
cobra_plot@plotcolors <- cobra_plot@plotcolors[1:6]
(gg_real <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                             pointsize = 0, linewidth = 2)+
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.8)) +
    scale_color_manual(values = rep(cbbPalette[c(5,6,1,2,3,4)],4),
                       name = "",
                       breaks=c('2','4','6','8','10','12'))+
    guides(colour = guide_legend(ncol = 6, byrow = FALSE,
                                 override.aes = list(shape = shape_fill,
                                                     fill = cbbPalette[c(5,6,1,2,3,4)]) ) ) +
    geom_point(size = 10, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_border[sel_shape],4), 
    stroke = 2, 
    alpha = 1) + # stroke = line width
    scale_fill_manual(values =rep(cbbPalette,4), 
                      name = "",
                      breaks = c('2','4','6','8','10','12'))+
    geom_point(size = 10, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_fill[sel_fill],4), 
    stroke = 2, alpha = 0.25)+
    labs(title = "Reference") + 
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          title = element_text(size = rel(3)),
          axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(3)),
          axis.text.y=element_text(size=rel(3)),
          axis.title.y = element_text(size=rel(1)),
          axis.title.x = element_text(size=rel(1)),
          legend.title=element_text(size=rel(1)),
          legend.text=element_text(size=rel(3)),
          legend.key.width=unit(2.5, "cm"),
          aspect.ratio = 1, legend.position="bottom",
          legend.box="vertical", legend.margin=margin()) )


################################ Aggregate figures ########################

gg <- ggarrange(gg_real,gg_BayesSpace, gg_stLearn, 
                common.legend = TRUE, 
                ncol = 3, legend="bottom")
path_save = "./Figures/Figures/Supplementary/"
ggsave(filename = paste0(path_save,"SuppFig_FDR_varing_clusters.pdf"),
       plot = gg,
       device = "pdf",
       width = 30,
       height = 10,
       units = "in",
       dpi = 300,
       limitsize = TRUE)




