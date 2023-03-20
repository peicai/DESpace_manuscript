library(RColorBrewer);
all_colours = c(
  "#08306B", # "Manual_DESpace"
  "#A6CEE3", # "BayesSpace_DESpace"
  "#1F78B4", # "StLearn_DESpace"
  "#FFD92F", # "SPARK"
  "#A6761D", # "SPARK-X"
  "#B2DF8A", # "SpatialDE"
  "#33A02C", # "SpatialDE2"
  "#F4A582", # "MERINGUE"
  "#666666", # "SpaGCN"
  "#C51B7D"  # "nnSVG"
)
colours_CompuCost = c(
   "DESpace" = "#08306B",#"#41B6C4",
   "BayesSpace" = "#A6CEE3",
   "StLearn" = "#1F78B4"
)
#"#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6" "#4292C6" "#2171B5" "#08519C" "#08306B"
#"#FFF7FB" "#ECE7F2" "#D0D1E6" "#A6BDDB" "#74A9CF" "#3690C0" "#0570B0" "#045A8D" "#023858"
#"#FFFFD9" "#EDF8B1" "#C7E9B4" "#7FCDBB" "#41B6C4" "#1D91C0" "#225EA8" "#253494" "#081D58"
#"#FFF7FB" "#ECE7F2" "#D0D1E6" "#A6BDDB" "#74A9CF" "#3690C0" "#0570B0" "#045A8D" "#023858"
#"#A63603" "#E6550D" "#FD8D3C" "#FDBE85"
#"#006837" "#31A354" "#78C679" "#C2E699"
#"#E66101" "#FDB863" "#F7F7F7" "#B2ABD2" "#5E3C99"
#"#CA0020" "#F4A582" "#FFFFFF" "#BABABA" "#404040"
#brewer.pal(5, "RdGy")[1:5]    # scDD-perm 4 methods
methods_all = c(
  "Manual_DESpace",
  "BayesSpace_DESpace",
  "StLearn_DESpace",
  "SPARK",
  "SPARK-X",
  "SpatialDE",
  "SpatialDE2",
  "MERINGUE",
  "SpaGCN",
  "nnSVG"
)

colors_method = all_colours
names(colors_method) = methods_all
# methods names:
methods_order = c(
  "SV_edgeR_counts",
  "StLearn_edgeR",
  "spark",
  "spark_x",
  "spatialDE",
  "spatialDE2",
  "meringue",
  "SpaGCN",
  "nnSVG"
)

# points, borders:
shape_border = c(0, 1, 2, 5, 6, 3, 4, 8, 11)
# points, fill:
shape_fill = c(15, 16, 17, 23, 25, 3, 4, 8, 11)


# ordering = order(methods_order)
# 
# all_methods = all_methods[ordering]
# methods_names = methods_names[ordering]
# all_colours = all_colours[ordering]

#' path_data = paste0(path, dataset, patterns[i], "/probs_", spatial_probs[1], "_", spatial_probs[2], "FALSE/")
#' load(paste0(path_data, "probs_", spatial_probs[1], "_", spatial_probs[2], "_final_object",sample_names[1], '.rda'))
#' if(!exists("combined_object")){combined_object = final_object; remove(final_object)}
#' sce1 = combined_object; remove(combined_object)
#' load(paste0(path_data, "probs_", spatial_probs[1], "_", spatial_probs[2], "_final_object",sample_names[2], '.rda'))
#' if(!exists("combined_object")){combined_object = final_object; remove(final_object)}
#' sce2 = combined_object; remove(combined_object)
#' load(paste0(path_data, "probs_", spatial_probs[1], "_", spatial_probs[2], "_final_object",sample_names[3], '.rda'))
#' if(!exists("combined_object")){combined_object = final_object; remove(final_object)}
#' sce3 = combined_object; remove(combined_object)
#' gene1 <- rownames(sce1);gene2 <- rownames(sce2);gene3 <- rownames(sce3)
#' sum(gene1 == gene2)== 5000; sum(gene1 == gene3)== 5000
#' 
#' load(paste0(path_data,"result_SV_edgeR_counts.rda"))
#' result = result1 = spatial_gene_results
#' file = paste0(path_data, 'probs_0.8_selected_genes.txt')
#' genes <- read.csv(file,header = TRUE,sep =  '\t')
#' all_genes <- as.data.frame((unique(result$genes)))
#' 
#' all_genes_id <- apply(all_genes,1,function(x) substr(x, 10,18))
#' genes_id <- apply(genes,1,function(x) substr(x, 10,18))
#' status <- ifelse(all_genes_id %in% genes_id, 1, 0)
#' truth <- cbind(all_genes, status)
#' # check
#' sum(truth$status) == dim(genes)[1]
#' 
#' result <- result %>% dplyr::distinct()
#' setDT(result)
#' ################################
#' methods_order = c('Original_single_DESpace','Original_multi_DESpace'#,
#'                   #'stLearn_single_DESpace','stLearn_multi_DESpace'
#' )
#' single_Bayes = as.data.table(result1 %>% filter(method == "SV_single_edgeR" & cluster == "Original"))
#' multi_Bayes = as.data.table(result1 %>% filter(method == "SV_multi_edgeR" & cluster == "Original"))
#' single_Bayes[, methods := 'Original_single_DESpace']
#' multi_Bayes[, methods := 'Original_multi_DESpace']
#' 
#' single_Bayes = single_Bayes[,.(genes, FDR, PValue, sample, methods)]
#' multi_Bayes = multi_Bayes[,.(genes, FDR, PValue, sample, methods)]
#' 
#' samples = unique(single_Bayes$sample)
#' combined_data <- data.frame(0)
#' combined_data2 <- data.frame(0)
#' combined_df2 <- NULL
#' combined_result <- NULL
#' combined_genes <- NULL
#' combined_all_genes <- NULL
#' 
#' for(ii in c(1:length(samples))){
#'   combine_results = rbind(
#'     
#'     single_Bayes %>% filter(sample == samples[ii]),
#'     multi_Bayes
#'   )
#'   setDT(combine_results)
#'   data <- dcast(combine_results, methods ~ genes,value.var = 'FDR')
#'   df2 = truth[match(colnames(data)[-1], truth$`(unique(result$genes))`),]
#'   df2 <- as.data.frame(df2[,2])
#'   rownames(df2) <- colnames(data)[-1]
#'   colnames(df2) <- 'status'
#'   data2 <- dcast(combine_results, methods ~ genes,value.var = 'PValue')
#'   
#'   colnames(data) <- paste0(colnames(data), ".",ii)
#'   colnames(data2) <- paste0(colnames(data2),".", ii)
#'   rownames(df2) <- paste0(rownames(df2), ".",ii)
#'   
#'   colnames(data)[1] <- "method"
#'   colnames(data2)[1] <- "method"
#'   
#'   #data <- data[,1:(dim(data)[2]-1)]
#'   #data2 <- data2[,1:(dim(data2)[2]-1)]
#'   combined_data <- cbind(combined_data,data)
#'   combined_data2 <- cbind(combined_data2,data2)
#'   combined_df2 <- rbind(combined_df2,df2)
#' }
#' 
#' combined_data <- combined_data[,-1]
#' combined_data2 <- combined_data2[,-1]
#' 
#' combined_data <- t(combined_data)
#' methods <- combined_data[1,]
#' combined_data <- combined_data[-1,]
#' colnames(combined_data) <- methods
#' 
#' combined_data2 <- t(combined_data2)
#' combined_data2 <- combined_data2[-1,]
#' colnames(combined_data2) <- methods
#' 
#' genes_id <- rownames(combined_data2)
#' combined_data2 <- apply(combined_data2, 2,            # Specify own function within apply
#'                         function(x) as.numeric(as.character(x)))
#' rownames(combined_data2) <- genes_id
#' 
#' combined_data <- apply(combined_data, 2,            # Specify own function within apply
#'                        function(x) as.numeric(as.character(x)))
#' rownames(combined_data) <- genes_id
#' combined_data <- combined_data[,match(methods_order, colnames(combined_data))]
#' combined_data2 <- combined_data2[,match(methods_order, colnames(combined_data2))]
#' 
#' DF_COBRA <- COBRAData(pval = data.frame(combined_data2
#' ),
#' padj = data.frame(combined_data
#' ),
#' truth = data.frame(combined_df2))
#' # truth = 1 for SV genes and 0 for uniform genes.
#' 
#' #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#' # ROC, FDR plots:
#' #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#' perf <- calculate_performance(DF_COBRA, binary_truth = "status",
#'                               thrs = c(0.01, 0.05, 0.1, 0.2))
#' cobra_plot <- prepare_data_for_plot(perf, colorscheme = colours,incloverall = FALSE,
#'                                     facetted = TRUE,conditionalfill = FALSE)
#' cobra_plot@plotcolors <- cobra_plot@plotcolors[1:2]
#' 
#' (gg_roc <- plot_roc(cobra_plot,linewidth=2)  + 
#'     #scale_fill_brewer(breaks=methods, type="qual", palette=3)+
#'     # scale_color_manual(values = cobra_plot@plotcolors, name = "", labels = methods_order,limits = force) +
#'     theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
#'                                      hjust = 1, size = 10),
#'           axis.text.y = element_text(size = 10),
#'           axis.title.x = element_text(size = 10),
#'           axis.title.y = element_text(size = 10),
#'           legend.position = "top") )
#' (gg_roc <- gg_roc +
#'     guides(color=guide_legend(nrow=1,byrow=TRUE)))
#' 
#' #ggsave(paste0(path,roc_save_name),width = 7, height = 5)
#' col = colours[-length(colours)]
#' names(col) = methods_order
#' # plot FDR/TPR curve
#' gg_fdr <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
#'                            pointsize = 0, linewidth = 2)+ 
#'   scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1), limits = c(0,0.2)) +
#'   scale_color_manual(values = col,
#'                      name = "",
#'                      breaks=methods_all,
#'                      labels=methods_all)+
#'   theme(strip.background = element_blank(),
#'         strip.text = element_blank(),
#'         axis.text.x = element_text(angle = 45, hjust=0.5, size = rel(3)),
#'         axis.text.y=element_text(size=rel(3)),
#'         axis.title.y = element_text(size=rel(3)),
#'         axis.title.x = element_text(size=rel(3)),
#'         legend.title=element_text(size=rel(2)),
#'         legend.text=element_text(size=rel(3)),
#'         legend.key.width=unit(2, "cm"),
#'         aspect.ratio = 1, legend.position="bottom",
#'         legend.box="vertical", legend.margin=margin())  +
#'   guides(colour = guide_legend(ncol = 2, byrow = FALSE,
#'                                override.aes = list(shape = shape_fill,
#'                                                    fill = colours[-length(colours)]) ) ) +
#'   geom_point(size = 6, aes(fill = method, colour =method, shape = method), 
#'              shape = rep(shape_border,4), stroke = 2, alpha = 1) + # stroke = line width
#'   scale_fill_manual(values =col, guide = "none") +
#'   geom_point(size = 6, aes(fill = method, colour = method, shape = method), 
#'              shape = rep(shape_fill,4), stroke = 2, alpha = 0.25)
#' 
