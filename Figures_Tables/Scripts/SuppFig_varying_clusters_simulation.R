rm(list = ls())
library(spdep)
## load data:
LIBD_files <- list.files(path = "./DESpace_data/Data/LIBD",
                         pattern = "^15 *.*rda$", full.names = TRUE)

sample_names <- substr(gsub(".rda.*$", "", gsub("^.*/", "", LIBD_files)),1,6)
List_sce_one = list()

source("./Analyses/03_simulations/02_simulation_function_data_specific.R")

## load filtered melanoma data
set.seed(169612)
for(ii in 1:length(sample_names)){
  load(LIBD_files[ii])
  sample <- sample_names[ii]
  layer = colData(sce_one)$layer_guess_reordered
  List_sce_one[[ii]] = sce_one
} 
## pick the first sample
i = 1
set.seed(123)
sce_one = List_sce_one[[i]]
sce_one = sce_one[, !is.na(sce_one$layer_guess_reordered)]
exprs <- assays(sce_one)$counts
sce_one <- addPerCellQCMetrics(sce_one)
sce_one <- logNormCounts(sce_one)
CD <- colData(sce_one)
cell_ID = rownames(colData(sce_one))
colData(sce_one) <- cbind(CD, cell_ID)
sample_genes = rownames(exprs)

object_subset = sce_one[sample_genes,]
coordinate_name <- c("imagerow", "imagecol")
covariates = cbind.data.frame(x=as.numeric(as.character(colData(object_subset)[[coordinate_name[1]]])),
                              y=as.numeric(as.character(colData(object_subset)[[coordinate_name[2]]])),
                              cell_ID = rownames(colData(object_subset))
)
colnames(covariates) <- c("sdimx","sdimy","cell_ID")
CD <- as.data.frame(CD)

## split the original 7 clusters (manual annotations) into 2, 4, 6, 8, 10 and 12 clusters

CD$pattern_ids <- with(CD, sign(imagerow + 1.1*imagecol - 635))
CD$pattern_ids2 <- with(CD, sign(imagerow + 1.1*imagecol - 480))
CD$pattern_ids3 <- with(CD, sign(imagerow + 1.1*imagecol - 700))
creat_clus <- function(n_clus){
  CD$layer_23 <- with(CD, sign(imagerow - 0.7*imagecol + 120))
  CD$layer_new <- as.character(CD$layer_guess_reordered)
  CD$layer_new <- ifelse(CD$layer_guess_reordered %in% c("Layer2", "Layer3") &
                           CD$layer_23 == 1,
                         paste0(CD$layer_guess_reordered, "_2nd"),
                         as.character(CD$layer_guess_reordered))
  if(n_clus == 12){
    CD$layer_new[CD$layer_new %in% c("Layer3_2nd", "Layer2_2nd")] <- "Layer32_2nd"
    CD$layer_new[CD$layer_new %in% c("Layer3", "Layer2")] <- "Layer32"
    CD$layer_new[CD$layer_new %in% c("Layer5", "Layer4")] <- "Layer45"
    CD$layer_new[CD$layer_new %in% c("WM", "Layer6")] <- "WM6"
    CD$pattern <- CD$layer_new
    CD$pattern <- ifelse(CD$pattern %in% c("Layer32_2nd", "Layer45",
                                           "Layer1") &
                           CD$pattern_ids3 == 1,
                         as.character(paste0(CD$pattern, "1")),
                         CD$pattern)
    CD$pattern <- ifelse(CD$pattern %in% c("Layer32_2nd", "Layer45",
                                           "Layer1") &
                           CD$pattern_ids3 == -1&
                           CD$pattern_ids2 == 1,
                         as.character(paste0(CD$pattern, "0")),
                         CD$pattern)
    CD$pattern <- ifelse(CD$pattern %in% c("Layer32_2nd", "Layer45",
                                           "Layer1") &
                           CD$pattern_ids2 == 1,
                         as.character(paste0(CD$pattern, "-1")),
                         CD$pattern)
    CD$pattern <- ifelse(CD$pattern %in% c("WM6") &
                           CD$pattern_ids3 == 1,
                         "WM6-1",
                         CD$pattern)
  }else if(n_clus == 10){
    CD$layer_new[CD$layer_new %in% c("Layer3_2nd", "Layer2_2nd")] <- "Layer32_2nd"
    CD$layer_new[CD$layer_new %in% c("Layer3", "Layer2")] <- "Layer32"
    CD$layer_new[CD$layer_new %in% c("Layer5", "Layer4")] <- "Layer45"
    CD$layer_new[CD$layer_new %in% c("WM", "Layer6")] <- "WM6"
    CD$pattern <- paste0(CD$layer_new, CD$pattern_ids)
    CD$pattern[CD$pattern %in% c("WM61", "WM6-1")] <- "WM6"
    CD$pattern[CD$pattern %in% c("Layer32-1", "Layer321")] <- "Layer32"
    CD$pattern[CD$pattern %in% c("Layer32_2nd-1", "Layer32_2nd1")] <- "Layer32_2nd"
    CD$pattern <- ifelse(CD$pattern %in% c("Layer32_2nd") &
                           CD$pattern_ids3 == 1,
                         "Layer32_2nd1",
                         CD$pattern)
    CD$pattern <- ifelse(CD$pattern %in% c("Layer32_2nd") &
                           CD$pattern_ids3 == -1&
                           CD$pattern_ids2 == 1,
                         "Layer32_2nd-1",
                         CD$pattern)
    CD$pattern <- ifelse(CD$pattern %in% c("Layer32_2nd") &
                           CD$pattern_ids2 == 1,
                         "Layer32_2nd1",
                         CD$pattern)
    CD$pattern <- ifelse(CD$pattern %in% c("WM6") &
                           CD$pattern_ids3 == 1,
                         "WM6-1",
                         CD$pattern)
  }else if(n_clus == 8){
    CD$layer_new[CD$layer_new %in% c("Layer3_2nd", "Layer2_2nd")] <- "Layer32_2nd"
    CD$layer_new[CD$layer_new %in% c("Layer3", "Layer2")] <- "Layer32"
    CD$layer_new[CD$layer_new %in% c("Layer5", "Layer4")] <- "Layer45"
    CD$layer_new[CD$layer_new %in% c("WM", "Layer6")] <- "WM6"
    CD$pattern <- paste0(CD$layer_new, CD$pattern_ids)
    CD$pattern[CD$pattern %in% c("WM6-1", "WM61")] <- "WM6"
    CD$pattern[CD$pattern %in% c("Layer32-1", "Layer321")] <- "Layer32"
  }else if(n_clus == 6){
    CD$pattern <- paste0(CD$layer_new, CD$pattern_ids)
    CD$pattern[CD$pattern %in% c("WM1","Layer61", "Layer51","Layer41")] <- "A"
    CD$pattern[CD$pattern %in% c("WM-1","Layer6-1", "Layer5-1","Layer4-1")] <- "B"
    CD$pattern[CD$pattern %in% c("Layer3_2nd-1", "Layer2_2nd-1")] <- "C"
    CD$pattern[CD$pattern %in% c("Layer3_2nd1", "Layer2_2nd1")] <- "D"
    CD$pattern[CD$pattern %in% c("Layer31", "Layer21","Layer11")] <- "E"
    CD$pattern[CD$pattern %in% c("Layer3-1", "Layer2-1","Layer1-1")] <- "F"
  }else if(n_clus == 4){
    CD$pattern <- paste0(CD$layer_guess_reordered, CD$pattern_ids)
    CD$pattern[CD$pattern %in% c("WM-1", "Layer6-1", "Layer5-1","Layer4-1")] <- "A"
    CD$pattern[CD$pattern %in% c("WM1", "Layer61", "Layer51","Layer41")] <- "B"
    CD$pattern[CD$pattern %in% c("Layer1-1", "Layer2-1","Layer3-1")] <- "C"
    CD$pattern[CD$pattern %in% c("Layer11", "Layer21","Layer31")] <- "D"
  }else if(n_clus == 2){
    CD$pattern <- as.character(CD$pattern_ids)
  }
  
  ## figure
  
  mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))
  (p <- ggplot(CD, aes(x = imagecol, imagerow)) +
      scale_x_reverse() +
      geom_point(aes(colour = factor(pattern)),size = 0.5)+
      
      theme_classic() +
      scale_color_manual(values = mycolors) +
      labs(title = paste0("Number of Clusters: ", n_clus)) +
      theme(legend.position = "none",
            line = element_blank(),
            text=element_text(size=10),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
      ))
  return(p)
}

p_all <- lapply(c(2,4,6,8,10,12), creat_clus)
pp <- gridExtra::grid.arrange(grobs = p_all,
                              ncol = 3)

## save the figure
path_save = "./Figures/Figures/Simulated/"
ggsave(paste0(path_save,"SuppFig6_varing_clusters_example.pdf"),
       plot = pp,
       width = 6,
       height = 5)
