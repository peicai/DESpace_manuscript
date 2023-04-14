rm(list = ls())
path = "./Real/"
path_save = "./Figures/Figures/Real/"
source("./Figures/Scripts/Plots_Function.R")
############################ For LIBD ###################################
data = 'LIBD'
sample_names = c("151507", "151508", "151509", "151510",
                 "151669", "151670", "151671", "151672",
                 "151673", "151674", "151675", "151676")
# list_result <- process_data(
#   path = path, Manual = TRUE,
#   data = 'LIBD', threshold = 0.01,
#   name_head = -10, name_tail = -5
# )

list_result_LIBD <- list()
list_result_LIBD <- process_data(
  path = path,
  data = 'LIBD', Manual = TRUE, 
  threshold = 2, # indicates that we do not set threshold
  name_head = -10,name_tail = -5
)

gg_Bar_LIBD <- ave_jaccard_btw_rep_LIBD(
  list_test_all=list_result_LIBD,
  path = path_save,
  sample_names=sample_names,
  top_genes = 1000,
  colors_method = colors_method,
  data = 'LIBD')

####################### For computational cost ##########################
library(lubridate)
setwd(paste0(path, data,"/"))
computational_cost <- read_csv("computational_cost.csv")
df = computational_cost
(LIBD_files <- list.files(pattern = "*SpaGCN_runtime.csv$", full.names = TRUE))
SpaGCN <- do.call(rbind,lapply(LIBD_files,read.csv))
SpaGCN$times <- period_to_seconds(hms(gsub("0 days ", "", SpaGCN$X0)))
(LIBD_files <- list.files(pattern = "*spatialDE_runtime.csv$", full.names = TRUE))
spatialDE <- do.call(rbind,lapply(LIBD_files,read.csv))
spatialDE$times <- period_to_seconds(hms(gsub("0 days ", "", spatialDE$X0)))
(LIBD_files <- list.files(pattern = "*spatialDE2_runtime.csv$", full.names = TRUE))
spatialDE2 <- do.call(rbind,lapply(LIBD_files,read.csv))
spatialDE2$times <- period_to_seconds(hms(gsub("0 days ", "", spatialDE2$X0)))
(LIBD_files <- list.files(pattern = "*stLearn_runtime.csv$", full.names = TRUE))
stLearn <- do.call(rbind,lapply(LIBD_files,read.csv))
stLearn$times <- period_to_seconds(hms(gsub("0 days ", "", stLearn$X0)))
computational_cost_spark_ <- read_csv("computational_cost_spark_.csv")
spark <- c(rep(computational_cost_spark_$V1,12))
computational_cost_stLearn_manual_SV <- read_csv("computational_cost_stLearn_manual_SV.csv")
computational_cost_nnSVG <- read_csv("computational_cost_nnSVG.csv")
computational_cost_nnSVG2 <- read_csv("computational_cost_nnSVG2.csv")
df <- cbind(df,SpaGCN$times,spatialDE$times,spatialDE2$times,stLearn$times+computational_cost_stLearn_manual_SV$V1,
            computational_cost_nnSVG$V1,spark)
colnames(df) <- c("times" ,"Manual_DESpace","BayesSpace_DESpace","SPARK-X","MERINGUE","sample",
                  "SpaGCN","SpatialDE","SpatialDE2" ,                               
                  "StLearn_DESpace","nnSVG","SPARK")
df = as.data.frame(df)
df2 <- t(as.data.frame(colMeans(df[,c(2:5,7:12)])))
rownames(df2) <- c("Average")
df2 <- reshape::melt(df2)
df2$X1 <- as.factor(df2$X1)
df2$value <- df2$value/60
gg_data <- df2
colnames(gg_data) <- c("Sample_id", "Method", "Time_Minutes")

gg_data$Method <- factor(gg_data$Method,levels=c("Manual_DESpace","BayesSpace_DESpace","StLearn_DESpace",
                                                 "MERINGUE", "nnSVG",
                                                 "SPARK","SPARK-X",
                                                 "SpatialDE", "SpatialDE2",
                                                 "SpaGCN"))
Time <- round(gg_data$Time_Minutes)
gg_data$Method <- factor(gg_data$Method, levels=c(methods_all), labels=c(methods_all))
## Common legend
gg_legend <- ggplot(gg_data,
       #aes_string(y = "Time_Minutes", fill = "Method"),
       aes(x = reorder(Method,Time_Minutes), y = Time_Minutes, fill = Method)) +  
  geom_bar( stat = "identity",position = position_dodge(0.8), width = 0.7) +  theme_bw() +   xlab("") +  ylab("Time (Minutes)") + 
  scale_y_sqrt( breaks = c(5,10,20,50,100,500,1000,2000) ) + 
  geom_text(aes(label=Time), vjust=-0.5,size = 3.5, fontface = "bold") +
  scale_fill_manual("Method", values =colors_method)+
  theme(axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1),   
        axis.text.y= element_text(size=rel(0.8),angle = 0, hjust = 1),  
        axis.title = element_text(size=rel(1)),      
        panel.grid.minor = element_blank(),         
        panel.grid.major.x = element_blank(),       
        panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),        
        aspect.ratio = 1,        
        legend.position = "bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=30),
        legend.key.width=unit(2, "cm"),
        legend.box="vertical",
        legend.key.height=unit(2.25,"line"),
        legend.margin=margin()) +
  guides(colour = guide_legend(ncol = 10, byrow = FALSE,
                               override.aes = list(fill = colors_method) ) ) 
## Main plots
colors = as.data.frame(colors_method)
colors$Method = rownames(colors)
dd_gg_data <- gg_data %>% left_join(colors)

dd_gg_data$Method_sub <- dd_gg_data$Method
## split BayesSpace_DESpace and StLearn_DESpace into 2 parts, respectively: clustering tool and DESpace parts.
dd_gg_data <- rbind(dd_gg_data, dd_gg_data[c(2,8),])
dd_gg_data[c(11,12),]$Method_sub <- "Clustering tool"
dd_gg_data[c(11,12),]$Time_Minutes <- c(dd_gg_data[2,]$Time_Minutes - dd_gg_data[1,]$Time_Minutes,
                                        mean(stLearn$times/60))
dd_gg_data$Time_All <- dd_gg_data$Time_Minutes
dd_gg_data[c(2,8),]$Time_Minutes <- c(dd_gg_data[1,]$Time_Minutes,
                                        mean(computational_cost_stLearn_manual_SV$V1/60))
color_all <- c(colors_method, colours_CompuCost[1])
colors = as.data.frame(color_all)
colors$Method_sub = rownames(colors)
dd_data <- rbind(dd_gg_data) %>% left_join(colors)
Time <- round(dd_data$Time_Minutes)
(gg_Time_LIBD <- dd_data %>% 
  ggplot(mapping = aes(x = reorder(Method,Time_All), Time_Minutes, fill = color_all)) + 
  geom_bar( stat = "identity",position = "stack", width = 0.7) +  
  theme_bw() +   xlab("") +  ylab("Time (Minutes)") + 
    geom_text(aes(y = Time_Minutes, label=Time), 
              vjust=1.6,size = 3.5, fontface = "bold") +
    # geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
    #           color="white", size=3.5)+
    # vjust=-0.5,size = 3.5,
  scale_y_sqrt( breaks = c(5,10,20,50,100,500,1000,2000,5000) ) + 
  scale_fill_identity() +
    theme(axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1),   
          axis.text.y= element_text(size=rel(0.8),angle = 0, hjust = 1),  
          axis.title = element_text(size=rel(1)),      
          panel.grid.minor = element_blank(),         
          panel.grid.major.x = element_blank(),       
          panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),        
          aspect.ratio = 1,        
          legend.position = "none",
          legend.text=element_text(size=15)
    ) )

combined_legend<- dd_data[c(2,8,11,12),c("Method", "Method_sub", "color_all")] %>% 
  ggplot(mapping = aes(y = Method, x = factor(Method_sub), fill = color_all)) + 
  geom_tile() + 
  scale_fill_identity() + 
  coord_fixed() +
  theme_bw()


###################################### For melanoma ################
sample_names = c("mel1_rep1", "mel1_rep2",
                 "mel2_rep1", "mel2_rep2",
                 "mel3_rep1", "mel3_rep2",
                 "mel4_rep1", "mel4_rep2")
list_result_melanoma <- list()
list_result_melanoma <- process_data(
  path = path,
  data = 'melanoma', Manual = FALSE, 
  threshold = 2, # indicates that we do not set threshold
  name_head = -13,name_tail = -5
)

gg_Bar_melanoma <- ave_jaccard_btw_rep(
  list_test_all=list_result_melanoma,
  path = path_save,
  sample_names=sample_names,
  top_genes = 200,
  colors_method = colors_method,
  data = 'melanoma')

gg_Bar_LIBD2 <- gg_Bar_LIBD +
  theme(axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1),   
        axis.text.y= element_text(size=30,angle = 0, hjust = 1), 
  ) 
gg_Bar_melanoma2 <- gg_Bar_melanoma +
  theme(axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1),   
        axis.text.y= element_text(size=30,angle = 0, hjust = 1), 
  ) 
(AA = egg::ggarrange( plots = list(gg_Bar_LIBD2 +ylab("")+ labs(title = "LIBD") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size=35)),
                            gg_Bar_melanoma2 +  labs(title = "melanoma") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size=35))
                            ),
                     bottom = ggpubr::get_legend(gg_legend),
                     ncol = 2, nrow = 1))
ggsave(filename = paste0(path_save,"Barplot_Jaccrd.pdf"),
       plot = AA,
       device = "pdf",
       width = 25,
       height = 15,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
