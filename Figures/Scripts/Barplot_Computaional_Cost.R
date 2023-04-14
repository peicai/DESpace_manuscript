rm(list = ls())
library(ggplot2)
library(dplyr)
library(forcats)
library(readr)
library(lubridate)
library(ggpattern)
path = "./Real/"
path_save = "./Figures/Figures/Real/"
source("./Figures/Scripts/Plots_Function.R")
source("./Figures/Scripts/All_methods.R")
color_bar = c("#A63603","#A6CEE3","#1F78B4","#08306B","#FFD92F",
                       "#A6761D","#F4A582","#B2DF8A",
                       "#33A02C","#666666","#C51B7D")
                       ############################ For LIBD ###################################
data = 'LIBD'
sample_names = c("151507", "151508", "151509", "151510",
                 "151669", "151670", "151671", "151672",
                 "151673", "151674", "151675", "151676")
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
computational_cost_BayesSpace_DESpace_seperate <- read_csv("computational_cost_BayesSpace_DESpace_seperate.csv")
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
Time <- round(gg_data$Time_Minutes)
gg_data$Method <- factor(gg_data$Method, levels=c(methods_all), labels=c(methods_all))

## Main plots
colors = as.data.frame(colors_method)
colors$Method = rownames(colors)
dd_gg_data <- gg_data %>% left_join(colors)
Time <- round(dd_gg_data$Time_Minutes)
label_ypos = dd_gg_data$Time_Minutes
Methods <- arrange(dd_gg_data, Time_Minutes)$Method
dd_gg_data$Method_sub <- dd_gg_data$Method
## split BayesSpace_DESpace and StLearn_DESpace into 2 parts, respectively: clustering tool and DESpace parts.
dd_gg_data <- rbind(dd_gg_data, dd_gg_data[c(2,8),])
dd_gg_data[c(11,12),]$Time_Minutes <- c(mean(computational_cost_BayesSpace_DESpace_seperate$BayesSpace/60),
                                             mean(stLearn$times/60))
#dd_gg_data[c(2,8),]$Method_sub <- "DESpace"
dd_gg_data[c(2,8),]$Time_Minutes <- c(mean(computational_cost_BayesSpace_DESpace_seperate$edgeR/60),
                                     mean(computational_cost_stLearn_manual_SV$V1/60))
dd_gg_data$Time_All <- dd_gg_data$Time_Minutes
dd_gg_data[c(11,12),]$Method_sub <- c("BayesSpace", "StLearn")
color_all <- c(colors_method, colours_CompuCost)
colors = as.data.frame(color_all)
colors$Method_sub = rownames(colors)
dd_data <- rbind(dd_gg_data) %>% left_join(colors)
d3 <- dd_data %>% 
mutate(color_all = factor(color_all, levels = color_bar))

`%notin%` <- Negate(`%in%`)
d3$Cluster <- as.factor(ifelse(d3$Method_sub %in% c('BayesSpace','StLearn'), 'Cluster', 'Others'))  
d3 %>% mutate(highligth = case_when(d3$Method_sub %in% c('BayesSpace','StLearn') ~ "yes", 
                                    TRUE ~ "no"))                      
#remotes::install_github("coolbutuseless/ggpattern")
d3$Time = c(Time, NA, NA)
d3$label_ypos = c(label_ypos, NA, NA)
#d3[11,3] <- d3[11,3] - 2*d3[2,3]
d3_LIBD = d3
(gg_Time_LIBD <- d3_LIBD %>%
    ggplot(aes(x = Method, Time_Minutes, fill = color_all)) +
    geom_bar( stat = "identity",position = "stack",width = 0.7) + 
    theme_bw() +   xlab("") +  ylab("Time (Minutes)") + 
    # geom_text(aes(x = Method,y = reorder(label_ypos), label=Time), 
    #           vjust=-0.5,size = 3.5, 
    #           fontface = "bold") + 
    scale_x_discrete(limits = Methods) +
    geom_text(aes(y = label_ypos, label=Time), vjust=-0.5, 
              color="black", size=3.5)+
    coord_trans(y = "sqrt") +
    scale_y_continuous(breaks = c(5,10,20,50,100,500,1000,2000),
                       limits = c(0, 5000))+
    scale_fill_identity() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),   
          axis.text.y= element_text(size=rel(0.8),angle = 0, hjust = 1),  
          axis.title = element_text(size=rel(1)),      
          panel.grid.minor = element_blank(),         
          panel.grid.major.x = element_blank(),       
          panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),        
          aspect.ratio = 0.8,        
          legend.position = c(0.2, 0.85),
          legend.text=element_text(size=10)
    ) ) #+
  #guides(pattern = guide_legend(title = "",label.position = "bottom",
         #                       label.hjust = 0.5, label.vjust = 1, override.aes = list(fill = "white"), order = 2),
         #fill = guide_legend(override.aes = list(pattern = c("none", order = 1)))) 

############################ For melanoma ###################################
data = 'melanoma'
sample_names = c("mel1_rep1", "mel1_rep2",
                 "mel2_rep1", "mel2_rep2",
                 "mel3_rep1", "mel3_rep2",
                 "mel4_rep1", "mel4_rep2")
computational_cost <- read_csv("./Real/melanoma/computational_cost.csv")
setwd("./Real/melanoma")
df <- computational_cost
colnames(df) <- c("times","SPARK-X","BayesSpace_edgeR","nnSVG","SPARK","MERINGUE")
## SPARK wrong time
(melanoma_files <- list.files(pattern = "*SpaGCN_runtime.csv$", full.names = TRUE))
SpaGCN <- do.call(rbind,lapply(melanoma_files,read.csv))
SpaGCN$times <- period_to_seconds(hms(gsub("0 days ", "", SpaGCN$X0)))
(melanoma_files <- list.files(pattern = "*spatialDE_runtime.csv$", full.names = TRUE))
spatialDE <- do.call(rbind,lapply(melanoma_files,read.csv))
spatialDE$times <- period_to_seconds(hms(gsub("0 days ", "", spatialDE$X0)))
(melanoma_files <- list.files(pattern = "*spatialDE2_runtime.csv$", full.names = TRUE))
spatialDE2 <- do.call(rbind,lapply(melanoma_files,read.csv))
spatialDE2$times <- period_to_seconds(hms(gsub("0 days ", "", spatialDE2$X0)))
(melanoma_files <- list.files(pattern = "*stLearn_runtime.csv$", full.names = TRUE))
stLearn <- do.call(rbind,lapply(melanoma_files,read.csv))
stLearn$times <- period_to_seconds(hms(gsub("0 days ", "", stLearn$X0)))
computational_cost_SPARK <- read_csv("computational_cost_SPARK.csv")
computational_cost_stLearn_SV <- read_csv("computational_cost_stLearn_SV.csv")
computational_cost_BayesSpace_DESpace_seperate <- read_csv("computational_cost_BayesSpace_DESpace_seperate.csv")
df <- cbind(df,SpaGCN$times,spatialDE$times,spatialDE2$times,stLearn$times+computational_cost_stLearn_SV$V1)
colnames(df) <- c("times" ,"SPARK-X"  ,"BayesSpace_DESpace","nnSVG"  ,"SPARK",
                  "MERINGUE","sample","SpaGCN","SpatialDE", "SpatialDE2",                               
                  "StLearn_DESpace")
df = as.data.frame(df)
df$SPARK = computational_cost_SPARK$V1
df2 <- t(as.data.frame(colMeans(df[,c(2:6,8:11)])))
rownames(df2) <- c("Average")
df2 <- reshape::melt(df2)
df2$X1 <- as.factor(df2$X1)
df2$value <- df2$value/60
gg_data <- df2
colnames(gg_data) <- c("Sample_id", "Method", "Time_Minutes")
Time <- round(gg_data$Time_Minutes)
gg_data$Method <- factor(gg_data$Method, levels=c(methods_all), labels=c(methods_all))
## Main plots
colors = as.data.frame(colors_method)
colors$Method = rownames(colors)
dd_gg_data <- gg_data %>% left_join(colors)
Time <- round(dd_gg_data$Time_Minutes)
label_ypos = dd_gg_data$Time_Minutes
Methods <- arrange(dd_gg_data, Time_Minutes)$Method
dd_gg_data$Method_sub <- dd_gg_data$Method
## split BayesSpace_DESpace and StLearn_DESpace into 2 parts, respectively: clustering tool and DESpace parts.
dd_gg_data <- rbind(dd_gg_data, dd_gg_data[c(2,9),])
dd_gg_data[c(10,11),]$Time_Minutes <- c(mean(computational_cost_BayesSpace_DESpace_seperate$BayesSpace/60),
                                        mean(stLearn$times/60))
#dd_gg_data[c(2,9),]$Method_sub <- "DESpace"
dd_gg_data[c(2,9),]$Time_Minutes <- c(mean(computational_cost_BayesSpace_DESpace_seperate$edgeR/60),
                                      mean(computational_cost_stLearn_SV$V1/60))
dd_gg_data$Time_All <- dd_gg_data$Time_Minutes
dd_gg_data[c(10,11),]$Method_sub <- c("BayesSpace", "StLearn")
color_all <- c(colors_method, colours_CompuCost)
colors = as.data.frame(color_all)
colors$Method_sub = rownames(colors)
dd_data <- rbind(dd_gg_data) %>% left_join(colors)
d3 <- dd_data %>%
  mutate(
    color_all = factor(color_all, levels = c(color_bar)))
`%notin%` <- Negate(`%in%`)
d3$Cluster <- as.factor(ifelse(d3$Method_sub %in% c('BayesSpace','StLearn'), 'Cluster', 'Others'))  
d3 %>% mutate(highligth = case_when(d3$Method_sub %in% c('BayesSpace','StLearn') ~ "yes",
                                    TRUE ~ "no"))
#remotes::install_github("coolbutuseless/ggpattern")
d3$Time = c(Time, NA, NA)
d3$label_ypos = c(label_ypos, NA, NA)
d3_melanoma <- d3
(gg_Time_melanoma <- d3_melanoma %>%
    ggplot(aes(x = Method, Time_Minutes, fill = color_all)) +
    geom_bar( stat = "identity",position = "stack",width = 0.7) +
    theme_bw() +   xlab("") +  ylab("Time (Minutes)") +
    # geom_text(aes(x = Method,y = reorder(label_ypos), label=Time),
    #           vjust=-0.5,size = 3.5,
    #           fontface = "bold") +
    scale_x_discrete(limits = Methods) +
    geom_text(aes(y = label_ypos, label=Time), vjust=-0.5,
              color="black", size=3.5)+
    coord_trans(y = "sqrt") +
    scale_y_continuous(breaks = c(5,10,20,50,100,200),
                       limits = c(0, 250))+
    #scale_y_sqrt( breaks = c(5,8,10,20,50,100,200) ) +
    scale_fill_identity() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y= element_text(size=rel(0.8),angle = 0, hjust = 1),
          axis.title = element_text(size=rel(1)),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),
          aspect.ratio = 0.8,
          legend.position = c(0.2, 0.85),
          legend.text=element_text(size=10)
    ) ) #+
  #guides(pattern = guide_legend(title = "",label.position = "bottom",
       #                         label.hjust = 0.5, label.vjust = 1, override.aes = list(fill = "white"), order = 2),
       #  fill = guide_legend(override.aes = list(pattern = c("none", order = 1))))


############################ For mouse ###################################
data = 'mouse cerebellum'
computational_cost1 <- read_csv("./Real/mouse_cerebellum/computational_cost1.csv")
computational_cost1$...1 <- c("none","BayesSpace_DESpace","none","SPARK-X", "MERINGUE")
computational_cost2 <- read_csv("./Real/mouse_cerebellum/computational_cost2.csv")
df <- as.data.frame(computational_cost2)
computational_cost_stLearn_edgeR <- read_csv("./Real/mouse_cerebellum/slideseq2/computational_cost_stLearn_edgeR.csv")
colnames(df)[1] <- "times"
df$BayesSpace_DESpace <- as.numeric(computational_cost1[2,2])
df$`SPARK-X` <- as.numeric(computational_cost1[4,2])
df$MERINGUE <- as.numeric(computational_cost1[5,2])
setwd("./Real/mouse_cerebellum")
(slideseq2_files <- list.files(pattern = "*spatialDE_runtime.csv$", full.names = TRUE))
spatialDE <- do.call(rbind,lapply(slideseq2_files,read.csv))
spatialDE$times <- period_to_seconds(hms(gsub("0 days ", "", spatialDE$X0)))
(slideseq2_files <- list.files(pattern = "*spatialDE2_runtime.csv$", full.names = TRUE))
spatialDE2 <- do.call(rbind,lapply(slideseq2_files,read.csv))
spatialDE2$times <- period_to_seconds(hms(gsub("0 days ", "", spatialDE2$X0)))
(slideseq2_files <- list.files(pattern = "*stLearn_runtime.csv$", full.names = TRUE))
stLearn <- do.call(rbind,lapply(slideseq2_files,read.csv))
stLearn$times <- period_to_seconds(hms(gsub("0 days ", "", stLearn$X0)))
df$StLearn_DESpace<- as.numeric(stLearn$times+computational_cost_stLearn_edgeR$`a[3]`)

df <- cbind(df,spatialDE$times,spatialDE2$times)
df = df[-c(3,4)]
colnames(df) <- c("times" ,"nnSVG","SPARK-X",
                  "MERINGUE","SPARK","sample","BayesSpace_DESpace",
                  "StLearn_DESpace","SpatialDE",                                 
                  "SpatialDE2")
df = as.data.frame(df)
df2 <- t(as.data.frame(colMeans(df[,c(2:5,7:10)])))
rownames(df2) <- c("Average")
df2 <- reshape::melt(df2)
df2$X1 <- as.factor(df2$X1)
df2$value <- df2$value/60
gg_data <- df2
colnames(gg_data) <- c("Sample_id", "Method", "Time_Minutes")
gg_data$Method <- factor(gg_data$Method, levels=c(methods_all), labels=c(methods_all))
## Main plots
colors = as.data.frame(colors_method)
colors$Method = rownames(colors)
dd_gg_data <- gg_data %>% left_join(colors)
dd_gg_data[c(5,6),]$Time_Minutes <- c(sum(computational_cost_BayesSpace_DESpace_seperate[,2])/60,
                                      (stLearn$times+computational_cost_stLearn_edgeR$`a[3]`)/60)
Methods <- arrange(dd_gg_data, Time_Minutes)$Method
Time <- round(dd_gg_data$Time_Minutes)
label_ypos = dd_gg_data$Time_Minutes
dd_gg_data$Method_sub <- dd_gg_data$Method
## split BayesSpace_DESpace and StLearn_DESpace into 2 parts, respectively: clustering tool and DESpace parts.
dd_gg_data <- rbind(dd_gg_data, dd_gg_data[c(5,6),])
computational_cost_BayesSpace_DESpace_seperate <- as.data.frame(read_csv("computational_cost_BayesSpace_DESpace_seperate.csv"))
dd_gg_data[c(9,10),]$Time_Minutes <- c(computational_cost_BayesSpace_DESpace_seperate[1,2]/60,
                                       stLearn$times/60)
#dd_gg_data[c(5,6),]$Method_sub <- "DESpace"
dd_gg_data[c(5,6),]$Time_Minutes <- c(computational_cost_BayesSpace_DESpace_seperate[2,2]/60,
                                      computational_cost_stLearn_edgeR$`a[3]`/60)
dd_gg_data[c(9,10),]$Method_sub <- c("BayesSpace","StLearn")
dd_gg_data$Time_All <- dd_gg_data$Time_Minutes
color_all <- c(colors_method, colours_CompuCost)
colors = as.data.frame(color_all)
colors$Method_sub = rownames(colors)
dd_data <- rbind(dd_gg_data) %>% left_join(colors)
d3 <- dd_data %>%
  mutate(
    color_all = factor(color_all, levels = c(color_bar)))
`%notin%` <- Negate(`%in%`)
d3$Cluster <- as.factor(ifelse(d3$Method_sub %in% c('BayesSpace','StLearn'), 'Cluster', 'Otheres'))
d3 %>% mutate(highligth = case_when(d3$Method_sub %in% c('BayesSpace','StLearn') ~ "yes",
                                    TRUE ~ "no"))
#remotes::install_github("coolbutuseless/ggpattern")
label_ypos[5] <- d3$Time_Minutes[5] + d3$Time_Minutes[9]
label_ypos[6] <- d3$Time_Minutes[6] + d3$Time_Minutes[10]

Time[5] <- round(d3$Time_Minutes[5] + d3$Time_Minutes[9])
Time[6] <- round(d3$Time_Minutes[6] + d3$Time_Minutes[10])
d3$Time = c(Time, NA, NA)
d3$label_ypos = c(label_ypos, NA, NA)
d3_mouse <- d3
(gg_Time_mouse <- d3_mouse %>%
    ggplot(aes(x = Method, Time_Minutes, fill = color_all)) +
    geom_bar( stat = "identity",position = "stack",width = 0.7) +
    #geom_bar_pattern(stat='identity', position = "stack")+
    theme_bw() +   xlab("") +  ylab("Time (Minutes)") +
    # geom_text(aes(x = Method,y = reorder(label_ypos), label=Time),
    #           vjust=-0.5,size = 3.5,
    #           fontface = "bold") +
    scale_x_discrete(limits = Methods) +
    geom_text(aes(y = label_ypos, label=Time), vjust=-0.5,
              color="black", size=3.5)+
    coord_trans(y = "sqrt") +
    scale_y_continuous(breaks = c(30,100,1000, 20000),
                       limits = c(0, 21000))+
    #scale_y_sqrt( breaks = c(5,8,10,20,50,100,200) ) +
    scale_fill_identity() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y= element_text(size=rel(0.8),angle = 0, hjust = 1),
          axis.title = element_text(size=rel(1)),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),
          aspect.ratio = 0.8,
          legend.position = c(0.2, 0.7),#c(0.2, 0.6),
          legend.text=element_text(size=10),
          #legend.key.width=unit(2, "cm"),
          legend.box="horizontal",
         # legend.key.height=unit(1.25,"line"),
          legend.margin=margin()
    ))


(AA = egg::ggarrange( plots = list(gg_Time_LIBD + labs(title = "LIBD") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)),
                                   gg_Time_melanoma +  labs(title = "melanoma") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)),
                                   gg_Time_mouse +  labs(title = "mouse cerebellum")+ theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12))
),
#bottom = ggpubr::get_legend(gg_legend),
ncol = 3, nrow = 1))
ggsave(filename = paste0(path_save,"Barplot_Computational_Cost.pdf"),
       plot = AA,
       device = "pdf",
       width = 10,
       height = 3.8,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
