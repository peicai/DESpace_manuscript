rm(list = ls())
library(ggplot2)
library(dplyr)
library(forcats)
library(readr)
library(lubridate)
library(ggpattern)
path = "./Real_data/results/"
path_save = "./Figures/Figures/Supplementary/"
source("./Figures/Scripts/All_methods.R")
source("./Figures/Scripts/Plots_Functions.R")
color_bar = c("#A63603","#A6CEE3","#1F78B4","#08306B","#FFD92F",
              "#A6761D","#F4A582","#B2DF8A",
              "#33A02C","#BABABA","#C51B7D","plum2","mediumpurple","orchid",
              "turquoise","turquoise4", "aquamarine")
colors_method2 <- c(colors_method[1:10], "orchid",colors_method[11:12],"aquamarine",colors_method[13:14])
names(colors_method2)[c(11,14)] <- c("Manual_findMarkers", "Manual_FindAllMarkers")
#"#666666"
############################ For LIBD ###################################
data = 'LIBD'
d3_LIBD <- read.csv("./Real_data/results/LIBD/computation_time_Final.csv")
Methods <- arrange(d3_LIBD[1:16,], label_ypos)$Method
ord <- c("#A63603" ,unique(d3_LIBD$color_all)[1:16])
d3_LIBD_ord <- d3_LIBD[c(17:22, 1:16), ]
(gg_Time_LIBD <- d3_LIBD_ord %>%
    ggplot(aes(x = Method, Time_Minutes, fill = factor(color_all, levels = ord))) +
    geom_bar( stat = "identity",position = "stack",width = 0.7,
              fill = factor(d3_LIBD_ord$color_all, levels = ord)) + 
    theme_bw() +   xlab("") +  ylab("Time (Minutes)") + 
    # geom_text(aes(x = Method,y = reorder(label_ypos), label=Time), 
    #           vjust=-0.5,size = 3.5, 
    #           fontface = "bold") + 
    scale_x_discrete(limits = Methods) +
    geom_text(aes(y = label_ypos, label=Time), vjust= 0.5, hjust = -0.3, 
              color="black", size=3.5)+
    coord_trans(y = "sqrt") +
    scale_y_continuous(breaks = c(5,10,20,50,100,500,1000,2000,4000),
                       limits = c(0, 5000))+
    # scale_y_sqrt(breaks = c(20,50,100,500,1000,2000,4000),
    #                                limits = c(0, 5000))+
    #scale_fill_identity() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),   
          axis.text.y= element_text(size=rel(0.8),angle = 0, hjust = 1),  
          axis.title = element_text(size=rel(1)),      
          panel.grid.minor = element_blank(),         
          panel.grid.major.y = element_blank(),       
          panel.grid.major.x = element_line(size = 0.2, color = "grey", linetype = 2),        
          aspect.ratio = 0.5,        
          legend.position = c(0.2, 0.85),
          legend.text=element_text(size=10)
    ) )


#################### melanoma ################
data = 'melanoma'
d3_melanoma <- read.csv("./Real_data/results/melanoma/computation_time_Final.csv")
Methods <- arrange(d3_melanoma[1:13,], label_ypos)$Method
ord <- c("#A63603" ,unique(d3_melanoma$color_all)[1:13])
d3_melanoma_ord <- d3_melanoma[c(14:19, 1:13), ]
(gg_Time_melanoma <- d3_melanoma_ord %>%
    ggplot(aes(x = Method, Time_Minutes, fill = color_all)) +
    geom_bar( stat = "identity",position = "stack",width = 0.7,
              fill = factor(d3_melanoma_ord$color_all, levels = ord)) + 
    theme_bw() +   xlab("") +  ylab("Time (Minutes)")+
    # geom_text(aes(x = Method,y = reorder(label_ypos), label=Time),
    #           vjust=-0.5,size = 3.5,
    #           fontface = "bold") +
    scale_x_discrete(limits = Methods) +
    geom_text(aes(y = label_ypos, label=Time), vjust= 0.5, hjust = -0.3,
              color="black", size=3.5)+
     coord_trans(y = "sqrt") +
     scale_y_continuous(breaks = c(5,10,20,50,100,200),
                        limits = c(0, 250))+
    #scale_y_sqrt(breaks = c(5,10,20,50,100,200),
    #             limits = c(0, 250))+
    #scale_fill_identity() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y= element_text(size=rel(0.8),angle = 0, hjust = 1),
          axis.title = element_text(size=rel(1)),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(size = 0.2, color = "grey", linetype = 2),
          aspect.ratio = 0.5,
          legend.position = c(0.2, 0.85),
          legend.text=element_text(size=10)
    ) ) 
#################### mouse cerebellum ################
data = 'slideseq2'
d3_mouse <- read.csv("./Real_data/results/slideseq2/computation_time_Final.csv")
Methods <- arrange(d3_mouse[1:12,], label_ypos)$Method
ord <- c("#A63603" ,unique(d3_mouse$color_all)[1:12])
d3_mouse_ord <- d3_mouse[c(13:18, 1:12), ]
(gg_Time_mouse <- d3_mouse_ord %>%
    ggplot(aes(x = Method, Time_Minutes, fill = color_all)) +
    geom_bar( stat = "identity",position = "stack",width = 0.7,
              fill = factor(d3_mouse_ord$color_all, levels = ord)) + 
    #geom_bar_pattern(stat='identity', position = "stack")+
    theme_bw() +   xlab("") +  ylab("Time (Minutes)") +
    # geom_text(aes(x = Method,y = reorder(label_ypos), label=Time),
    #           vjust=-0.5,size = 3.5,
    #           fontface = "bold") +
    scale_x_discrete(limits = Methods) +
    geom_text(aes(y = label_ypos, label=Time), vjust= 0.5, hjust = -0.3,
              color="black", size=3.5)+
     coord_trans(y = "sqrt") +
     scale_y_continuous(breaks = c(50, 100,1000, 20000),
                        limits = c(0, 21000))+
    #scale_y_sqrt(breaks = c(100,1000, 20000),
    #             limits = c(0, 21000))+
    #scale_y_sqrt( breaks = c(5,8,10,20,50,100,200) ) +
    #scale_fill_identity() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y= element_text(size=rel(0.8),angle = 0, hjust = 1),
          axis.title = element_text(size=rel(1)),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(size = 0.2, color = "grey", linetype = 2),
          aspect.ratio = 0.5,
          legend.position = c(0.2, 0.7),#c(0.2, 0.6),
          legend.text=element_text(size=10),
          #legend.key.width=unit(2, "cm"),
          legend.box="horizontal",
          # legend.key.height=unit(1.25,"line"),
          legend.margin=margin()
    ))
library(deeptime)
(AA = egg::ggarrange( plots = list(gg_Time_LIBD +coord_trans_flip(y = "sqrt")
                                   + labs(title = "LIBD") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)),
                                   gg_Time_melanoma +coord_trans_flip(y = "sqrt")
                                   +  labs(title = "melanoma") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12)),
                                   gg_Time_mouse +coord_trans_flip(y = "sqrt")
                                   +  labs(title = "mouse cerebellum")+ theme(plot.title = element_text(hjust = 0.5, face = "bold", size=12))
),
#bottom = ggpubr::get_legend(gg_legend),
ncol = 1, nrow = 3))
ggsave(filename = paste0(path_save,"Barplot_Computational_Cost.pdf"),
       plot = AA,
       device = "pdf",
       width = 8,
       height = 11,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

