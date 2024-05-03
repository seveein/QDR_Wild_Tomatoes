# Fith HABROCHAITES SCREENING


# OVERVIEW SCRIPT FOR Slopes - ID'ing and rought figures
rm(list=ls())
library(tidyverse)
library(ggpubr)

data_0 <- read.csv("/mnt/PHYTOP/SEVERIN/navautron/2023_05_30/0/slopes_1.txt", sep="\t", header=F)

data_0 <- data_0 %>% 
  mutate(Box = 1)
#data_0$V1 <- as.integer(data_0$V1)
#data_0$V3 <- as.numeric(data_0$V3)

data_2 <- read.csv("/mnt/PHYTOP/SEVERIN/navautron/2023_05_30/2/slopes_1.txt", sep="\t", header=F)

data_2 <- data_2 %>% 
  mutate(Box = 2)

data_4 <- read.csv("/mnt/PHYTOP/SEVERIN/navautron/2023_05_30/4/slopes_1.txt", sep="\t", header=F)

data_4 <- data_4 %>% 
  mutate(Box = 3)

#data_6 <- read.csv

id <- read.csv2("/mnt/PHYTOP/SEVERIN/navautron/2023_05_30/layout_scheme.csv", header=T)

head(id)

data_merge <- data_0 %>% 
  bind_rows(data_2, data_4) 

colnames(data_merge)=c("ID", "lag", "slope", "comm", "LDT", "Box")

data_merge_id <- id %>% 
  left_join(data_merge, by=c("Box", "ID"))

write.csv( data_merge_id, "/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_03_24_shabro_sclero_ann_rep5.csv")

(p1 <- data_merge_id %>% 
    filter(comm == "ok") %>%
    group_by(genotype) %>% 
    mutate(mean_ldt=mean(LDT, na.rm = T)) %>% 
    mutate(mean_lag=mean(lag, na.rm = T)) %>% 
    ungroup() %>% 
    ggplot(aes(x=reorder(factor(genotype), mean_lag), y=lag, fill=genotype))+
    geom_boxplot()+
    geom_jitter(alpha=.6, height = 0)+
    geom_point(aes(y=mean_lag), size=3, col="Red", show.legend = F)+
    xlab("")+
    ylab("Latency phase \n [min]")+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11)
    )+
    
    guides(fill=guide_legend(title="Plant Genotype"))
)

# extract order from lag phase

order <- data_merge_id %>% 
  filter(comm == "ok") %>%
  group_by(genotype) %>% 
  mutate(mean_lag=mean(lag, na.rm = T)) %>%
  ungroup() %>% 
  dplyr::select(c("genotype", "mean_lag")) %>% 
  unique()
order_order <- levels(reorder(order$genotype, order$mean_lag))
levels(order_order)

(p2 <- data_merge_id %>% 
    filter(comm == "ok") %>%
    group_by(genotype) %>% 
    mutate(mean_ldt=mean(LDT, na.rm = T)) %>% 
    mutate(mean_lag=mean(lag, na.rm = T)) %>% 
    ungroup() %>% 
    ggplot(aes(x=factor(genotype, levels=order_order), y=LDT, fill=genotype))+
    geom_boxplot()+
    geom_jitter(alpha=.6, height = 0)+
    geom_point(aes(y=mean_ldt), size=3, col="Red", show.legend = F)+
    ylab("Lesion doubling time \n [log(LDT)]")+
    xlab("")+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 00, l = 0)),
          legend.text = element_text(size=11)
    )+ 
    guides(fill=guide_legend(title="Plant Genotype"))
)


# Infection freq

(p3 <- data_merge_id %>% 
    mutate(inf= ifelse(comm == "ok", 1,
                       ifelse(comm=="f", 0, NA))) %>%
    group_by(genotype) %>%
    filter(comm!="s") %>% 
    mutate(n=n()) %>% 
    mutate(sum_inf=sum(inf)) %>% 
    mutate(ratio=sum_inf/n) %>% 
    ggplot(aes(x=factor(genotype, levels=order_order), y=ratio*100, fill=genotype))+
    geom_point(size=6, shape=23, show.legend = F)+
    ylim(0,100)+
    ylab("Infection frequency \n [%]")+
    xlab("")+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 00, l = 0)),
          legend.text = element_text(size=11)
    )
)

(p_tot <- ggarrange(p1, p2, p3, nrow = 3, common.legend = T))
ggsave("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/figures/2023_SolanumScleroReps/2023_06_21shabro_screen5.png", p_tot, dpi=600,
       height=6, width=8, bg = "white")


##### COUNTS annotation
data_0 <- read.csv("/mnt/PHYTOP/SEVERIN/navautron/2023_05_30/0/counts1.txt", sep="\t", header=T)

data_0 <- data_0 %>% 
  mutate(Box = 1)
#data_0$V1 <- as.integer(data_0$V1)
#data_0$V3 <- as.numeric(data_0$V3)

data_2 <- read.csv("/mnt/PHYTOP/SEVERIN/navautron/2023_05_30/2/counts1.txt", sep="\t", header=T)

data_2 <- data_2 %>% 
  mutate(Box = 2)

data_4 <- read.csv("/mnt/PHYTOP/SEVERIN/navautron/2023_05_30/4/counts1.txt", sep="\t", header=T)

data_4 <- data_4 %>% 
  mutate(Box = 3)

data_6 <- read.csv("/mnt/PHYTOP/SEVERIN/navautron/2023_05_30/6/counts1.txt", sep="\t", header=T)

data_6 <- data_6 %>% 
  mutate(Box = 4)

id <- read.csv2("/mnt/PHYTOP/SEVERIN/navautron/2023_05_30/layout_scheme.csv", header=T)

head(id)

data_merge <- data_0 %>% 
  bind_rows(data_2, data_4, data_6) 

#colnames(data_merge)=c("ID", "lag", "slope", "comm", "LDT", "Box")

data_merge_id <- id %>% 
  left_join(data_merge, by=c("Box", "ID"))



write.csv( data_merge_id, "/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16_shabro_sclero_COUNTS_ann_rep5_wthMOCK.csv")

