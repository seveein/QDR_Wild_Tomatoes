# LYCOPERSICUM ALL TOGEHTER 

rm(list=ls())

library(tidyverse)

rep_1 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_06_14_slycop_sclero_ann_rep1.csv")
rep_2 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_06_14_slycop_sclero_ann_rep2.csv")
rep_3 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_06_14_slycop_sclero_ann_rep3.csv")
rep_4 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_06_14_slycop_sclero_ann_rep4.csv")

data_merge <- rep_1 %>% 
  bind_rows(rep_2, rep_3, rep_4)

write.csv(x= data_merge, file= "/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_06_19_slycop_sclero_ann_ALL.csv")


(p1 <- data_merge %>% 
    filter(comm == "ok") %>%
    group_by(genotype) %>% 
    mutate(mean_ldt=mean(LDT, na.rm = T)) %>% 
    mutate(mean_lag=mean(lag, na.rm = T)) %>% 
    ungroup() %>% 
    ggplot(aes(x=reorder(factor(genotype), mean_lag), y=lag, fill=genotype))+
    geom_boxplot()+
    #geom_jitter(alpha=.6, height = 0)+
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

order <- data_merge %>% 
  filter(comm == "ok") %>%
  group_by(genotype) %>% 
  mutate(mean_lag=mean(lag, na.rm = T)) %>%
  ungroup() %>% 
  dplyr::select(c("genotype", "mean_lag")) %>% 
  unique()
order_order <- levels(reorder(order$genotype, order$mean_lag))
levels(order_order)

(p2 <- data_merge %>% 
    filter(comm == "ok") %>%
    group_by(genotype) %>% 
    mutate(mean_ldt=mean(LDT, na.rm = T)) %>% 
    mutate(mean_lag=mean(lag, na.rm = T)) %>% 
    ungroup() %>% 
    ggplot(aes(x=factor(genotype, levels=order_order), y=LDT, fill=genotype))+
    geom_boxplot()+
    #geom_jitter(alpha=.6, height = 0)+
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

(p3 <- data_merge %>% 
    filter(!(genotype=="LA4322" & start_date=="2023_04_17")) %>% 
    filter(comm!="s") %>%
    mutate(inf= ifelse(comm == "ok", 1,
                       ifelse(comm=="f", 0, NA))) %>%
    group_by(genotype, start_date) %>%
    mutate(n=n()) %>% 
    mutate(sum_inf=sum(inf)) %>% 
    mutate(ratio=sum_inf/n) %>%
    ungroup() %>% 
    group_by(genotype) %>% 
    mutate(mean=mean(ratio),
           sd= sd(ratio)) %>% 
    ggplot(aes(x=factor(genotype, levels=order_order), y=mean, fill=genotype))+
    
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.5)+
    geom_point(size=6, shape=23, show.legend = F)+
    ylim(0,1.1)+
    ylab("Infection frequency \n [%]")+
    xlab("")+
    theme_bw()+
    ggpubr::stat_compare_means(size=5)+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 00, l = 0)),
          legend.text = element_text(size=11)
    )
)
library(ggpubr)

(p_tot <- ggarrange(p1, p2, p3, nrow = 3, common.legend = T))

ggsave("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/figures/2023_SolanumScleroReps/2023_06_15slycop_screenALL.png", p_tot, dpi=600,
       height=6, width=8, bg = "white")

#COUNTS


rep_1 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/lycopersicoides/2023_08_16_slycop_sclero_COUNTS_ann_wthMOCK_rep1.csv")
rep_2 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/lycopersicoides/2023_08_16_slycop_sclero_COUNTS_ann_wthMOCK_rep2.csv")
rep_3 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/lycopersicoides/2023_08_16_slycop_sclero_COUNTS_ann_wthMOCK_rep3.csv")
rep_4 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/lycopersicoides/2023_08_16_slycop_sclero_COUNTS_ann_wthMOCK_rep4.csv")
#rep_5 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_03_24_spen_sclero_COUNTS_ann_rep5.csv")

data_merge <- rep_1 %>% 
  bind_rows(rep_2, rep_3, rep_4)

write.csv(x= data_merge, file= "/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/lycopersicoides/2023_08_16_slycop_sclero_COUNTS_wthMOCK_ALL.csv")

## REMOVE CONTAMINANTS
rm(list=ls())
data_slope <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_06_19_slycop_sclero_ann_ALL.csv", header=T)
data_slope$species[data_slope$species=="lycopersicoides"]="S. lycopersicoides"


data_count <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_06_14_slycop_sclero_COUNTS_ALL.csv", header=T)

exclude <- read.csv2("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/2023_SolSclerSkip.csv", header=T)
exclude$species[exclude$species=="S.pimpinellifolium"]="S. pimpinellifolium"
exclude$species[exclude$species=="S.lycopersicoides"]="S. lycopersicoides"

data_cleaned_slope <- data_slope %>% 
  filter(!(species%in%exclude$species & start_date %in% exclude$experiment & genotype%in%exclude$exclude))  

write.csv(x= data_cleaned_slope, file= "/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_06_19_slycop_sclero_ann_ALL_filtered.csv")

data_count$species[data_count$species=="lycopersicoides"]="S. lycopersicoides"

data_cleaned_counts <- data_count %>% 
  filter(!(species%in%exclude$species & start_date %in% exclude$experiment & genotype%in%exclude$exclude))

write.csv(x= data_cleaned_counts, file= "/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_06_14_slycop_sclero_COUNTS_ALL_filtered.csv")

