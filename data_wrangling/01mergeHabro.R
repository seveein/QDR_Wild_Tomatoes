# Habrochaites ALL TOGEHTER 

rm(list=ls())

library(tidyverse)
rep_1 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_03_24_shabro_sclero_ann_rep1.csv")
rep_2 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_03_24_shabro_sclero_ann_rep2.csv")
rep_3 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_03_24_shabro_sclero_ann_rep3.csv")
rep_4 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_03_24_shabro_sclero_ann_rep4.csv")
rep_5 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_03_24_shabro_sclero_ann_rep5.csv")

rep_1 <- read.csv("C:/Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/2023_03_24_shabro_sclero_ann_rep1.csv")
rep_2 <- read.csv("C:/Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/2023_03_24_shabro_sclero_ann_rep2.csv")
rep_3 <- read.csv("C:/Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/2023_03_24_shabro_sclero_ann_rep3.csv")
rep_4 <- read.csv("C:/Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/2023_03_24_shabro_sclero_ann_rep4.csv")
rep_5 <- read.csv("C:/Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/2023_03_24_shabro_sclero_ann_rep5.csv")



#rep_4 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_03_2023_spen_sclero_ann_rep4.csv")
#rep_5 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_06_spen_sclero_ann_rep5.csv", header=T)
data_merge <- rep_5 %>% 
  bind_rows(rep_1, rep_2, rep_3, rep_4)

write.csv(x= data_merge, file= "C:/Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/2023_06_shabro_sclero_ann_ALL.csv")

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
    filter(comm!="s") %>%
    mutate(inf= ifelse(comm == "ok", 1,
                       ifelse(comm=="f", 0, NA))) %>%
    group_by(genotype, start_date) %>%
    mutate(n=n()) %>% 
    mutate(sum_inf=sum(inf)) %>% 
    summarize(ratio=sum_inf/n) %>% 
    group_by(genotype) %>% 
    mutate(mean = mean(ratio)) %>% 
    mutate(sd = sd(ratio)) %>% 
    ggplot(aes(x=factor(genotype, levels=order_order), fill=genotype))+
    #geom_point(size=6, shape=23, show.legend = F)+
    geom_point(aes(y=mean), shape=21, size=3) +
    geom_errorbar(aes(ymax=mean+sd, ymin=mean-sd), width=.4)+
    #ylim(0,100)+
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

ggsave("/ResiDEvo/figures/2023_SolanumScleroReps/2023_06_15shabro_screenALL.png", p_tot, dpi=600,
       height=6, width=8, bg = "white")


#COUNTS


rep_1 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/habrochaites/2023_08_16_shabro_sclero_COUNTS_ann_rep1_wthMOCK.csv")
rep_2 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/habrochaites/2023_08_16_shabro_sclero_COUNTS_ann_rep2_wthMOCK.csv")
rep_3 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/habrochaites/2023_08_16_shabro_sclero_COUNTS_ann_rep3_wthMOCK.csv")
rep_4 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/habrochaites/2023_08_16_shabro_sclero_COUNTS_ann_rep4_wthMOCK.csv")
rep_5 <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/habrochaites/2023_08_16_shabro_sclero_COUNTS_ann_rep5_wthMOCK.csv")

data_merge <- rep_1 %>% 
  bind_rows(rep_2, rep_3, rep_4, rep_5)

check_mock <- data_merge %>% 
  filter(inoculum=="mock")

write.csv(x= data_merge, file= "/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/habrochaites/2023_08_16_shabro_sclero_COUNTS_ann_ALL_wthMOCK.csv")

