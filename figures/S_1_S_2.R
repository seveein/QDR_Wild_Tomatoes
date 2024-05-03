# Focus on the measures: Lag, ldt etc- but comparing species | s. habro | s. pimp
# suppl. Fig 1 + 2
# Severin Einspanier
# 2024_04_25

rm(list=ls())
library(tidyverse)
library(ggpubr)
library(ochRe)
data_slopes  <- read.csv("C://Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/AllSlopes_Ann_wthMock_noContaminants_threeReps_filtered.csv")


(p1 <- data_slopes %>% 
    filter(comm == "ok"&
             genotype!="C32") %>%
    mutate(lag_days = lag/(60*24)) %>% 
    group_by(species) %>% 
    mutate(mean_ldt=mean(LDT, na.rm = T),
           mean_lag=mean(lag_days, na.rm = T),
           n=n()) %>% 
    ungroup() %>% 
    summarise(species=species,
              genotype=genotype,
              lag=lag,
              lag_days=lag_days,
              LDT=LDT,
              SampleID=SampleID,
              start_date=start_date,
              mean_lag=mean_lag,
              mean_ldt = mean_ldt,
              n=n) %>% 
    unique() %>% 
    ggplot(aes(x=reorder(factor(species), mean_lag), y=lag_days, fill=species))+
    geom_boxplot(show.legend = F, size=1, col="black")+
    #geom_jitter(alpha=.6, height = 0)+
    #geom_point(aes(y=mean_lag), size=3, col="Red", show.legend = F)+
    #geom_text(aes(y=0, label=n))+
    xlab("")+
    ylab("Lag phase \n [days]")+
    scale_y_continuous(breaks = c(2,4,6))+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(size=11, angle=45,
                                     hjust=1, vjust=1, color="black",
                                     face="italic"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm')
    )+
    scale_fill_ochre(palette="parliament") +
    guides(fill=guide_legend(title="Plant Species",
                             label.theme = element_text(angle = 0, face = "italic")))
)


cols_pen <- c("#2061d0", "#3c75d6", "#5889dc", "#749ce2", "#90b0e8", "#abc4ed", "#c7d8f3","#e3ebf9", "#ffffff" )
cols_pen <- c("#266F6F", "#2C8484", "#339A9A", "#39B0B0", "#43C3C3", "#58CACA", "#5ACBCB", "#8ED8DB", "#ffffff")
cols_pen <- c("#787249", "#948c5a", "#a9a272", "#bbb48e", "#ccc7ab", "#dddac7", "#eeece3","#ffffff")


(p2 <- data_slopes %>% 
    filter(comm == "ok", 
           species=="S. habrochaites") %>%
    mutate(lag_days = lag/(60*24)) %>% 
    group_by(genotype) %>% 
    mutate(mean_ldt=mean(LDT, na.rm = T),
           mean_lag=mean(lag_days, na.rm = T),
           n=n()) %>% 
    ungroup() %>% 
    summarise(species=species,
              genotype=genotype,
              lag=lag,
              lag_days=lag_days,
              LDT=LDT,
              SampleID=SampleID,
              start_date=start_date,
              mean_lag=mean_lag,
              mean_ldt = mean_ldt,
              n=n) %>% 
    unique() %>% 
    ggplot(aes(x=reorder(factor(genotype), -mean_lag), y=lag_days, fill=genotype))+
    geom_boxplot(show.legend = F, size=1, col="black",
                 fill=cols_pen)+
    #geom_jitter(alpha=.6, height = 0)+
    #geom_point(aes(y=mean_lag), size=3, col="Red", show.legend = F)+
    geom_text(aes(y=0, label=n))+
    xlab("")+
    ylab("")+
    scale_y_continuous(limits=c(0,6), breaks = c(2,4,6))+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(angle=45, vjust=1, 
                                     hjust=1, size=11, color="black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm')
    )+
    
    #    geom_text(data=groups, aes(x=genotype, label=groups), y=7)+
    
    guides(fill=guide_legend(title="Plant Species",
                             label.theme = element_text(angle = 0, face = "italic"))
    )
)


# VIZ Pimpi

cols_lyco = c("#9b8a20",  "#b4a758", "#c1b674", "#cdc590",  "#e6e2c7", "#f3f0e3", "#ffffff")
cols_lyco = c("#287235", "#369A5C", "#3EB66C", "#53C57F", "#6ECF93", "#89D9A7", "#ffffff")
cols_lyco = c("#851e21", "#a52528", "#c42c30", "#d44346", "#dc6265", "#e38184", "#eaa1a3","#f1c0c1", "#f8e0e0", "#ffffff")

(p3 <- data_slopes %>% 
    filter(comm == "ok", 
           species=="S. pimpinellifolium") %>%
    mutate(lag_days = lag/(60*24)) %>% 
    group_by(genotype) %>% 
    mutate(mean_ldt=mean(LDT, na.rm = T),
           mean_lag=mean(lag_days, na.rm = T),
           n=n()) %>% 
    ungroup() %>% 
    summarise(species=species,
              genotype=genotype,
              lag=lag,
              lag_days=lag_days,
              LDT=LDT,
              SampleID=SampleID,
              start_date=start_date,
              mean_lag=mean_lag,
              mean_ldt = mean_ldt,
              n=n) %>% 
    unique() %>% 
    #arrange(mean_lag) %>% 
    ggplot(aes(x=reorder(factor(genotype), -mean_lag), y=lag_days, fill=genotype))+
    geom_boxplot(show.legend = F, size=1, col="black",
                 fill=cols_lyco)+
    #geom_jitter(alpha=.6, height = 0)+
    #geom_point(aes(y=mean_lag), size=3, col="Red", show.legend = F)+
    geom_text(aes(y=0, label=n))+
    xlab("")+
    ylab("")+
    scale_y_continuous(breaks = c(2,4,6), limits=c(0,6))+
    #scale_fill_manual(values=cols_pen)+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(size=11, angle=45,
                                     hjust=1, vjust=1, color="black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm')
    )+
    #geom_text(data=groups,aes(x=genotype, label=groups), y=5)+
    
    guides(fill=guide_legend(title="Plant Species",
                             label.theme = element_text(angle = 0, face = "italic")))
)


(p_tot <- ggarrange(p1, p2, p3, nrow=1,
                    labels = c("A", "B", "C"),
                    align = "hv")
)

ggsave("C://Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/figures/S_1.png", p_tot, dpi=600,
       width=12, height=4)


#### LDT


(p1 <- data_slopes %>% 
    filter(comm == "ok" & 
             genotype != "C32") %>%
    mutate(lag_days = lag/(60*24)) %>% 
    group_by(species) %>% 
    mutate(mean_ldt=mean(LDT, na.rm = T),
           mean_lag=mean(lag_days, na.rm = T),
           n=n()) %>% 
    ungroup() %>% 
    summarise(species=species,
              genotype=genotype,
              lag=lag,
              lag_days=lag_days,
              LDT=LDT,
              SampleID=SampleID,
              start_date=start_date,
              mean_lag=mean_lag,
              mean_ldt = mean_ldt,
              n=n) %>% 
    unique() %>% 
    ggplot(aes(x=reorder(factor(species), mean_ldt), y=LDT, fill=species))+
    geom_boxplot(show.legend = F, size=1)+
    #geom_jitter(alpha=.6, height = 0)+
    #geom_point(aes(y=mean_lag), size=3, col="Red", show.legend = F)+
    #geom_text(aes(y=4.5, label=n))+
    xlab("")+
    ylab("log(LDT)")+
    #scale_y_continuous(breaks = c(1,2,3,4,5,6))+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(size=11, angle=45,
                                     hjust=1, vjust=1, color="black",
                                     face="italic"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm')
    )+
    
    scale_fill_ochre(palette="parliament") +
    guides(fill=guide_legend(title="Plant Species",
                             label.theme = element_text(angle = 0, face = "italic")))+
    scale_y_continuous(limit=c(4,10.7), breaks=c(4.094345, 5.192957, 6.291569,6.984716,7.677864, 8.371011, 9.064158,10.6),labels = c(1, 3,9, 18, 36, 72, 144, 720))
  
)

(p2 <- data_slopes %>% 
    filter(comm == "ok", 
           species=="S. habrochaites") %>%
    mutate(lag_days = lag/(60*24)) %>% 
    group_by(genotype) %>% 
    mutate(mean_ldt=mean(LDT, na.rm = T),
           mean_lag=mean(lag_days, na.rm = T),
           n=n()) %>% 
    ungroup() %>% 
    summarise(species=species,
              genotype=genotype,
              lag=lag,
              lag_days=lag_days,
              LDT=LDT,
              SampleID=SampleID,
              start_date=start_date,
              mean_lag=mean_lag,
              mean_ldt = mean_ldt,
              n=n) %>% 
    unique() %>% 
    ggplot(aes(x=reorder(factor(genotype), -mean_ldt), y=LDT, fill=genotype))+
    geom_boxplot(show.legend = F, size=1,
                 fill=cols_pen)+
    #geom_jitter(alpha=.6, height = 0)+
    #geom_point(aes(y=mean_lag), size=3, col="Red", show.legend = F)+
    geom_text(aes(y=4.5, label=n))+
    xlab("")+
    ylab("")+
    ylim(c(4.5, 10))+
    #scale_y_continuous(breaks = c(1,2,3,4,5,6))+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(angle=45, vjust=1, 
                                     hjust=1, size=11, color="black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm')
    )+
    
    guides(fill=guide_legend(title="Plant Species",
                             label.theme = element_text(angle = 0, face = "italic")))+
    scale_y_continuous(limit=c(4,10.7), breaks=c(4.094345, 5.192957, 6.291569,6.984716,7.677864, 8.371011, 9.064158,10.6),labels = c(1, 3,9, 18, 36, 72, 144, 720))
  
  )

(p3 <- data_slopes %>% 
    filter(comm == "ok", 
           species=="S. pimpinellifolium") %>%
    mutate(lag_days = lag/(60*24)) %>% 
    group_by(genotype) %>% 
    mutate(mean_ldt=mean(LDT, na.rm = T),
           mean_lag=mean(lag_days, na.rm = T),
           n=n()) %>% 
    ungroup() %>% 
    summarise(species=species,
              genotype=genotype,
              lag=lag,
              lag_days=lag_days,
              LDT=LDT,
              SampleID=SampleID,
              start_date=start_date,
              mean_lag=mean_lag,
              mean_ldt = mean_ldt,
              n=n) %>% 
    unique() %>% 
    ggplot(aes(x=reorder(factor(genotype), -mean_ldt), y=LDT, fill=genotype))+
    geom_boxplot(show.legend = F, size=1, 
                 fill=cols_lyco)+
    #geom_jitter(alpha=.6, height = 0)+
    #geom_point(aes(y=mean_lag), size=3, col="Red", show.legend = F)+
    geom_text(aes(y=4.5, label=n))+
    xlab("")+
    ylab("")+
    #scale_y_continuous(breaks = c(1,2,3,4,5,6))+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(size=11, angle=45,
                                     hjust=1, vjust=1, color="black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm')
    )+
    
    guides(fill=guide_legend(title="Plant Species",
                             label.theme = element_text(angle = 0, face = "italic")))+
    scale_y_continuous(limit=c(4,10.7), breaks=c(4.094345, 5.192957, 6.291569,6.984716,7.677864, 8.371011, 9.064158,10.6),labels = c(1, 3,9, 18, 36, 72, 144, 720))
  
)


(p_tot <- ggarrange(p1, p2, p3, nrow=1,
                    labels = c("A", "B", "C"),
                    align = "hv"
))

ggsave("C://Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/figures/S_2.png", p_tot, dpi=600,
       width=12, height=4)
