# Figure 3 of Solanum x Sclerotinia phenotyping paper
# Deep-dive into LAG+LDT of two focus species (S. pennellii and S. lycopersicoides)
# as well as statistical tables 


rm(list=ls())
options(scipen=999)
library(tidyverse)
library(ggpubr)
library(ochRe)

#setwd("C://Users/Sever/Nextcloud/ResiDEvo/")
#setwd("C://Users/suaph281/Nextcloud/ResiDEvo/")
#data_full <- read.csv("data/2023_ScleroPhenotypes/AllSlopes_Ann_wthMock_coContaminants_threeReps_filtered.csv")

# Viz LDT and LAG for Populations 
# Focus on the measures: Lag, ldt etc- but comparing species | s.penn | s. lycoper
### EXCLUDED THE SPECIES OVERVIEW 
### former code at the end of this script. 

rm(list=ls())
library(tidyverse)
library(ggpubr)

# import basic data for LDT and Lag
data_slopes  <- read.csv("C://Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/AllSlopes_Ann_wthMock_noContaminants_threeReps_filtered.csv")
                         
                         
#### S. pennellii ####

# get individual population estimates of lag phase duration
pop_mean_lag <- read.csv("C://Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/2023_10_09_lag_Pops_Means.csv")

cols_pen <- c("#5c7052", "#6f8762", "#829b76", "#829b76", "#97ac8d", "#acbca4", "#c1cdba","#d5ded1", "#eaeee8")

(p1_lag <- data_slopes %>% 
    filter(comm == "ok", 
           species=="S. pennellii") %>%
    mutate(lag_days = lag/(60*24)) %>% 
    group_by(genotype) %>% 
    mutate(n=n()) %>% 
    ungroup() %>% 
    left_join(pop_mean_lag, by="genotype") %>% 
    summarise(species=species,
              genotype=genotype,
              lag=lag,
              lag_days=lag_days,
              SampleID=SampleID,
              start_date=start_date,
              mean_lag=lsmean,
              n=n) %>% 
    unique() %>% 
    ggplot()+
    geom_boxplot(aes(x=reorder(factor(genotype), -mean_lag), y=lag_days, fill=genotype),
                 show.legend = F, size=1, col="black", fill=cols_pen
    )+
    geom_text(aes(x=reorder(factor(genotype), -mean_lag), y=0, label=n))+    
    xlab("")+
    ylab("Lag phase [days]")+
    scale_y_continuous(limits=c(0,6.4), breaks = c(2,4,6))+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(angle=45, vjust=1, 
                                     hjust=1, size=11, color="black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm'),
          panel.grid.major.x =  element_blank()
    )+
    guides(fill=guide_legend(title="Plant Species",
                             label.theme = element_text(angle = 0, face = "italic"))
    )
)


p1_data <- read.table("C:/Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/tables/tab4/2023_10_09_lag_CompPops.csv")
p1_means <- read.csv("C:/Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/tables/tab4/2023_10_09_lag_Pops_Means.csv")

p1_data$contrast <- as.factor(p1_data$contrast)
# shape statistics


p1_complete_table <- p1_data %>% 
  separate(col=contrast, into = c("group1", "group2"), sep = ", ", remove = F) %>% 
  filter(group1 %in% c("S. pennellii")) %>% 
  mutate(group1=factor(group1),
         group2=factor(group2)) %>% 
  dplyr::select(group1, group2, estimate, std.error, adj.p.value)  %>% 
  mutate(p.signif = symnum(adj.p.value, corr = FALSE, na = FALSE, 
                           cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", "ns"))
  ) %>% 
  summarise(species=group1,
            genotype = group2,
            estimate=estimate,
            std.error =std.error,
            adj.p.value=adj.p.value,
            p.signif = p.signif) %>% 
  left_join(p1_means, by="genotype" ) %>% 
  arrange(species) %>% 
  group_by(species) %>% 
  arrange(lsmean, .by_group = T) %>%
  mutate( adj.p.value = format.pval(pv = adj.p.value,
                                    # digits : number of digits, but after the 0.0
                                    digits = 2, 
                                    # eps = the threshold value above wich the 
                                    # function will replace the pvalue by "<0.0xxx"
                                    eps = 0.001,
                                    # nsmall = how much tails 0 to keep if digits of 
                                    # original value < to digits defined
                                    nsmall = 3),
          estimate = round(estimate/(60*24), digits=2),
          lsmean = round(lsmean/(60*24), digits=2),
          SE = round(SE/(60*24), digits=2),
          std.error = round(std.error/(60*24), digits=2)
          ) %>% 
  ungroup() %>% 
  dplyr::select(
    genotype, lsmean, SE, estimate,std.error, adj.p.value, p.signif
  )

colnames(p1_complete_table) = c("Genotype", "lsmean", "SE", "Estimate comp.", "SE comp.", "adj. p-value", "signif.")



(p1_lag_stats <- ggtexttable(p1_complete_table, theme = ttheme("blank", 
                                                  padding = unit(c(4, 4), "mm"),
                                                  base_size = 7.5), 
                   rows = NULL) %>% 
    table_cell_font( row=2:7, column = 1, face = "italic", size=7.5)%>% 
    tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 2, linetype = 1)%>% 
    tab_add_hline(at.row = c(10), row.side = "bottom", linewidth = 2, linetype = 1)
)



#### S. lycopersicoides ####

data_lycop <- data_slopes %>% 
  filter(comm == "ok", 
         species=="S. lycopersicoides")


cols_lyco = c("#537e8c", "#6997a6", "#87acb8", "#a5c1ca", "#c3d5db", "#e1eaed", "#e9f0f2")

(p2_lag <- data_slopes %>% 
    filter(comm == "ok", 
           species=="S. lycopersicoides") %>%
    mutate(lag_days = lag/(60*24)) %>% 
    group_by(genotype) %>% 
    mutate(n=n()) %>% 
    ungroup() %>% 
    left_join(pop_mean_lag, by="genotype") %>% 
    summarise(species=species,
              genotype=genotype,
              lag=lag,
              lag_days=lag_days,
              SampleID=SampleID,
              start_date=start_date,
              mean_lag=lsmean,
              n=n) %>% 
    unique() %>% 
    ggplot()+
    geom_boxplot(aes(x=reorder(factor(genotype), -mean_lag), y=lag_days, fill=genotype),
                 show.legend = F, size=1, col="black",
                 fill=cols_lyco)+    
    geom_text(aes(x=reorder(factor(genotype), -mean_lag), y=0, label=n))+
    xlab("")+
    ylab("")+
    scale_y_continuous(breaks = c(2,4,6), limits=c(0,6.4))+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(size=11, angle=45,
                                     hjust=1, vjust=1, color="black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm'),
          panel.grid.major.x =  element_blank()
    )+
    
    guides(fill=guide_legend(title="Plant Species",
                             label.theme = element_text(angle = 0, face = "italic")))
)


p2_complete_table <- p1_data %>% 
  separate(col=contrast, into = c("group1", "group2"), sep = ", ", remove = F) %>% 
  filter(group1 %in% c("S.lycopersicoides")) %>% 
  mutate(group1=factor(group1),
         group2=factor(group2)) %>% 
  dplyr::select(group1, group2, estimate, std.error, adj.p.value)  %>% 
  mutate(p.signif = symnum(adj.p.value, corr = FALSE, na = FALSE, 
                           cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", "ns"))
  ) %>% 
  summarise(species=group1,
            genotype = group2,
            estimate=estimate,
            std.error =std.error,
            adj.p.value=adj.p.value,
            p.signif = p.signif) %>% 
  left_join(p1_means, by="genotype" ) %>% 
  arrange(species) %>% 
  group_by(species) %>% 
  arrange(lsmean, .by_group = T) %>%
  mutate( adj.p.value = format.pval(pv = adj.p.value,
                                    digits = 2, 
                                    eps = 0.001,
                                    nsmall = 3),
          estimate = round(estimate/(60*24), digits=2),
          lsmean = round(lsmean/(60*24), digits=2),
          SE = round(SE/(60*24), digits=2),
          std.error = round(std.error/(60*24), digits=2)
  ) %>% 
  ungroup() %>% 
  dplyr::select(
    genotype, lsmean, SE, estimate,std.error, adj.p.value, p.signif
  )

colnames(p2_complete_table) = c("Genotype", "lsmean", "SE", "Estimate comp.", "SE comp.", "adj. p-value", "signif.")



(p2_lag_stats <- ggtexttable(p2_complete_table, theme = ttheme("blank", 
                                                               padding = unit(c(4, 4), "mm"),
                                                               base_size = 7.5), 
                             rows = NULL) %>% 
    table_cell_font( row=2:7, column = 1, face = "italic", size=7.5)%>% 
    tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 2, linetype = 1)%>% 
    tab_add_hline(at.row = c(8), row.side = "bottom", linewidth = 2, linetype = 1)
)

#### ARRANGE plots lag ####
(p_tot_lag <- ggarrange(p1_lag, p2_lag,
                        p1_lag_stats, p2_lag_stats,
                        nrow = 2, 
                        ncol=2, 
                        heights = c(4,3,4,3),
                        labels=c("A", "B", "C", "D"),
                        align = "v"))



######################## LDT #############################


#### S. pennellii #####

#
#  IMPORTANT NOTE:Mean from individual  
# estimates in LDT is NOT correct


# import adjusted estimates of each population
pop_mean <- read.csv("C://Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/tables/tab5/2024_01_19_ldt_Populations_Means.csv")

(p3 <- data_slopes %>% 
    filter(comm == "ok", 
           species=="S. pennellii") %>%
    group_by(genotype) %>% 
    mutate(mean_ldt=mean(LDT, na.rm = T),
           n=n()) %>% 
    ungroup() %>% 
    left_join(pop_mean, by="genotype") %>% 
    summarise(species=species,
              genotype=genotype,
              LDT=LDT,
              SampleID=SampleID,
              start_date=start_date,
              mean_ldt = lsmean,
              n=n) %>% 
    unique() %>% 
    ggplot()+
    
    geom_boxplot(aes(x=reorder(factor(genotype), -mean_ldt), y=LDT, fill=genotype),
                 show.legend = F, size=1,
                 fill=cols_pen)+    
    geom_text(aes(x=reorder(factor(genotype), -mean_ldt), y=4.5, label=n))+
    xlab("")+
    ylab("Lesion Doubling Time \nLDT [h]")+
    scale_y_continuous(limit=c(4.4,10.7), 
                       breaks=c(5.192957, 6.291569,6.984716,7.677864, 8.371011, 9.064158,10.6),
                       labels = c(3,9, 18, 36, 72, 144, 720))+
    #scale_y_continuous(breaks = c(1,2,3,4,5,6))+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(angle=45, vjust=1, 
                                     hjust=1, size=11, color="black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, 
                                      margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm'),
          panel.grid.major.x =  element_blank()
    )+
    
    guides(fill=guide_legend(title="Plant Species",
                             label.theme = element_text(angle = 0, face = "italic")))
)

p3_data <- read.table("C:/Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/tables/tab5/2024_01_19_ldt_CompPops.csv")
p3_means <- read.csv("C:/Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/tables/tab5/2024_01_19_ldt_Populations_Means.csv")

p3_data$contrast <- as.factor(p3_data$contrast)
# shape statistics


p3_complete_table <- p3_data %>% 
  separate(col=contrast, into = c("group1", "group2"), sep = ", ", remove = F) %>% 
  filter(group1 %in% c("S. pennellii")) %>% 
  mutate(group1=factor(group1),
         group2=factor(group2)) %>% 
  dplyr::select(group1, group2, estimate, std.error, adj.p.value)  %>% 
  mutate(p.signif = symnum(adj.p.value, corr = FALSE, na = FALSE, 
                           cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", "ns"))
  ) %>% 
  summarise(species=group1,
            genotype = group2,
            estimate=estimate,
            std.error =std.error,
            adj.p.value=adj.p.value,
            p.signif = p.signif) %>% 
  left_join(p3_means, by="genotype" ) %>% 
  arrange(species) %>% 
  group_by(species) %>% 
  arrange(lsmean, .by_group = T) %>%
  mutate( adj.p.value = format.pval(pv = adj.p.value,
                                    digits = 2, 
                                    eps = 0.001,
                                    nsmall = 3),
          estimate = round(estimate, digits=2),
          lsmean = round(lsmean, digits=2),
          SE = round(SE, digits=2),
          std.error = round(std.error, digits=2)
  ) %>% 
  ungroup() %>% 
  dplyr::select(
    genotype, lsmean, SE, estimate,std.error, adj.p.value, p.signif
  )

colnames(p3_complete_table) = c("Genotype", "lsmean", "SE", "Estimate comp.", "SE comp.", "adj. p-value", "signif.")



(p3_ldt_stats <- ggtexttable(p3_complete_table, theme = ttheme("blank", 
                                                               padding = unit(c(4, 4), "mm"),
                                                               base_size = 7.5), 
                             rows = NULL) %>% 
    table_cell_font( row=2:7, column = 1, face = "italic", size=7.5)%>% 
    tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 2, linetype = 1)%>% 
    tab_add_hline(at.row = c(10), row.side = "bottom", linewidth = 2, linetype = 1)
)


#### S. lycopersicoides #################


(p4 <- data_slopes %>% 
   filter(comm == "ok", 
          species=="S. lycopersicoides") %>%
   group_by(genotype) %>% 
   mutate(mean_ldt=mean(LDT, na.rm = T),
          n=n()) %>% 
   ungroup() %>% 
   left_join(pop_mean, by="genotype") %>% 
   summarise(species=species,
             genotype=genotype,
             LDT=LDT,
             SampleID=SampleID,
             start_date=start_date,
             mean_ldt = lsmean,
             n=n) %>% 
   unique() %>% 
   ggplot()+
   geom_boxplot(aes(x=reorder(factor(genotype), -mean_ldt), y=LDT, fill=genotype),
                show.legend = F, size=1, 
                fill=cols_lyco)+
   geom_text(aes(x=reorder(factor(genotype), -mean_ldt), y=4.5, label=n))+
   xlab("")+
   ylab("")+
   scale_y_continuous(limit=c(4.4,10.7), 
                      breaks=c(5.192957, 6.291569,6.984716,7.677864, 8.371011, 9.064158,10.6),
                      labels = c(3,9, 18, 36, 72, 144, 720))+
   theme_bw()+
   theme(axis.text = element_text(size=11, color="Black"),
         axis.text.x = element_text(size=11, angle=45,
                                    hjust=1, vjust=1, color="black"),
         axis.title = element_text(size=11, color="Black"),
         legend.title = element_text(size=13, 
                                     margin = margin(t = 0, r = 10, b = 0, l = 0)),
         legend.text = element_text(size=11), 
         legend.spacing.x = unit(.5, 'cm'),
         panel.grid.major.x =  element_blank()
   )+
   guides(fill=guide_legend(title="Plant Species",
                            label.theme = element_text(angle = 0, face = "italic")))

)

p4_complete_table <- p3_data %>% 
  separate(col=contrast, into = c("group1", "group2"), sep = ", ", remove = F) %>% 
  filter(group1 %in% c("S. lycopersicoides")) %>% 
  mutate(group1=factor(group1),
         group2=factor(group2)) %>% 
  dplyr::select(group1, group2, estimate, std.error, adj.p.value)  %>% 
  mutate(p.signif = symnum(adj.p.value, corr = FALSE, na = FALSE, 
                           cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", "ns"))
  ) %>% 
  summarise(species=group1,
            genotype = group2,
            estimate=estimate,
            std.error =std.error,
            adj.p.value=adj.p.value,
            p.signif = p.signif) %>% 
  left_join(p3_means, by="genotype" ) %>% 
  arrange(species) %>% 
  group_by(species) %>% 
  arrange(lsmean, .by_group = T) %>%
  mutate( adj.p.value = format.pval(pv = adj.p.value,
                                    digits = 2, 
                                    eps = 0.001,
                                    nsmall = 3),
          estimate = round(estimate, digits=2),
          lsmean = round(lsmean, digits=2),
          SE = round(SE, digits=2),
          std.error = round(std.error, digits=2)
  ) %>% 
  ungroup() %>% 
  dplyr::select(
    genotype, lsmean, SE, estimate,std.error, adj.p.value, p.signif
  )

colnames(p4_complete_table) = c("Genotype", "lsmean", "SE", "Estimate comp.", "SE comp.", "adj. p-value", "signif.")



(p4_ldt_stats <- ggtexttable(p4_complete_table, theme = ttheme("blank", 
                                                               padding = unit(c(4, 4), "mm"),
                                                               base_size = 7.5), 
                             rows = NULL) %>% 
    table_cell_font( row=2:7, column = 1, face = "italic", size=7.5)%>% 
    tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 2, linetype = 1)%>% 
    tab_add_hline(at.row = c(8), row.side = "bottom", linewidth = 2, linetype = 1)
)


#### ARRANGE plots LDT ####

(p_tot_ldt <- ggarrange(p3, p4,
                        p3_ldt_stats, p4_ldt_stats,
                        nrow = 2, 
                        ncol=2, 
                        heights = c(4,3,4,3),
                        labels=c("E", "F", "G", "H"),
                        align = "v"))


# massive plot 

(p_tot <- ggarrange(p_tot_lag, p_tot_ldt, 
                    nrow = 2))
ggsave("C://Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/figures/fig_3.svg",
       p_tot,
       dpi=600, 
       width = 10, height=12)
