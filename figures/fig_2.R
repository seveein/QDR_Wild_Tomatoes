# Figure 2 of Solanum x Sclerotinia phenotyping paper
# Overview on QDR-parameters, Species LAG+LDT, IF
# as well as statistical tables 

rm(list=ls())
options(scipen=999)
library(tidyverse)
library(ggsignif)
library(ggpubr)
library(ochRe)
library(png)

setwd("C://Users/Sever/Nextcloud/ResiDEvo/")
setwd("C://Users/suaph281/Nextcloud/ResiDEvo/")
data_full <- read.csv("data/2023_ScleroPhenotypes/AllSlopes_Ann_wthMock_noContaminants_threeReps_filtered.csv")

########################
# P1 - Overview figure #
########################

img1 <- readPNG("figures/Disease_Parameters.png")


(p1 <- ggplot()+
    background_image(img1)+
    theme(plot.margin = margin(t=.0, l=1, r=1, b=.0, unit = "cm"),
          panel.background =  element_rect(fill = "transparent"))
)


###########
# P2 - IF #
###########

img2 <- readPNG("C:/Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/tables/tab1/tab_1.png")
(p2 <- ggplot()+
    background_image(img2)+
    theme(plot.margin = margin(t= 10, l=0, r=0, b=10, unit = "mm"),
          panel.background =  element_rect(fill = "transparent")))


#################
#   P3          #
# lag PHASE #
#################

(p3 <- data_full %>% 
    filter(comm == "ok"&
             genotype != "C32") %>%
    drop_na(lag) %>%
    group_by(species) %>% 
    mutate(lag_days = lag/(60*24)) %>% 
    mutate(mean_ldt=mean(LDT, na.rm = T),
           mean_lag=mean(lag_days, na.rm = T),
           n=n()) %>% 
    ungroup() %>% 
    summarise(species=species,
              genotype=genotype,
              lag=lag_days,
              LDT=LDT,
              SampleID=SampleID,
              start_date=start_date,
              mean_lag=mean_lag,
              mean_ldt = mean_ldt,
              n=n) %>% 
    unique() %>% 
    ggplot()+
    geom_boxplot(aes(x=reorder(factor(species), mean_lag), y=lag, fill=species), 
                 show.legend = F, size=1, col="black")+
    geom_text(aes(x=reorder(factor(species), mean_lag), y=-0.1, label=n))+
    xlab("")+
    ylab("Lag phase \n [days]")+
    theme_bw()+
    #ylim(c(-0.1, 10))+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_blank(),#element_text(size=11, angle=45,
          #hjust=1, vjust=1, color="black",
          #face="italic"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm')
    )+
    scale_y_continuous(limits=c(-.10,7), breaks = seq(0,8, by=2))+
    scale_fill_ochre(palette="parliament") +
    
    guides(fill=guide_legend(title="Plant Species",
                             label.theme = element_text(angle = 0, face = "italic")))
  
)
#############
#     P4    #
# lag-Stats #
#############

p4_data <- read.table("C:/Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/tables/tab2/2023_10_09_lag_CompSpecies.csv")
p4_data$contrast <- as.factor(p4_data$contrast)
p4_data_edited <- p4_data %>% 
  arrange(desc(contrast)) %>% 
  mutate(
    estimate = estimate/(60*24),
    std.error = std.error/(60*24),
           signif. = symnum(adj.p.value, corr = FALSE, na = FALSE, 
                            cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                            symbols = c("***", "**", "*", ".", "ns")),
         adj.p.value = format.pval(pv = adj.p.value,
                                   # digits : number of digits, but after the 0.0
                                   digits = 2, 
                                   # eps = the threshold value above wich the 
                                   # function will replace the pvalue by "<0.0xxx"
                                   eps = 0.001,
                                   # nsmall = how much tails 0 to keep if digits of 
                                   # original value < to digits defined
                                   nsmall = 3),
         estimate = round(estimate, digits=2),
         std.error = round(std.error, digits=2),
         statistic = round(statistic, digits=2)) %>% 
  select(!null.value & !adj.p.value & !statistic)
  



colnames(p4_data_edited) = c("Contrast", "Estimate", "Std. Error", "Signif.")

(p4 <- ggtexttable(p4_data_edited, theme = ttheme("blank", 
                                                  padding = unit(c(4, 4), "mm"),
                                                  base_size = 7), 
                   rows = NULL) %>% 
    table_cell_font( row=2:7, column = 1, face = "italic", size=7)%>% 
    tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 2, linetype = 1)%>% 
    tab_add_hline(at.row = c(7), row.side = "bottom", linewidth = 2, linetype = 1)
)


#######
# P5  #
# LDT #
#######

order <- data_full%>% 
  filter(comm == "ok" & 
           genotype != "C32") %>%
  group_by(species) %>% 
  mutate(mean_lag=mean(lag, na.rm = T)) %>%
  ungroup() %>% 
  dplyr::select(c("species", "mean_lag")) %>% 
  unique()
order_order <- levels(reorder(order$species, order$mean_lag))



(p5 <- data_full %>% 
    filter(comm == "ok" & 
             genotype != "C32") %>%
    drop_na(LDT) %>%
    group_by(species) %>% 
    mutate(mean_ldt=mean(LDT, na.rm = T),
           mean_lag=mean(lag, na.rm = T),
           n=n()) %>% 
    ungroup() %>% 
    summarise(species=species,
              genotype=genotype,
              lag=lag,
              LDT=(LDT),
              SampleID=SampleID,
              start_date=start_date,
              mean_lag=mean_lag,
              mean_ldt = mean_ldt,
              n=n) %>% 
    unique() %>% 
    #filter(LDT < 10000) %>% 
    ggplot()+
    geom_boxplot(aes(x=factor(species, levels=order_order), y=LDT, fill=species),
                 show.legend = F, size=.71, col="black")+
    geom_text(aes(x=factor(species, levels=order_order),
                  y=4, label=n))+
    #geom_jitter(alpha=.6, height = 0)+
    #geom_point(aes(y=mean_ldt), size=3, col="Red", show.legend = F)+
    ylab("Lesion doubling time \n LDT [h]")+
    #ylim(c(4,13))+
    xlab("")+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(size=11, angle=45,
                                     hjust=1, vjust=1, color="black",
                                     face="italic"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm'),
          panel.grid.minor.y = element_blank()
    )+ 
    scale_fill_ochre(palette="parliament") +
    scale_y_continuous(limit=c(4,10.7), breaks=c(4.094345, 5.192957, 6.291569,6.984716,7.677864, 8.371011, 9.064158,10.6),labels = c(1, 3,9, 18, 36, 72, 144, 720))
)

#############
#    p6     #
# LDT stats #
#############

p6_data <- read.table("C:/Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/tables/tab3/2024_01_19_ldt_CompSpecies.csv")
p6_data$contrast <- as.factor(p6_data$contrast)
p6_data_edited <- p6_data %>% 
  arrange(desc(contrast)) %>% 
  mutate(
    #estimate=exp(estimate), 
    #std.error=exp(std.error),
    signif. = symnum(adj.p.value, corr = FALSE, na = FALSE, 
                          cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                          symbols = c("***", "**", "*", ".", "ns")),
         adj.p.value = format.pval(pv = adj.p.value,
                                   # digits : number of digits, but after the 0.0
                                   digits = 2, 
                                   # eps = the threshold value above wich the 
                                   # function will replace the pvalue by "<0.0xxx"
                                   eps = 0.001,
                                   # nsmall = how much tails 0 to keep if digits of 
                                   # original value < to digits defined
                                   nsmall = 3),
         estimate = round(estimate, digits=2),
         std.error = round(std.error, digits=2),
         statistic = round(statistic, digits=2)) %>% 
  select(!null.value & !adj.p.value & !statistic) 
colnames(p6_data_edited) = c("Contrast", "Estimate", "Std. Error", "Signif.")

(p6 <- ggtexttable(p6_data_edited, theme = ttheme("blank", 
                                                  padding = unit(c(4, 4), "mm"),
                                                  base_size = 7), 
                   rows = NULL) %>% 
    table_cell_font( row=2:7, column = 1, face = "italic", size=7) %>% 
    tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 2, linetype = 1)%>% 
    tab_add_hline(at.row = c(7), row.side = "bottom", linewidth = 2, linetype = 1)
)



#########
# P_TOT #
######### 


(p_tot <- ggpubr::ggarrange(p1, p2, p3, p4,p5, p6,
                            labels=c("A", "B", "C", "D", "E", "F"), 
                            heights = c(2,2,3),
                            widths = c(4, 4.6),
                            #align = "v",
                            #hjust=-0.5,
                            #vjust=1.,
                            nrow=3, 
                            ncol=2,
                            common.legend = T
                            #legend.args = list(hjust = -1, vjust = 0)
))


ggsave("C://Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/figures/fig_2.png",
       p_tot, dpi=600, units = "in",
       height=7, width=7, bg = "white")


