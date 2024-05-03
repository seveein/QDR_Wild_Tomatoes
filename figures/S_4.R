# Pooled correlation analysis of all sampels
# suppl. fig. 4
# Severin Einspanier
# 2024_04_25


rm(list=ls())

library(tidyverse)
library(ochRe)
library(ggpubr)

# Correlation analysis based on the estimates from the models.

#setwd("C:/Users/Sever/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/")
setwd("C:/Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/")
#meta_data <- read.csv("C:/Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/AllSlopes_Ann_wthMock_coContaminants_threeReps_filtered.csv")
meta_data <- read.csv("C:/Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/AllSlopes_Ann_wthMock_noContaminants_threeReps_filtered.csv")
meta_data <- meta_data %>% 
  dplyr::select(species, genotype) %>% 
  unique()

estimates_lag  <- read.csv("2023_10_09_lag_Pops_Means.csv")
estimates_if <- read.table(file = "2023_10_27_infection_frequency_Pops_estimates.csv")
estimates_ldt  <- read.csv("2023_10_24_ldt_Populations_Means.csv")

combined <- estimates_lag %>% 
  rename(lsman_lag = lsmean) %>% 
  left_join(estimates_if, by="genotype") %>%
  rename(lsman_if = lsmean) %>% 
  left_join(estimates_ldt, by="genotype") %>%
  rename(lsman_ldt = lsmean) %>% 
  left_join(meta_data, by="genotype") %>% 
  filter(genotype != "C32") %>% 
  dplyr::select(species, genotype, lsman_lag, lsman_ldt, p.IF.)



(p_1_pooled <- combined %>% 
  ggplot(aes(x=lsman_lag/(60*24), y=lsman_ldt))+
  geom_point(size=0, alpha=0, show.legend = F)+
  ylim(c(5,8.4))+
  ylab("adj. log(LDT)")+
  xlab("adj. Lag Phase Duration [days]")+
  geom_rect(xmin=1, xmax=2.1, ymin=8, ymax=8.2, fill="white", col="black")+
  geom_point(size=3, alpha=.7, show.legend = F)+
  geom_smooth(size=2, method="lm", show.legend = F, se=F)+
  ggpubr::stat_cor(method="pearson",  size=5, show.legend = F,
                   p.accuracy = 0.01, r.accuracy = 0.01,
                   label.x = 1.1, label.y = 8.05,alpha=1,fontface="bold")+
  #stat_cor(method = "pearson", cor.coef.name = "r", ) +
  theme_bw()+
  theme(axis.text = element_text(size=15, color="Black"),
        axis.text.x = element_text(size=15, angle=0,
                                   hjust=0.5, vjust=0, color="black",
        ),
        axis.title = element_text(size=15, color="Black"),
        legend.title = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size=15), 
        legend.spacing.x = unit(.5, 'cm')
  )+
  #scale_line_continuous(guide=guide_legend(override.aes = list(shape=18)))+
  guides(color=guide_legend(title="Plant Species",
                            label.theme = element_text(angle = 0, face = "italic"),
                            #keywidth = 2,
                            #keyheight = 2,
                            override.aes = list(size=10, alpha=1)
  ))+
  scale_color_ochre(palette = "parliament")
  #scale_fill_ochre(palette = "parliament")
)
(p_2_pooled <- ggplot(combined, aes(x=lsman_lag/(60*24), y=p.IF.))+
    geom_point(size=0, alpha=0)+
    geom_rect(xmin=2.05, xmax=3.4, ymin=1, ymax=1.10, fill="white", col="black")+
    geom_point(size=3, alpha=.7, show.legend = F)+
    geom_smooth(size=2, method="lm", show.legend = F,
                alpha = 0.1,
                se=F)+
    
    ggpubr::stat_cor(method="pearson", size=5, show.legend = F,
                     label.x = 2.1, label.y = 1.05)+
    theme_bw()+
    ylab("adj. Infection frequency [%]")+
    xlab("adj. Mean Lag Phase Duration [days]")+
    scale_y_continuous(breaks=seq(0,1, by= .25), limits = c(0,1.15))+
    scale_x_continuous(limits=c(0.8,3.6))+
    theme(axis.text = element_text(size=15, color="Black"),
          axis.text.x = element_text(size=15, angle=0,
                                     hjust=.5, vjust=1, color="black",
          ),
          axis.title = element_text(size=15, color="Black"),
          legend.title = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=15), 
          legend.spacing.x = unit(.5, 'cm')
    )+
    #scale_line_continuous(guide=guide_legend(override.aes = list(shape=18)))+
    guides(color=guide_legend(title="Plant Species",
                              label.theme = element_text(angle = 0, face = "italic"),
                              #keywidth = 2,
                              #keyheight = 2,
                              override.aes = list(size=10, alpha=1)
    ))+
    scale_color_ochre(palette = "parliament")+
    scale_fill_ochre(palette = "parliament")
)

(p_tot_pooled <- ggarrange(p_1_pooled, p_2_pooled, nrow=1))
ggsave("C://Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/figures/S_4.png",
       dpi=600, width=11, height=6)
