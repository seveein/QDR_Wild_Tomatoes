## Cor. between IF and LDT
# suppl. fig. 5
# Severin Einspanier
# 2024_04_25

rm(list=ls())

library(tidyverse)
library(ochRe)
library(ggpubr)

# Correlation analysis based on the estimates from the models.

setwd("C:/Users/Sever/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/")
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




(p_3 <- ggplot(combined, aes(x=lsman_ldt, y=p.IF., col=species))+
  geom_point(size=0, alpha=0)+
  xlim(c(5.5,8))+
  geom_rect(xmin=7, xmax=8, ymin=1.05, ymax=1.50, fill="white", col="black")+
  geom_point(size=3, alpha=.7, show.legend = F)+
  geom_smooth(aes(fill=species), size=2, method="lm", show.legend = F,
              alpha = 0.1,
              se=F)+
  
  ggpubr::stat_cor(method="pearson", size=5, show.legend = F,
                   label.x = 7.1)+
  theme_bw()+
  ylab("Mean Estimates IF [%]")+
  xlab("Mean Estimates LDT")+
  scale_y_continuous(breaks=seq(0,1, by= .25), limits = c(0,1.5))+
  #scale_x_continuous(limits=c(0.8,3.6))+
  theme(axis.text = element_text(size=15, color="Black"),
        axis.text.x = element_text(size=15, angle=0,
                                   hjust=.5, vjust=1, color="black",
        ),
        axis.title = element_text(size=15, color="Black"),
        #legend.title = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size=15), 
        legend.spacing.x = unit(.5, 'cm'),
        legend.position = "top",
        legend.box.background = element_rect(colour = "black")
  )+
  #scale_line_continuous(guide=guide_legend(override.aes = list(shape=18)))+
  guides(color=guide_legend(title="Plant Species",
                            label.theme = element_text(angle = 0, face = "italic", 
                                                       size=13),
                            title.theme = element_text(size=15),
                            #keywidth = 2,
                            #keyheight = 2,
                            nrow=2,
                            override.aes = list(size=3, alpha=1)
  ))+
  scale_color_ochre(palette = "parliament")+
  scale_fill_ochre(palette = "parliament")
)

ggsave("C:/Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/figures/S_5.png",
       p_3, dpi=600, width=8, height=5)

