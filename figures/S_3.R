# suppl. Fig for Infection frequency on pop-level
# Severin Einspanier
# 2024_04_25

rm(list=ls())
library(tidyverse)
library(ggrepel)
data_pop_IF <- read.csv("C:/Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/2024_04_15_IF_Populations_lsMeans.csv", header=T
                      )

# set positions 

pos <- position_jitter(width = 0.2, seed = 2)

(p1 <- data_pop_IF %>% 
    filter(genotype != "C32") %>% #mutate(species = sub(",.*", "", data_pop_IF$contrast),
    #       genotype = sub(".*?,", "", data_pop_IF$contrast)) %>% 
    ggplot(aes(x=species, y= p.IF., col=species, label=as.character(genotype)))+
    geom_point(position=pos, show.legend = F, size=2)+
    theme_bw()+
    xlab("")+ 
    ylab("p(IF)")+ 
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(size=11, angle=45,
                                     hjust=1, vjust=1, color="black", 
                                     face="italic"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm'),
          panel.grid.major.x =  element_blank()
          )+
    geom_text_repel(position=pos, 
                    color="grey50",  
                    size=3)
    )

ggsave("C:/Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/figures/S_3.png", 
       p1, height=5, width=5, dpi=600)

# % SE of logit 

(data_pop_IF %>% mutate(share = round((SE/logit.IF.), digits=2))
)
