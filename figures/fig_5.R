# Figure S. pennellii Leaf size, AUDPC, tt100, 1/tt100 
# SEVERITY

rm(list=ls())
library(tidyverse)
library(agricolae)
library(ggpubr)
library(slider)
library(ochRe)

setwd("C://Users/suaph281/Nextcloud/ResiDEvo/")
data <- read.csv("data/2023_ScleroPhenotypes/AllCountsSlopes_Ann_wthMock_noContaminants_threeReps_filtered.csv",
                 header=T)

stats_leaf_size <- read.table("data/2023_ScleroPhenotypes/stats_tables/2023_10_26_spen_leafsize_pairwise.txt")
groups <- as.data.frame(stats_leaf_size) %>% 
  mutate(genotype = trt, 
         groups=label) %>% 
  dplyr::select(genotype, groups)

library(ggforce)

(p1 <- data %>%
    filter(species == "S. pennellii") %>% 
    group_by(SampleID) %>%
    mutate(mean_leaf =  mean(leaf_area)) %>% 
    ungroup() %>% 
    group_by(genotype) %>%
    dplyr::select(genotype, start_date, SampleID, mean_leaf) %>% 
    unique() %>%
    ggplot(aes(x=fct_reorder(genotype, -mean_leaf), y=mean_leaf, fill=genotype))+
    scale_fill_ochre(palette="parliament") +
    geom_rect(xmin=.5, xmax=3.5, ymin=500, 
              ymax=12000, alpha=.3, fill=NA, col="lightgrey",
              linejoin = "round",
              size=1.5, linetype="dashed")+
    #geom_mark_rect(xmin="LA1941", xmax="LA1303", ymin=0, 
    #               ymax=12500, alpha=.2)+
    geom_boxplot(show.legend = F)+
    #geom_point(aes(y=group_mean_leaf))+
    theme_bw()+
    xlab("")+
    #ylim(c(0,12000))+
    scale_y_continuous(limits = c(0, 13000), breaks = seq(0, 12000, by = 3000))+
    ylab("mean leaf area [px]")+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(size=11, angle=45,
                                     hjust=1, vjust=1, color="black",
                                     #face="italic"
          ),
          axis.title = element_text(size=11, color="Black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm')
    )+
    #scale_fill_discrete(guide = FALSE) +
    #scale_x_continuous(expand = c(0, 0.5),
    #                   labels = levels(data_pen_leaf_filtered_EXcluded$genotype),
    #                   breaks = 1:length(levels(data_pen_leaf_filtered_EXcluded$genotype))) +
    
    # add grouping labels 
    geom_text(data=groups, aes(x=genotype, label=groups), y=13000)
)


# Plot AUDPC


stats_audpc <- read.table("data/2023_ScleroPhenotypes/stats_tables/2023_10_30_spen_audpc_pairwise_signCLD.txt")
stats_audpc_adj <- stats_audpc %>% 
  rownames_to_column("gentoype") %>% 
  mutate(groups = output_sign.Letters)
stats_audpc_adj$genotype <- gsub('genotype','',stats_audpc_adj$gentoype)

genotypes <- c("LA1282", "LA1941", "LA1809")

data_filter_audpc <- data %>%
  filter(inoculum=="ss",
         comm=="ok",
         time < 890 & time > 30, 
         genotype %in% genotypes) %>%
  group_by(SampleID) %>%
  mutate(share=lesion_area/leaf_area) %>%
  #filter(leaf_area < p97.5 & leaf_area >p2.5) %>%
  summarize(audpc = agricolae::audpc(share, time, type="absolute"), 
            species = species,
            genotype=genotype,
            SampleID=SampleID,
            start_date = start_date) %>% 
  unique() 

(p_audpc <- data_filter_audpc %>% 
    ggplot(aes(x=genotype, y=audpc, fill=genotype))+
    geom_boxplot(show.legend = F)+ 
    theme_bw()+
    xlab("")+
    ylab("AUDPC")+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(size=11,  color="black",
                                     #face="italic"
          ),
          axis.title = element_text(size=11, color="Black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm')
    )+
    scale_y_continuous(limits = c(0, 600), breaks = seq(0, 750, by = 200))+
    scale_fill_ochre(palette="parliament")+
    geom_text(data=stats_audpc_adj, aes(x=genotype, label=groups), y=600)
  )


##############
# Plot tt100 #
##############
stats_tt100 <- read.table("data/2023_ScleroPhenotypes/stats_tables/2023_10_30_spen_tt100_pairwise_signCLD.txt")
stats_tt100_adj <- stats_tt100 %>% 
  rownames_to_column("gentoype") %>% 
  mutate(groups = output_cld.Letters)
stats_audpc_adj$genotype <- gsub('genotype','',stats_audpc_adj$gentoype)

data_pen_subset_tt100 <- data %>% 
  filter(inoculum=="ss",
         comm=="ok",
         genotype %in% genotypes) %>%
  group_by(SampleID) %>% 
  mutate(share=lesion_area/leaf_area, 
         sliding_share = slide_dbl(share, ~ mean(.x), 
                                   .before = 3, .after = 3, 
                                   .complete=T)) %>% 
  filter(sliding_share > 0.95) %>% 
  slice(1)%>% 
  mutate(tt100 = time) %>% 
  dplyr::select(genotype, species, SampleID, tt100 ) 

(p3 <- data_pen_subset_tt100 %>%
    mutate(tt100=(tt100*10/60)/24) %>% 
    ggplot(aes(x=genotype, y=tt100, fill=genotype))+
    geom_boxplot(show.legend = F)+
    theme_bw()+
    xlab("")+
    ylab("Time till 100 % severity \n[days]")+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(size=11, color="black",
                                     #face="italic"
          ),
          axis.title = element_text(size=11, color="Black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm')
    )+
    geom_text(data=stats_audpc_adj, aes(x=genotype, label=groups), y=7.6)+
    scale_fill_ochre(palette="parliament")+
    scale_y_continuous(limits = c(0, 8), breaks = seq(0, 10, by = 2.5))
)


# plot 4 How many samples are not infected fully at all?
# 
# p4 ggtexttable

tab_100f <- read.csv("C:\\Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/stats_tables/2023_10_27_estimates_reach100_bin.csv")
cld_100f <- read.csv("C:\\Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/stats_tables/2023_10_30_spen_reach100_pairwiseCLD.csv")
cld_100f$genotype <- gsub('genotype','',cld_100f$X)
cld_100f$groups <- cld_100f$x
tab_100f_complete <- tab_100f %>% 
  left_join(cld_100f, by="genotype") %>% 
  dplyr::select(genotype, lsmean, SE, estimate, groups)

library(ggpubr)
(p_100f <- ggpubr::ggtexttable(tab_100f_complete, theme = ttheme("blank"), rows = NULL) )

(p_100_ff <- p_100f %>% 
    tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
    tab_add_hline(at.row = c(4), row.side = "bottom", linewidth = 3, linetype = 1) %>%
    tab_add_vline(at.column = 2:tab_ncol(p_100f), column.side = "left", from.row = 2, linetype = 2) %>% 
  tab_add_title("100 % severity frequency", face = "bold"))


### Infection frequency for subset

tab_if <- read.csv("C:\\Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/stats_tables/2023_10_30_spen_subset_IF_estimates.csv")
cld_if <- read.csv("C:\\Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/stats_tables/2023_10_30_spen_subset_IF_pairwiseCLD.csv")
cld_if$genotype <- gsub('genotype','',cld_100f$X)
cld_if$groups <- cld_if$x
tab_if_complete <- tab_if %>% 
  left_join(cld_if, by="genotype") %>% 
  dplyr::select(genotype, lsmean, SE, estimate, groups)


(p_if <- ggpubr::ggtexttable(tab_if_complete, theme = ttheme("blank"), rows = NULL) )

(p_iff <- p_if %>% 
    tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
    tab_add_hline(at.row = c(4), row.side = "bottom", linewidth = 3, linetype = 1) %>%
    tab_add_vline(at.column = 2:tab_ncol(p_if), column.side = "left", from.row = 2, linetype = 2)%>% 
  tab_add_title("Infection Frequency", face = "bold")
)



# Plot 5

(p5  <- data_filter_audpc %>%
    left_join(data_pen_subset_tt100, by="SampleID") %>%
    mutate(tt100=(tt100*10/60)/24) %>% 
       ggplot(aes(x=tt100, y=audpc, add = "reg.line")) +
    xlab("tt100 [days]")+
    ylab("AUDPC")+
       stat_cor(
         aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
         label.x = 4,
         size=4
       )+
       geom_smooth(method="lm", col="#A86048", size=2, fill="#D8D8D8")+
    theme_bw()+
       geom_point(size=2, alpha=.5)+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(size=11, 
                                     hjust=1, vjust=1, color="black",
                                     #face="italic"
          ),
          axis.title = element_text(size=11, color="Black"),
          #panel.grid.major.x = element_blank(),
          #panel.grid.minor.y = element_blank(),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm')
    )+
    scale_color_ochre(palette="parliament")
     
)

library(png)

img1 <- readPNG("figures/2023_SolanumScleroAll/paper_figures/2023_10_09_panel.png")

(p_png <- ggplot()+
    background_image(img1)+
    theme(plot.margin = margin(t=.1, l=.3, r=.3, b=.1, unit = "cm")))
  
## big plot

(p_tot <- ggarrange(ggarrange(p1, p_png, 
                                nrow=1,
                                align = "h", widths = c(4,3,3), 
                                labels=c("A", "B")),
                    ggarrange(p_iff, p_100_ff, nrow=1,
                              align="hv"),
                    ggarrange(p3, p_audpc, p5, nrow=1, labels=c("D", "E","F"),
                              align = 'hv', widths = c(3,3,4)),
                    
                    nrow=3, heights = c(10,5,9), labels=c("", "C", ""))
)
  

ggsave("figures/2023_SolanumScleroAll/paper_figures/fig6_severity.png", 
       p_tot, dpi=600, height=8.5, width=10)

