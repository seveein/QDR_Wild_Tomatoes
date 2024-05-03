#### Statistics_bin (1/100F)
# 2023_10_26


rm(list=ls())
setwd("C:/Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/")
#data_full_adj <- read.csv("C:/Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/2023_09_08AllSlopes_Ann_wthMock_coContaminants_threeReps_filtered_adjusted.csv")
data_full <- read.csv("2023_09_08AllSlopesCounts_Ann_wthMock_coContaminants_threeReps_filtered.csv", header=T)



library(gplots)
library(lsmeans)
library(nlme)
library(tidyverse)
library(ggpubr)
library(slider)

data_filter <- data_full %>% filter(comm=="ok", 
                                    species=="S. pennellii")

data_filter$genotype <- as.factor(data_filter$genotype)
data_filter$species <- as.factor(data_filter$species)

genotypes <- c("LA1282", "LA1941", "LA1809")

data_pen_subset_bin100 <- data_full %>% 
  dplyr::filter(comm=="ok", 
                genotype %in% genotypes,
                inoculum=="ss",
                comm=="ok",
  ) %>%
  group_by(SampleID) %>% 
  mutate(share=lesion_area/leaf_area, 
         sliding_share = slide_dbl(share, ~ mean(.x), 
                                   .before = 3, .after = 3, 
                                   .complete=T),
         mean_leaf_size = mean(leaf_area)) %>% 
  dplyr::filter(!is.na(sliding_share)) %>%  
  #left_join(data_pen_subset_tt100, by="SampleID") %>% 
  mutate(inTime = ifelse(sliding_share< .95, 0, 1)) %>% 
  dplyr::select(genotype, SampleID, start_date, mean_leaf_size, inTime) %>% 
  unique()


### GET AND SHAPE DATA

#data_filter = data_pen_subset_bin100
data_pen_subset_bin100$genotype <- as.factor(data_pen_subset_bin100$genotype)

### GET BIN VALUES
data_bin<- data_pen_subset_bin100
#matrix_bin <- as.matrix(data_pen_subset_bin100)

# Generalized linear model 

mod_bin <- glm(inTime ~  genotype + start_date, data=data_bin,
               family=binomial(link="logit"))
mod_bin2 <- update(mod_bin, . ~ 0 + genotype + start_date)

# pseudo R^2:
#install.packages("piecewiseSEM")
library(piecewiseSEM)
rsquared(mod_bin2)


# anova:
library(car)
Anova(mod_bin)


library(multcomp)

###
# pairwise comparison 
###


## NOTE 
# we cannot compare on whole data set since other leaves have different leaf sizes!

#CoMaPops <- rbind(
#  "LA1282 - LA1809"=c(0,0,0, 1/2,  0,0,0,0,0,0,0,0,0,0,0,-1/2,   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,numeric(11)),
#  "LA1282 - LA1941"=c(0,0,0, 1/2,  0,0,0,0,0,0,0,0,0,0,0,   0,   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,numeric(11)),
#  "LA1809 - LA1941"=c(0,0,0,   0,  0,0,0,0,0,0,0,0,0,0,0, 1/2,-1/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,numeric(11)),

CoMaPops <- rbind(
  "genotypeLA1282 - genotypeLA1809"=c(1/2,-1/2,    0,numeric(2)),
  "genotypeLA1282 - genotypeLA1941"=c(1/2,   0, -1/2,numeric(2)),
  "genotypeLA1809 - genotypeLA1941"=c(  0,1/2, -1/2,numeric(2))
)
colnames(CoMaPops) <- names(coef(mod_bin))

compPops <- glht(mod_bin2, linfct=CoMaPops)
summary(compPops)

output_cld <- compPops %>% 
  broom::tidy() %>% 
  mutate(contrast = str_replace_all(contrast, "\\s", "")) %>% 
  pull(adj.p.value, contrast)  %>% 
  multcompView::multcompLetters()

write.csv(output_cld$Letters, "C:\\Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/stats_tables/2023_10_30_spen_reach100_pairwiseCLD.csv")


## GET MEAN ESTIMATES 

bin_dataframe_lsmeans <- as.data.frame(lsmeans(mod_bin , specs="genotype")) %>% 
  mutate(estimate = exp(lsmean)/ (1 + exp(lsmean))) %>% 
  dplyr::select(!df)


##
# convert logits:

logit <- lag_dataframe_lsmeans$lsmean
names(logit) <- lag_dataframe_lsmeans$genotype
exp(logit) / (1 + exp(logit))

###

#### Export stats to table


write.csv(bin_dataframe_lsmeans, "C:\\Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/stats_tables/2023_10_27_estimates_reach100_bin.csv")