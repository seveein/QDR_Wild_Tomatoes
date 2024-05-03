# get the strongest influence on lesion area (AUDPC)
# inspired from Caseys et al 2021

library(gplots)
library(lsmeans)
library(nlme)
library(tidyverse)
setwd("C://Users/suaph281/Nextcloud/ResiDEvo/")
data <- read.csv("data/2023_ScleroPhenotypes/2023_09_08AllSlopesCounts_Ann_wthMock_coContaminants_threeReps_filtered.csv",
                 header=T)


genotypes <- c("LA1282", "LA1941", "LA1809")

data_filter_audpc <- data %>%
  filter(inoculum=="ss",
         comm=="ok",
         time < 890 & time > 30, 
         genotype %in% genotypes) %>%
  group_by(SampleID) %>%
  mutate(share=lesion_area/leaf_area) %>%
  summarize(audpc = agricolae::audpc(share, time, type="absolute"), 
            species = species,
            genotype=genotype,
            SampleID=SampleID,
            LDT=LDT,
            lag=lag,
            start_date = start_date) %>% 
  unique() 

write.csv(data_filter_audpc, "C:/Users/suaph281/Desktop/local_data/2024_01_29_contribution_variance_mario.csv")

plot(data_filter_audpc$audpc)

data_filter_audpc$genotype <- as.factor(data_filter_audpc$genotype)
data_filter_audpc$species <- as.factor(data_filter_audpc$species)


mod_contr <- gls(audpc ~  genotype + lag + LDT + start_date+genotype*start_date, data=data_filter_audpc,
                      weights = varIdent(form=~1|genotype),
                      control=list(maxIter =500, msMaxIter= 500, 
                                   niterEM=500, msMaxEval=500, opt="nlminb"))
mod_contr2 <- update(mod_contr, . ~ 0 + genotype + lag + LDT + start_date)
mod_contr2
anova(mod_contr)
      