### STATSTICS AUDPC  
### 2023_10_30
rm(list=ls())

setwd("C:/Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/")
#data_full_adj <- read.csv("C:/Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/2023_09_08AllSlopes_Ann_wthMock_coContaminants_threeReps_filtered_adjusted.csv")
data_full <- read.csv("2023_09_08AllSlopesCounts_Ann_wthMock_coContaminants_threeReps_filtered.csv", header=T)



library(gplots)
library(lsmeans)
library(nlme)
library(tidyverse)

# get genotypes for analysis
genotypes <- c("LA1282", "LA1941", "LA1809")
# get data: derive audpc from leaf and lesion size, get unique values
data_filter <- data_full %>%
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

# set factors
data_filter$genotype <- as.factor(data_filter$genotype)
data_filter$species <- as.factor(data_filter$species)
###################################
######## Statistical model######### 
# generalized least squares model #
###################################


mod_audpc <- gls(audpc ~  genotype + start_date, data=data_filter, 
                 weights = varIdent(form=~1|genotype),
                 control=list(maxIter =500, msMaxIter= 500, 
                              niterEM=500, msMaxEval=500, opt="nlminb"))

mod_audpc2 <- update(mod_audpc, . ~ 0 + genotype + start_date)
# Sets first alphanumeric variant zero --> all following estimates are the actual mean values of each variant. 

windows(); plot(mod_audpc2)

mod_audpc2

######################
# Get mean estimates #
######################

audpc_dataframe_lsmeans <- as.data.frame(lsmeans(mod_audpc2, specs="genotype"))
write.csv(audpc_dataframe_lsmeans, "stats_tables/2023_10_30_spen_audpc_lsmeans.txt")

# pseudo R^2:
#install.packages("piecewiseSEM")
library(piecewiseSEM)
rsquared(mod_audpc2)


# anova:
anova(mod_audpc2)


library(multcomp)
###########################
# 1.: pairwise comparison #
###########################

coef(mod_audpc2)


CoMaPops <- rbind(
  "genotypeLA1282 - genotypeLA1809"=c(1/2,-1/2,    0,numeric(2)),
  "genotypeLA1282 - genotypeLA1941"=c(1/2,   0, -1/2,numeric(2)),
  "genotypeLA1809 - genotypeLA1941"=c(  0,1/2, -1/2,numeric(2))
)

colnames(CoMaPops) <- names(coef(mod_audpc2))

compPops <- glht(mod_audpc2, linfct=CoMaPops, df=mod_audpc$dims$N - mod_audpc$dims$p)
summary(compPops)

output_sign <- compPops %>% 
  broom::tidy() %>% 
  mutate(contrast = str_replace_all(contrast, "\\s", "")) %>% 
  pull(adj.p.value, contrast)  %>% 
  multcompView::multcompLetters()
write.table(data.frame(output_sign$Letters), "stats_tables/2023_10_30_spen_audpc_pairwise_signCLD.txt")