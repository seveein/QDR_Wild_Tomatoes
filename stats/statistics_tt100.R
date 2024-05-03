### STATSTICS tt100 duration 
### 2023_10_30
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

genotypes <- c("LA1282", "LA1941", "LA1809")

data_filter <- data_full %>%
  dplyr::filter(inoculum=="ss",
                comm=="ok",
                genotype %in% genotypes) %>%
  group_by(SampleID) %>% 
  mutate(share=lesion_area/leaf_area, 
         sliding_share = slide_dbl(share, ~ mean(.x), 
                                   .before = 3, .after = 3, 
                                   .complete=T)) %>% 
  dplyr::filter(sliding_share > 0.95) %>% 
  slice(1)%>% 
  mutate(tt100 = time) %>% 
  dplyr::select(genotype, start_date, SampleID, tt100)

data_filter$genotype <- as.factor(data_filter$genotype)

# boxplot:
#windows(14,7); par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
boxplot(tt100 ~ genotype, data=data_filter, main="messgroesse", xlab="", ylab="", las=3, 
        col=c(rainbow(1),rainbow(9),rainbow(8),rainbow(8),rainbow(7)))



mod_tt100 <- gls(tt100 ~  genotype + start_date, data=data_filter, 
                 weights = varIdent(form=~1|genotype),
                 control=list(maxIter =500, msMaxIter= 500, 
                              niterEM=500, msMaxEval=500, opt="nlminb"))

mod_tt1002 <- update(mod_tt100, . ~ 0 + genotype + start_date)
# Sets first alphanumeric variant zero --> all following estimates are the actual mean values of each variant. 

windows(); plot(mod_tt1002)

mod_tt1002
tt100_dataframe_lsmeans <- as.data.frame(lsmeans(mod_tt1002, specs="genotype"))


write.csv(tt100_dataframe_lsmeans, "C:/Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/2023_10_27_tt100_Pops_Means.csv")

# pseudo R^2:
#install.packages("piecewiseSEM")
library(piecewiseSEM)
rsquared(mod_tt100)


# anova:
anova(mod_tt100)


library(multcomp)

###########################
# 1.: pairwise comparison #
###########################

coef(mod_tt1002)


CoMaPops <- rbind(
  "genotypeLA1282 - genotypeLA1809"=c(1/2,-1/2,    0,numeric(2)),
  "genotypeLA1282 - genotypeLA1941"=c(1/2,   0, -1/2,numeric(2)),
  "genotypeLA1809 - genotypeLA1941"=c(  0,1/2, -1/2,numeric(2))
)

colnames(CoMaPops) <- names(coef(mod_tt1002))

compPops <- glht(mod_tt1002, linfct=CoMaPops, df=mod_tt100$dims$N - mod_tt100$dims$p)
summary(compPops)

output_cld <- compPops %>% 
  broom::tidy() %>% 
  mutate(contrast = str_replace_all(contrast, "\\s", "")) %>% 
  pull(adj.p.value, contrast)  %>% 
  multcompView::multcompLetters()


#### OUTPUT stats to table 

table_compPops <- broom::tidy(compPops)

write.table(table_compPops, "C:/Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/2023_10_27_tt100_CompPops.csv")

# write cld
write.table(data.frame(output_cld$Letters), "stats_tables/2023_10_30_spen_tt100_pairwise_signCLD.txt")
