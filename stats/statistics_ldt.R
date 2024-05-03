### STATSTICS LDT PHASE 
### 2023_10_05
## Stats with log(values)
## decided in agreement with mario on 2023_10_12

rm(list=ls())

setwd("C:/Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/")
#data_full_adj <- read.csv("C:/Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/2023_09_08AllSlopes_Ann_wthMock_coContaminants_threeReps_filtered_adjusted.csv")
data_full <- read.csv("AllSlopes_Ann_wthMock_noContaminants_threeReps_filtered.csv", header=T)


library(multcomp)
library(gplots)
library(lsmeans)
library(nlme)
library(tidyverse)

data_filter <- data_full %>% filter(comm=="ok")

data_filter$genotype <- as.factor(data_filter$genotype)
data_filter$species <- as.factor(data_filter$species)
data_filter$spgen <- as.factor(paste(data_filter$species, data_filter$genotype))

# boxplot:
windows(14,7); par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
boxplot(slope ~ genotype, data=data_filter, main="messgroesse", xlab="", ylab="", las=3, 
        col=c(rainbow(1),rainbow(9),rainbow(8),rainbow(8),rainbow(7)))

mod_ldt <- gls(LDT ~  genotype + start_date, data=data_filter, 
                      weights = varIdent(form=~1|genotype),
                      control=list(maxIter =500, msMaxIter= 500, 
                                   niterEM=500, msMaxEval=500, opt="nlminb"))

mod_ldt2 <- update(mod_ldt, . ~ 0 + genotype + start_date)
# Sets first alphanumeric variant zero --> all following estimates are the actual mean values of each variant. 

windows(); plot(mod_ldt)

mod_ldt2
as.data.frame(lsmeans(mod_ldt, specs="genotype"))


ldt_dataframe_lsmeans <- as.data.frame(lsmeans(mod_ldt, specs="genotype"))

write.csv(ldt_dataframe_lsmeans,"C:/Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/2024_01_19_ldt_Populations_Means.csv")


# pseudo R^2:
#install.packages("piecewiseSEM")
library(piecewiseSEM)
rsquared(mod_ldt2)


# anova:
anova(mod_ldt2)



#############################################
# 1.: Which Species differ from each other ?#
#############################################

coef(mod_ldt2)


CoMaSpecies <- rbind(
  "S. habrochaites - S. lycopersicoides"=    c(0,   0,    0,    0,    0,    0,    0,    0,  1/8,    0,    0,    0,  1/8,  1/8,  1/8,    0,    0, -1/7,  1/8,  1/8,    0,  1/8,    0,    0, -1/7, -1/7, -1/7,    0,  1/8,-1/7,   0,    0,-1/7,-1/7,    0,numeric(11)),
  "S. lycopersicoides - S. pimpinellifolium"=c(0,   0,-1/10,    0,    0,-1/10,-1/10,-1/10,    0,-1/10,    0,-1/10,    0,    0,    0,    0,    0,  1/7,    0,    0,-1/10,    0,    0,    0,  1/7,  1/7,  1/7,-1/10,    0, 1/7,   0,-1/10, 1/7, 1/7,-1/10,numeric(11)),
  "S. pimpinellifolium - S. habrochaites"   =c(0,   0, 1/10,    0,    0, 1/10, 1/10, 1/10, -1/8, 1/10,    0, 1/10, -1/8, -1/8, -1/8,    0,    0,    0, -1/8, -1/8, 1/10, -1/8,    0,    0,    0,    0,    0, 1/10, -1/8,   0,   0, 1/10,   0,   0, 1/10,numeric(11)),
  "S. pennellii - S. habrochaites"         =c(0, 1/9,    0,  1/9,  1/9,    0,    0,    0, -1/8,    0,  1/9,    0, -1/8, -1/8, -1/8,  1/9,    0,    0, -1/8, -1/8,    0, -1/8,  1/9,  1/9,    0,    0,    0,    0, -1/8,   0, 1/9,    0,   0,   0,    0,numeric(11)),
  "S. pennellii - S. lycopersicoides"      =c(0, 1/9,    0,  1/9,  1/9,    0,    0,    0,    0,    0,  1/9,    0,    0,    0,    0,  1/9,    0, -1/7,    0,    0,    0,    0,  1/9,  1/9, -1/7, -1/7, -1/7,    0,    0,-1/7, 1/9,    0,-1/7,-1/7,    0,numeric(11)),
  "S. pennellii - S. pimpinellifolium"     =c(0, 1/9,-1/10,  1/9,  1/9,-1/10,-1/10,-1/10,    0,-1/10,  1/9,-1/10,    0,    0,    0,  1/9,    0,    0,    0,    0,-1/10,    0,  1/9,  1/9,    0,    0,    0,-1/10,    0,   0, 1/9,-1/10,   0,   0,-1/10,numeric(11))
)
colnames(CoMaSpecies) <- names(coef(mod_ldt2))

compSpecies <- glht(mod_ldt2, linfct=CoMaSpecies, df=mod_ldt$dims$N - mod_ldt$dims$p)
compSpecies
summary(compSpecies)

##########################################################
# 2.: Which Populations differ (from species Grand Mean)?#
##########################################################

coef(mod_ldt2)

CoMaPops <- rbind(
  "S. pennellii, LA0716"=c(0, 8/9,0,-1/9,-1/9,0,0,0,0,0,-1/9,0,0,0,0,-1/9,-1/9,0,0,0,0,0,-1/9,-1/9,0,0,0,0,0,0,-1/9,0,0,0,0,numeric(11)),
  "S. pennellii, LA1282"=c(0,-1/9,0, 8/9,-1/9,0,0,0,0,0,-1/9,0,0,0,0,-1/9,-1/9,0,0,0,0,0,-1/9,-1/9,0,0,0,0,0,0,-1/9,0,0,0,0,numeric(11)),
  "S. pennellii, LA1303"=c(0,-1/9,0,-1/9, 8/9,0,0,0,0,0,-1/9,0,0,0,0,-1/9,-1/9,0,0,0,0,0,-1/9,-1/9,0,0,0,0,0,0,-1/9,0,0,0,0,numeric(11)),
  "S. pennellii, LA1656"=c(0,-1/9,0,-1/9,-1/9,0,0,0,0,0, 8/9,0,0,0,0,-1/9,-1/9,0,0,0,0,0,-1/9,-1/9,0,0,0,0,0,0,-1/9,0,0,0,0,numeric(11)),
  "S. pennellii, LA1809"=c(0,-1/9,0,-1/9,-1/9,0,0,0,0,0,-1/9,0,0,0,0, 8/9,-1/9,0,0,0,0,0,-1/9,-1/9,0,0,0,0,0,0,-1/9,0,0,0,0,numeric(11)),
  "S. pennellii, LA1941"=c(0,-1/9,0,-1/9,-1/9,0,0,0,0,0,-1/9,0,0,0,0,-1/9, 8/9,0,0,0,0,0,-1/9,-1/9,0,0,0,0,0,0,-1/9,0,0,0,0,numeric(11)),
  "S. pennellii, LA2657"=c(0,-1/9,0,-1/9,-1/9,0,0,0,0,0,-1/9,0,0,0,0,-1/9,-1/9,0,0,0,0,0, 8/9,-1/9,0,0,0,0,0,0,-1/9,0,0,0,0,numeric(11)),
  "S. pennellii, LA2719"=c(0,-1/9,0,-1/9,-1/9,0,0,0,0,0,-1/9,0,0,0,0,-1/9,-1/9,0,0,0,0,0,-1/9, 8/9,0,0,0,0,0,0,-1/9,0,0,0,0,numeric(11)),
  "S. pennellii, LA2963"=c(0,-1/9,0,-1/9,-1/9,0,0,0,0,0,-1/9,0,0,0,0,-1/9,-1/9,0,0,0,0,0,-1/9,-1/9,0,0,0,0,0,0, 8/9,0,0,0,0,numeric(11)),
  "S. lycopersicoides, LA1964"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 6/7,0,0,0,0,0,0,-1/7,-1/7,-1/7,0,0,-1/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S. lycopersicoides, LA2772"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0, 6/7,-1/7,-1/7,0,0,-1/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S. lycopersicoides, LA2776"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7, 6/7,-1/7,0,0,-1/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S. lycopersicoides, LA2777"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7,-1/7, 6/7,0,0,-1/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S. lycopersicoides, LA2951"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7,-1/7,-1/7,0,0, 6/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S. lycopersicoides, LA4123"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7,-1/7,-1/7,0,0,-1/7,0,0, 6/7,-1/7,0,numeric(11)),
  "S. lycopersicoides, LA4130"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7,-1/7,-1/7,0,0,-1/7,0,0,-1/7, 6/7,0,numeric(11)),
  "S. habrochaites, LA1559"= c(0,0,0,0,0,0,0,0, 7/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S. habrochaites, LA1721"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0, 7/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S. habrochaites, LA1731"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8, 7/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S. habrochaites, LA1753"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8, 7/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S. habrochaites, LA2128"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0, 7/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S. habrochaites, LA2167"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8, 7/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S. habrochaites, LA2409"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0, 7/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S. habrochaites, LA2864"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0, 7/8,0,0,0,0,0,0,numeric(11)),
  "S. pimpinellifolium, LA1261"=c(0,0, 9/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S. pimpinellifolium, LA1332"=c(0,0,-1/10,0,0, 9/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S. pimpinellifolium, LA1348"=c(0,0,-1/10,0,0,-1/10, 9/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S. pimpinellifolium, LA1374"=c(0,0,-1/10,0,0,-1/10,-1/10, 9/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S. pimpinellifolium, LA1593"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0, 9/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S. pimpinellifolium, LA1659"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0, 9/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S. pimpinellifolium, LA2347"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0, 9/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S. pimpinellifolium, LA2853"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0, 9/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S. pimpinellifolium, LA2983"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0, 9/10,0,0,-1/10,numeric(11)),
  "S. pimpinellifolium, LA4713"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0, 9/10,numeric(11))
  
)

colnames(CoMaPops) <- names(coef(mod_ldt2))

compPops <- glht(mod_ldt2, linfct=CoMaPops, df=mod_ldt$dims$N - mod_ldt$dims$p)
summary(compPops)

######################################
#### 3. Get adjusted SPECIES mean #### 
######################################
coef(mod_ldt2)

# WORKS

CoMaSpecies_2 <- rbind(
  "S. habrochaites"=    c(0,    0,    0,    0,    0,    0,    0,    0,  1/8,    0,    0,    0,  1/8,  1/8,  1/8,    0,   0,    0,  1/8,  1/8,    0,  1/8,    0,    0,    0,    0,    0,    0,  1/8,   0,   0,    0,   0,   0,    0,numeric(11)),
  "S. lycopersicoides"= c(0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,   0,  1/7,    0,    0,    0,    0,    0,    0,  1/7,  1/7,  1/7,    0,    0, 1/7,   0,    0, 1/7, 1/7,    0,numeric(11)),
  "S. pimpinellifolium"=c(0,    0, 1/10,    0,    0, 1/10, 1/10, 1/10,    0, 1/10,    0, 1/10,    0,    0,    0,    0,   0,    0,    0,    0, 1/10,    0,    0,    0,    0,    0,    0, 1/10,    0,   0,   0, 1/10,   0,   0, 1/10,numeric(11)),
  "S. pennellii"       =c(0,  1/9,     0, 1/9,  1/9,    0,    0,    0,    0,    0,  1/9,    0,    0,    0,    0,  1/9, 1/9,    0,    0,    0,    0,    0,  1/9,  1/9,    0,    0,    0,    0,    0,   0, 1/9,    0,   0,   0,    0,numeric(11))
)
colnames(CoMaSpecies_2) <- names(coef(mod_ldt2))

compSpecies_2 <- glht(mod_ldt2, linfct=CoMaSpecies_2, df=mod_ldt$dims$N - mod_ldt$dims$p)
compSpecies_2
summary(compSpecies_2)

data_id <- data_full %>%
  filter(genotype!="LA2727", 
         genotype!="C32") %>% 
  dplyr::select(species, genotype) %>% 
  unique() %>% 
  mutate(contrast = species) %>% 
  left_join(broom::tidy(compSpecies_2), by="contrast") %>% 
  dplyr::select(species, estimate, std.error) %>% 
  unique()

write.csv(data_id,"C:/Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/2024_04_24_ldt_Species_Means.csv")

#### OUTPUT stats to table 

table_compPops <- broom::tidy(compPops)
table_comp_Species <- broom::tidy(compSpecies)

table_compPops_stars <- table_compPops %>% 
  mutate(signif = symnum(adj.p.value, corr = FALSE, na = FALSE, 
                          cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                          symbols = c("***", "**", "*", ".", "ns")))

write.csv(table_compPops_stars, "C:/Users/suaph281/Nextcloud/AG_STAM/PapersProgress/Phenotyping_Solanum_Sclerotinia/tables/supp_table4/2024_04_25_LDT_GrandMean.csv")
write.table(table_comp_Species, "C:/Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/2024_01_19_ldt_CompSpecies.csv")
