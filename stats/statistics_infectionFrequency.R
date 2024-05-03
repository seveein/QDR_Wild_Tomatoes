#### Statistics_bin
# 2023_10_05


rm(list=ls())

setwd("C:/Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/")
setwd("C:/Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/")
#data_full_adj <- read.csv("C:/Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/2023_09_08AllSlopes_Ann_wthMock_coContaminants_threeReps_filtered_adjusted.csv")
data_full <- read.csv("AllSlopes_Ann_wthMock_coContaminants_threeReps_filtered.csv", header=T)


### IMPORT 

#library(gplots)
library(lsmeans)
library(nlme)
library(tidyverse)
library(piecewiseSEM)

### GET AND SHAPE DATA

data_filter = data_full
data_filter$genotype <- as.factor(data_filter$genotype)
data_filter$species <- as.factor(data_filter$species)
#data_filter$spgen <- as.factor(paste(data_filter$species, data_filter$genotype))

### GET BIN VALUES


data_bin <- data_filter %>% 
  filter(inoculum =="ss") %>% 
  mutate(infection = ifelse(comm=="ok", 1, 0)) %>% 
  dplyr::select(SampleID,species,genotype,start_date,infection)%>% 
  na.omit()

data_bin$species <- as.factor(data_bin$species)
data_bin$genotype <- as.factor(data_bin$genotype)

matrix_bin <- as.matrix(data_bin)

# Generalized linear model 

mod_bin <- glm(infection ~  genotype + start_date, data=data_bin,
               family=binomial(link="logit"))
mod_bin2 <- update(mod_bin, . ~ 0 + genotype + start_date)

# pseudo R^2:
#install.packages("piecewiseSEM")
rsquared(mod_bin)


# anova:
library(car)
Anova(mod_bin)


library(multcomp)

# 1.: Which Species differ?
###########################

coef(mod_bin)




CoMaSpecies <- rbind(
  "S. habrochaites - S. lycopersicoides"=    c(0,   0,    0,    0,    0,    0,    0,    0,  1/8,    0,    0,    0,  1/8,  1/8,  1/8,    0,    0, -1/7,  1/8,  1/8,    0,  1/8,    0,    0, -1/7, -1/7, -1/7,    0,  1/8,-1/7,   0,    0,-1/7,-1/7,    0,numeric(11)),
  "S. lycopersicoides - S. pimpinellifolium"=c(0,   0,-1/10,    0,    0,-1/10,-1/10,-1/10,    0,-1/10,    0,-1/10,    0,    0,    0,    0,    0,  1/7,    0,    0,-1/10,    0,    0,    0,  1/7,  1/7,  1/7,-1/10,    0, 1/7,   0,-1/10, 1/7, 1/7,-1/10,numeric(11)),
  "S. pimpinellifolium - S. habrochaites"   =c(0,   0, 1/10,    0,    0, 1/10, 1/10, 1/10, -1/8, 1/10,    0, 1/10, -1/8, -1/8, -1/8,    0,    0,    0, -1/8, -1/8, 1/10, -1/8,    0,    0,    0,    0,    0, 1/10, -1/8,   0,   0, 1/10,   0,   0, 1/10,numeric(11)),
  "S. pennellii - S. habrochaites"         =c(0, 1/9,    0,  1/9,  1/9,    0,    0,    0, -1/8,    0,  1/9,    0, -1/8, -1/8, -1/8,  1/9,    0,    0, -1/8, -1/8,    0, -1/8,  1/9,  1/9,    0,    0,    0,    0, -1/8,   0, 1/9,    0,   0,   0,    0,numeric(11)),
  "S. pennellii - S. lycopersicoides"      =c(0, 1/9,    0,  1/9,  1/9,    0,    0,    0,    0,    0,  1/9,    0,    0,    0,    0,  1/9,    0, -1/7,    0,    0,    0,    0,  1/9,  1/9, -1/7, -1/7, -1/7,    0,    0,-1/7, 1/9,    0,-1/7,-1/7,    0,numeric(11)),
  "S. pennellii - S. pimpinellifolium"     =c(0, 1/9,-1/10,  1/9,  1/9,-1/10,-1/10,-1/10,    0,-1/10,  1/9,-1/10,    0,    0,    0,  1/9,    0,    0,    0,    0,-1/10,    0,  1/9,  1/9,    0,    0,    0,-1/10,    0,   0, 1/9,-1/10,   0,   0,-1/10,numeric(11))
)
colnames(CoMaSpecies) <- names(coef(mod_bin))

compSpecies <- glht(mod_bin2, linfct=CoMaSpecies)
#output <- summary(compSpecies)
#cld(output, level = 0.05)

#clds <- symnum(output$test$pvalues, corr = FALSE, na = FALSE, 
#                                    cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", "ns"))
#library(multcompView)
#ätable <- broom::tidy(compSpecies)
#multcomp::cld(compSpecies)

# 2.: Which Populations differ (within Species)?
################################################

coef(mod_bin)

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
  "S. habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0, 7/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S. habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8, 7/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S. habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8, 7/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S. habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0, 7/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S. habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8, 7/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S. habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0, 7/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S. habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0, 7/8,0,0,0,0,0,0,numeric(11)),
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

colnames(CoMaPops) <- names(coef(mod_bin))

compPops <- glht(mod_bin2, linfct=CoMaPops)
summary(compPops)

## GET MEAN ESTIMATES per species
CoMaSpMeans <- CoMaPops
CoMaSpMeans[CoMaSpMeans==9/10] <- -1/10
CoMaSpMeans[CoMaSpMeans==8/9]  <- -1/9
CoMaSpMeans[CoMaSpMeans==7/8]  <- -1/8
CoMaSpMeans[CoMaSpMeans==6/7]  <- -1/7
CoMaSpMeans[CoMaSpMeans<0] <- CoMaSpMeans[CoMaSpMeans<0]*(-1)
CoMaSpMeans <- unique(CoMaSpMeans)
compSpMeans <- glht(mod_bin2, linfct=CoMaSpMeans)
summary(compSpMeans)
logitSp <- summary(compSpMeans)$test$coefficients

(exp_logitSp <- exp(logitSp) / (1 + exp(logitSp)))

means_sp <- data.frame(logitSp, 
                       summary(compSpMeans)$test$sigma,
                       exp_logitSp)
colnames(means_sp) = c("logit(IF)", "SE", "p(IF)")
row.names(means_sp) = c("S. pennellii", "S. lycopersicoides", "S. habrochaites", "S. pimpinellifolium")

write.csv(means_sp, "C:\\Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/stats_tables/2023_10_16_IF_species_estimates.csv")

##
# convert logits:

logit <- lag_dataframe_lsmeans$lsmean
names(logit) <- lag_dataframe_lsmeans$genotype
exp(logit) / (1 + exp(logit))

###

logitSp <- summary(compPops)$test$coefficients

(exp_logitSp <- exp(logitSp) / (1 + exp(logitSp)))

table_compPops <- broom::tidy(compPops)
means_sp <- data.frame(table_compPops$contrast,
                        logitSp, 
                       summary(compPops)$test$sigma,
                       exp_logitSp)
colnames(means_sp) = c("contrast", "logit(IF)", "SE", "p(IF)")
#row.names(means_sp) = c("S. pennellii", "S. lycopersicoides", "S. habrochaites", "S. pimpinellifolium")


#####
# get estimates per accession
#####

IF_dataframe_lsmeans <- as.data.frame(lsmeans(mod_bin, specs="genotype"))



IF_dataframe_lsmeans_output <- data.frame("genotype"=IF_dataframe_lsmeans$genotype,
                                          "logit(IF)"= IF_dataframe_lsmeans$lsmean,
                                          "SE"=IF_dataframe_lsmeans$SE,
                                          "p(IF)"=exp(IF_dataframe_lsmeans$lsmean) / (1 + exp(IF_dataframe_lsmeans$lsmean))
                                          )
IF_dataframe_lsmeans_output_print <- data_bin %>% 
  dplyr::select(genotype, species) %>% 
  unique() %>% 
  mutate(genotype = as.factor(genotype)) %>% 
  left_join(IF_dataframe_lsmeans_output, by="genotype")

write.csv(IF_dataframe_lsmeans_output_print,"C:/Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/2024_04_15_IF_Populations_lsMeans.csv")



#### Export stats to table


#table_compPops <- broom::tidy(compPops)
table_comp_Species <- broom::tidy(compSpecies)

write.table(means_sp, "C:/Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/2024_01_09_infection_frequency_CompPops.csv")
write.table(table_comp_Species, "C:/Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/2023_10_09_infection_frequency_CompSpecies.csv")

#####


# The statistical software R (2023) was used to evaluate the data. The data evaluation startet with 
# the definition of an appropriate statistical 

# - generalized linear model (McCullagh and Nelder, 1989). 

# The model included ... (Level 1, Level 2, ... ), ... (...) and ... (...), as well as all 
# their interaction terms (two-fold and three-fold) as fixed factors. Also, ... and ... were modelled 
# as covarites. The block (and the plots, nested in block) was (were) regarded as random factor(s). 
# Also, the correlations of the measurement values due to the several levels of ... were taken into 
# account. The residuals were assumed to be normally distributed and to be homo/heteroscedastic with 
# respect to the different levels of .... These assumptions are based on a graphical residual 
# analysis. 
# Based on this model, a Pseudo R^2 was calculated (Nakagawa and Schielzeth, 2013) and an analysis of 
# covariances (ANCOVA) was conducted (Cochran, 1957), followed by multiple contrast tests (Hothorn et 
# al., 2008; see also Bretz et al., 2011; Schaarschmidt and Vaas, 2009) in order to 
# compare the several levels of the influence factors, respectively. If a factor of interest had no 
# significant interactions with the remaining factors, then the levels of these remainig factors were 
# pooled.
#
#
# R Core Team (2023). R: A language and environment for statistical computing. R Foundation for 
# Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, URL http://www.R-project.org/
#
# Laird, N. M. and Ware, J. H. (1982). Random-Effects Models For Longitudinal Data. Biometrics, 
# 38(4), 963–974.
#
# Pinheiro, J.C., and Bates, D.M. (2000) "Mixed-Effects Models in S and S-PLUS", Springer, New York.
#
# Box, G. E. P., Jenkins, G. M., and Reinsel G. C. (1994). Time Series Analysis: Forecasting and 
# Control, 3rd Edition, Holden-Day.
#
# Carroll, R. J. and Ruppert, D. (1988). Transformation and Weighting in Regression, Chapman and Hall. 
#
# McCullagh P. and Nelder, J. A. (1989). Generalized Linear Models, Chapman and Hall.
#
# Schall, R. (1991). Estimation in generalized linear models with random effects. Biometrika 78, 719–727. 
#
# Gelman, A., Jakulin, A., Pittau, M. G. and Su, Y.-S. (2009). A Weakly Informative Default Prior 
# Distribution For Logistic And Other Regression Models. The Annals of Applied Statistics 2 (4): 1360–1383. 
#
# Nakagawa, S. and Schielzeth, H. (2013). A General and Simple Method for Obtaining R2 from Generalized 
# Linear Mixed-Effects Models. Edited by Robert B. O’Hara. Methods in Ecology and Evolution 4, no. 2: 
# 133–42. 
#
# Cochran, W. G. (1957). Analysis Of Covariance - Its Nature And Uses, Biometrics, 13(3), 261-281
#
# Hothorn, T., Bretz, F., and Westfall, P. (2008). Simultaneous Inference in General Parametric Models, 
# Biometrical Journal, 50(3), 346-363.
#
# Bretz, F., Hothorn, T., and Westfall, P. (2011). Multiple Comparisons Using R, Chapman and Hall/CRC, 
# London, isbn 9781584885740.
#
# Schaarschmidt, F. and Vaas, L. (2009). Analysis of Trials with Complex Treatment Structure Using 
# Multiple Contrast Tests, Hortscience, 44(1), 188-195.
#


