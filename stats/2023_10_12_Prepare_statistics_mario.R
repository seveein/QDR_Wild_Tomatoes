# Script stat-session 13.10.2023
rm(list=ls())

library(gplots)
library(lsmeans)
library(nlme)
library(multcomp)
library(piecewiseSEM)

setwd("C://Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/")
data_full <- read.csv("2023_09_08AllSlopes_Ann_wthMock_coContaminants_threeReps_filtered.csv", header=T)

data_filter <-  droplevels(subset(data_full, comm=="ok"))
data_filter$genotype <- as.factor(data_filter$genotype)
data_filter$species <- as.factor(data_filter$species)
data_filter$spgen <- as.factor(paste(data_filter$species, data_filter$genotype))


# boxplot:
windows(14,7); par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
boxplot(lag ~ genotype, data=data_filter, main="messgroesse", xlab="", ylab="", las=3, 
        col=c(rainbow(1),rainbow(9),rainbow(8),rainbow(8),rainbow(7)))

# Model

mod_latenzzeit <- gls(lag ~  genotype + start_date, data=data_filter, 
                      weights = varIdent(form=~1|genotype),
                      control=list(maxIter =500, msMaxIter= 500, 
                                   niterEM=500, msMaxEval=500, opt="nlminb"))

mod_latenzzeit2 <- update(mod_latenzzeit, . ~ 0 + genotype + start_date)
# Sets first alphanumeric variant zero --> all following estimates are the actual mean values of each variant. 

windows(); plot(mod_latenzzeit)

mod_latenzzeit2
lag_dataframe_lsmeans <- as.data.frame(lsmeans(mod_latenzzeit, specs="genotype"))
lag_dataframe_lsmeans <- as.data.frame(lsmeans(mod_latenzzeit, specs=c("genotype","start_date")))


# pseudo R^2:
#install.packages("piecewiseSEM")
rsquared(mod_latenzzeit)


# anova:
anova(mod_latenzzeit)

#############################
# 1.: Which Species differ? #
#############################

coef(mod_latenzzeit2)

CoMaSpecies <- rbind(
  "S. habrochaites - S. lycopersicoides"=    c(0,   0,    0,    0,    0,    0,    0,    0,  1/8,    0,    0,    0,  1/8,  1/8,  1/8,    0,    0, -1/7,  1/8,  1/8,    0,  1/8,    0,    0, -1/7, -1/7, -1/7,    0,  1/8,-1/7,   0,    0,-1/7,-1/7,    0,numeric(11)),
  "S. lycopersicoides - S. pimpinellifolium"=c(0,   0,-1/10,    0,    0,-1/10,-1/10,-1/10,    0,-1/10,    0,-1/10,    0,    0,    0,    0,    0,  1/7,    0,    0,-1/10,    0,    0,    0,  1/7,  1/7,  1/7,-1/10,    0, 1/7,   0,-1/10, 1/7, 1/7,-1/10,numeric(11)),
  "S. pimpinellifolium - S. habrochaites"   =c(0,   0, 1/10,    0,    0, 1/10, 1/10, 1/10, -1/8, 1/10,    0, 1/10, -1/8, -1/8, -1/8,    0,    0,    0, -1/8, -1/8, 1/10, -1/8,    0,    0,    0,    0,    0, 1/10, -1/8,   0,   0, 1/10,   0,   0, 1/10,numeric(11)),
  "S. pennellii - S. habrochaites"         =c(0, 1/9,    0,  1/9,  1/9,    0,    0,    0, -1/8,    0,  1/9,    0, -1/8, -1/8, -1/8,  1/9,    0,    0, -1/8, -1/8,    0, -1/8,  1/9,  1/9,    0,    0,    0,    0, -1/8,   0, 1/9,    0,   0,   0,    0,numeric(11)),
  "S. pennellii - S. lycopersicoides"      =c(0, 1/9,    0,  1/9,  1/9,    0,    0,    0,    0,    0,  1/9,    0,    0,    0,    0,  1/9,    0, -1/7,    0,    0,    0,    0,  1/9,  1/9, -1/7, -1/7, -1/7,    0,    0,-1/7, 1/9,    0,-1/7,-1/7,    0,numeric(11)),
  "S. pennellii - S. pimpinellifolium"     =c(0, 1/9,-1/10,  1/9,  1/9,-1/10,-1/10,-1/10,    0,-1/10,  1/9,-1/10,    0,    0,    0,  1/9,    0,    0,    0,    0,-1/10,    0,  1/9,  1/9,    0,    0,    0,-1/10,    0,   0, 1/9,-1/10,   0,   0,-1/10,numeric(11))
)
colnames(CoMaSpecies) <- names(coef(mod_latenzzeit2))

compSpecies <- glht(mod_latenzzeit2, linfct=CoMaSpecies, df=mod_latenzzeit$dims$N - mod_latenzzeit$dims$p)
summary(compSpecies)

##################################################
# 2.: Which Populations differ (within Species)? #
##################################################

coef(mod_latenzzeit2)

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
  "S.lycopersicoides, LA1964"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 6/7,0,0,0,0,0,0,-1/7,-1/7,-1/7,0,0,-1/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S.lycopersicoides, LA2772"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0, 6/7,-1/7,-1/7,0,0,-1/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S.lycopersicoides, LA2776"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7, 6/7,-1/7,0,0,-1/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S.lycopersicoides, LA2777"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7,-1/7, 6/7,0,0,-1/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S.lycopersicoides, LA2951"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7,-1/7,-1/7,0,0, 6/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S.lycopersicoides, LA4123"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7,-1/7,-1/7,0,0,-1/7,0,0, 6/7,-1/7,0,numeric(11)),
  "S.lycopersicoides, LA4130"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7,-1/7,-1/7,0,0,-1/7,0,0,-1/7, 6/7,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0, 7/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0, 7/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8, 7/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8, 7/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0, 7/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8, 7/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0, 7/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0, 7/8,0,0,0,0,0,0,numeric(11)),
  "S.pimpinellifolium, LA1261"=c(0,0, 9/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA1332"=c(0,0,-1/10,0,0, 9/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA1348"=c(0,0,-1/10,0,0,-1/10, 9/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA1374"=c(0,0,-1/10,0,0,-1/10,-1/10, 9/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA1593"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0, 9/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA1659"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0, 9/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA2347"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0, 9/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA2853"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0, 9/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA2983"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0, 9/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA4713"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0, 9/10,numeric(11))
  
)

colnames(CoMaPops) <- names(coef(mod_latenzzeit2))

compPops <- glht(mod_latenzzeit2, linfct=CoMaPops, df=mod_latenzzeit$dims$N - mod_latenzzeit$dims$p)
summary(compPops)


#####
CoMaSpMeans <- CoMaPops
CoMaSpMeans[CoMaSpMeans==9/10] <- -1/10
CoMaSpMeans[CoMaSpMeans==8/9]  <- -1/9
CoMaSpMeans[CoMaSpMeans==7/8]  <- -1/8
CoMaSpMeans[CoMaSpMeans==6/7]  <- -1/7
CoMaSpMeans[CoMaSpMeans<0] <- CoMaSpMeans[CoMaSpMeans<0]*(-1)
CoMaSpMeans <- unique(CoMaSpMeans)
compSpMeans <- glht(mod_latenzzeit2, linfct=CoMaSpMeans, df=mod_latenzzeit$dims$N - mod_latenzzeit$dims$p)
summary(compSpMeans)
compSpMeans
#####


######################################
#### 3. Get adjusted SPECIES mean #### 
######################################
coef(mod_latenzzeit2)

# maybe compare mean against everything else?
# same concept as above!

CoMaSpecies_2 <- rbind(
  "S. habrochaites"=    c(0,    0,    0,    0,    0,    0,    0,    0,  1/8,    0,    0,    0,  1/8,  1/8,  1/8,    0,   0,    0,  1/8,  1/8,    0,  1/8,    0,    0,    0,    0,    0,    0,  1/8,   0,   0,    0,   0,   0,    0,numeric(11)),
  "S. lycopersicoides"= c(0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,   0,  1/7,    0,    0,    0,    0,    0,    0,  1/7,  1/7,  1/7,    0,    0, 1/7,   0,    0, 1/7, 1/7,    0,numeric(11)),
  "S. pimpinellifolium"=c(0,    0, 1/10,    0,    0, 1/10, 1/10, 1/10,    0, 1/10,    0, 1/10,    0,    0,    0,    0,   0,    0,    0,    0, 1/10,    0,    0,    0,    0,    0,    0, 1/10,    0,   0,   0, 1/10,   0,   0, 1/10,numeric(11)),
  "S. pennellii"       =c(0,  1/9,     0, 1/9,  1/9,    0,    0,    0,    0,    0,  1/9,    0,    0,    0,    0,  1/9, 1/9,    0,    0,    0,    0,    0,  1/9,  1/9,    0,    0,    0,    0,    0,   0, 1/9,    0,   0,   0,    0,numeric(11))
)
colnames(CoMaSpecies_2) <- names(coef(mod_latenzzeit2))

compSpecies_2 <- glht(mod_latenzzeit2, linfct=CoMaSpecies_2, df=mod_latenzzeit$dims$N - mod_latenzzeit$dims$p)
compSpecies_2
summary(compSpecies_2)


#####################
## Binomial values ##
#####################

#data_bin <- read.csv("2023_09_data_bin.csv", header=T)

data_bin <- data_full %>% 
  filter(inoculum =="ss") %>% 
  mutate(infection = ifelse(comm=="ok", 1, 0)) %>% 
  dplyr::select(SampleID,species,genotype,start_date,infection)%>% 
  na.omit()

data_bin$species <- as.factor(data_bin$species)
data_bin$genotype <- as.factor(data_bin$genotype)
data_bin$start_date <- as.factor(data_bin$start_date)


# Model

mod_bin <- glm(infection ~  genotype + start_date, data=data_bin,
               family=binomial(link="logit"))
mod_bin2 <- update(mod_bin, . ~ 0 + genotype + start_date)
# Sets first alphanumeric variant zero --> all following estimates are the actual mean values of each variant. 

#windows(); plot(mod_bin)

mod_bin2
lag_dataframe_lsmeans <- as.data.frame(lsmeans(mod_bin, specs="genotype"))
lag_dataframe_lsmeans <- as.data.frame(lsmeans(mod_bin, specs=c("genotype","start_date")))

logit <- lag_dataframe_lsmeans$lsmean
names(logit) <- lag_dataframe_lsmeans$genotype
exp(logit) / (1 + exp(logit))


# pseudo R^2:
#install.packages("piecewiseSEM")
rsquared(mod_bin)


# anova:
library(car)
Anova(mod_bin)


#############################
# 1.: Which Species differ? #
#############################

coef(mod_bin2)

CoMaSpecies <- rbind(
  "S. habrochaites - S. lycopersicoides"=    c(0,   0,    0,    0,    0,    0,    0,    0,  1/8,    0,    0,    0,  1/8,  1/8,  1/8,    0,    0, -1/7,  1/8,  1/8,    0,  1/8,    0,    0, -1/7, -1/7, -1/7,    0,  1/8,-1/7,   0,    0,-1/7,-1/7,    0,numeric(11)),
  "S. lycopersicoides - S. pimpinellifolium"=c(0,   0,-1/10,    0,    0,-1/10,-1/10,-1/10,    0,-1/10,    0,-1/10,    0,    0,    0,    0,    0,  1/7,    0,    0,-1/10,    0,    0,    0,  1/7,  1/7,  1/7,-1/10,    0, 1/7,   0,-1/10, 1/7, 1/7,-1/10,numeric(11)),
  "S. pimpinellifolium - S. habrochaites"   =c(0,   0, 1/10,    0,    0, 1/10, 1/10, 1/10, -1/8, 1/10,    0, 1/10, -1/8, -1/8, -1/8,    0,    0,    0, -1/8, -1/8, 1/10, -1/8,    0,    0,    0,    0,    0, 1/10, -1/8,   0,   0, 1/10,   0,   0, 1/10,numeric(11)),
  "S. pennellii - S. habrochaites"         =c(0, 1/9,    0,  1/9,  1/9,    0,    0,    0, -1/8,    0,  1/9,    0, -1/8, -1/8, -1/8,  1/9,    0,    0, -1/8, -1/8,    0, -1/8,  1/9,  1/9,    0,    0,    0,    0, -1/8,   0, 1/9,    0,   0,   0,    0,numeric(11)),
  "S. pennellii - S. lycopersicoides"      =c(0, 1/9,    0,  1/9,  1/9,    0,    0,    0,    0,    0,  1/9,    0,    0,    0,    0,  1/9,    0, -1/7,    0,    0,    0,    0,  1/9,  1/9, -1/7, -1/7, -1/7,    0,    0,-1/7, 1/9,    0,-1/7,-1/7,    0,numeric(11)),
  "S. pennellii - S. pimpinellifolium"     =c(0, 1/9,-1/10,  1/9,  1/9,-1/10,-1/10,-1/10,    0,-1/10,  1/9,-1/10,    0,    0,    0,  1/9,    0,    0,    0,    0,-1/10,    0,  1/9,  1/9,    0,    0,    0,-1/10,    0,   0, 1/9,-1/10,   0,   0,-1/10,numeric(11))
)
colnames(CoMaSpecies) <- names(coef(mod_bin2))

compSpecies <- glht(mod_bin2, linfct=CoMaSpecies)
summary(compSpecies)

##################################################
# 2.: Which Populations differ (within Species)? #
##################################################

coef(mod_bin2)

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
  "S.lycopersicoides, LA1964"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 6/7,0,0,0,0,0,0,-1/7,-1/7,-1/7,0,0,-1/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S.lycopersicoides, LA2772"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0, 6/7,-1/7,-1/7,0,0,-1/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S.lycopersicoides, LA2776"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7, 6/7,-1/7,0,0,-1/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S.lycopersicoides, LA2777"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7,-1/7, 6/7,0,0,-1/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S.lycopersicoides, LA2951"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7,-1/7,-1/7,0,0, 6/7,0,0,-1/7,-1/7,0,numeric(11)),
  "S.lycopersicoides, LA4123"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7,-1/7,-1/7,0,0,-1/7,0,0, 6/7,-1/7,0,numeric(11)),
  "S.lycopersicoides, LA4130"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/7,0,0,0,0,0,0,-1/7,-1/7,-1/7,0,0,-1/7,0,0,-1/7, 6/7,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0, 7/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0, 7/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8, 7/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8, 7/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0, 7/8,-1/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8, 7/8,0,-1/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0, 7/8,0,0,0,0,0,0,-1/8,0,0,0,0,0,0,numeric(11)),
  "S.habrochaites, LA1559"= c(0,0,0,0,0,0,0,0,-1/8,0,0,0,-1/8,-1/8,-1/8,0,0,0,-1/8,-1/8,0,-1/8,0,0,0,0,0,0, 7/8,0,0,0,0,0,0,numeric(11)),
  "S.pimpinellifolium, LA1261"=c(0,0, 9/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA1332"=c(0,0,-1/10,0,0, 9/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA1348"=c(0,0,-1/10,0,0,-1/10, 9/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA1374"=c(0,0,-1/10,0,0,-1/10,-1/10, 9/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA1593"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0, 9/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA1659"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0, 9/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA2347"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0, 9/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA2853"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0, 9/10,0,0,0,-1/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA2983"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0, 9/10,0,0,-1/10,numeric(11)),
  "S.pimpinellifolium, LA4713"=c(0,0,-1/10,0,0,-1/10,-1/10,-1/10,0,-1/10,0,-1/10,0,0,0,0,0,0,0,0,-1/10,0,0,0,0,0,0,-1/10,0,0,0,-1/10,0,0, 9/10,numeric(11))
  
)

colnames(CoMaPops) <- names(coef(mod_bin2))

compPops <- glht(mod_bin2, linfct=CoMaPops)
summary(compPops)


#####
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
exp(logitSp) / (1 + exp(logitSp))
#####







# The statistical software R (2023) was used to evaluate the data. The data evaluation startet with 
# the definition of an appropriate statistical model based on generalized least squares (Carroll and Ruppert, 1988). (Messgrößen?)
# anda generalized linear model (McCullagh and Nelder, 1989) for (Messgröße?). 
# These models included the factors genotype and start date (without interaction effect). For (normalvert. 
# Messgrößen?), the residuals were assumed to be approximately normally distributed and to be heteroscedastic with 
# respect to the different genotypes. These assumptions are based on a graphical residual analysis. For (nicht-
# normalvert. Messgrößen?), the residuals were binomialy distributed.
# Based on these models, a Pseudo R^2 was calculated (Nakagawa and Schielzeth, 2013) and an analysis of 
# variances (ANOVA) was conducted, followed by multiple contrast tests (Hothorn et 
# al., 2008; see also Bretz et al., 2011). User-defined contrast matrices were used i) to compare the species means with each other, 
# and ii) to compare the population means within their specific species with the corresponding species mean.
#
#
# R Core Team (2023). R: A language and environment for statistical computing. R Foundation for 
# Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, URL http://www.R-project.org/
#
# Carroll, R. J. and Ruppert, D. (1988). Transformation and Weighting in Regression, Chapman and Hall. 
#
# McCullagh P. and Nelder, J. A. (1989). Generalized Linear Models, Chapman and Hall.
#
# Nakagawa, S. and Schielzeth, H. (2013). A General and Simple Method for Obtaining R2 from Generalized 
# Linear Mixed-Effects Models. Edited by Robert B. O’Hara. Methods in Ecology and Evolution 4, no. 2: 
# 133–42. 
#
# Hothorn, T., Bretz, F., and Westfall, P. (2008). Simultaneous Inference in General Parametric Models, 
# Biometrical Journal, 50(3), 346-363.
#
# Bretz, F., Hothorn, T., and Westfall, P. (2011). Multiple Comparisons Using R, Chapman and Hall/CRC, 
# London, isbn 9781584885740.
#

# Arguments for not plotting the IF
# data set is not crossvalidated, and highly unbalanced. Since 2 out of 10 is something differnet than 20 100.
# so in general it is difficult to compare. Further, cant we easily include the estimates from the model
# model outputs logits, but we can't re-convert the SE to a plausible size