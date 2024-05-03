### STATSTICS LEAF SIZE
### 2023_10_25
rm(list=ls())

setwd("C:/Users/suaph281/Nextcloud/ResiDEvo/")
#data_full_adj <- read.csv("C:/Users/Sever/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/2023_09_08AllSlopes_Ann_wthMock_coContaminants_threeReps_filtered_adjusted.csv")
data_full  <- read.csv("data/2023_ScleroPhenotypes/2023_09_08AllSlopesCounts_Ann_wthMock_coContaminants_threeReps_filtered.csv",
                       header=T)


library(gplots)
library(lsmeans)
library(nlme)
library(tidyverse)

data_filter <- data_full %>% filter(comm=="ok") %>% 
  group_by(SampleID) %>% 
  summarise(leaf_size = mean(leaf_area),
            SampleID = SampleID, 
            genotype = genotype, 
            start_date = start_date, 
            species= species) %>% 
  unique()

data_filter$genotype <- as.factor(data_filter$genotype)
data_filter$species <- as.factor(data_filter$species)

# boxplot:
windows(14,7); par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
boxplot(leaf_size ~ genotype, data=data_filter, main="messgroesse", xlab="", ylab="", las=3, 
        col=c(rainbow(1),rainbow(9),rainbow(8),rainbow(8),rainbow(7)))



mod_leaf_size <- gls(leaf_size ~  genotype + start_date, data=data_filter, 
                     weights = varIdent(form=~1|genotype),
                     control=list(maxIter =500, msMaxIter= 500, 
                                  niterEM=500, msMaxEval=500, opt="nlminb"))

mod_leaf_size2 <- update(mod_leaf_size, . ~ 0 + genotype + start_date)
# Sets first alphanumeric variant zero --> all following estimates are the actual mean values of each variant. 

windows(); plot(mod_leaf_size)

mod_leaf_size2

leaf_size_dataframe_lsmeans <- as.data.frame(lsmeans(mod_leaf_size2, specs="genotype"))

### get genotype meanestimate

write.csv(leaf_size_dataframe_lsmeans, "C:/Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/2023_10_25_leaf_size_Pops_Means.csv")

# pseudo R^2:
#install.packages("piecewiseSEM")
library(piecewiseSEM)
rsquared(mod_leaf_size)


# anova:
anova(mod_leaf_size)


library(multcomp)

# 1.: Which Species differ?
###########################

coef(mod_leaf_size2)




CoMaSpecies <- rbind(
  "S. habrochaites - S. lycopersicoides"=    c(0,   0,    0,    0,    0,    0,    0,    0,  1/8,    0,    0,    0,  1/8,  1/8,  1/8,    0,    0, -1/7,  1/8,  1/8,    0,  1/8,    0,    0, -1/7, -1/7, -1/7,    0,  1/8,-1/7,   0,    0,-1/7,-1/7,    0,numeric(11)),
  "S. lycopersicoides - S. pimpinellifolium"=c(0,   0,-1/10,    0,    0,-1/10,-1/10,-1/10,    0,-1/10,    0,-1/10,    0,    0,    0,    0,    0,  1/7,    0,    0,-1/10,    0,    0,    0,  1/7,  1/7,  1/7,-1/10,    0, 1/7,   0,-1/10, 1/7, 1/7,-1/10,numeric(11)),
  "S. pimpinellifolium - S. habrochaites"   =c(0,   0, 1/10,    0,    0, 1/10, 1/10, 1/10, -1/8, 1/10,    0, 1/10, -1/8, -1/8, -1/8,    0,    0,    0, -1/8, -1/8, 1/10, -1/8,    0,    0,    0,    0,    0, 1/10, -1/8,   0,   0, 1/10,   0,   0, 1/10,numeric(11)),
  "S. pennellii - S. habrochaites"         =c(0, 1/9,    0,  1/9,  1/9,    0,    0,    0, -1/8,    0,  1/9,    0, -1/8, -1/8, -1/8,  1/9,    0,    0, -1/8, -1/8,    0, -1/8,  1/9,  1/9,    0,    0,    0,    0, -1/8,   0, 1/9,    0,   0,   0,    0,numeric(11)),
  "S. pennellii - S. lycopersicoides"      =c(0, 1/9,    0,  1/9,  1/9,    0,    0,    0,    0,    0,  1/9,    0,    0,    0,    0,  1/9,    0, -1/7,    0,    0,    0,    0,  1/9,  1/9, -1/7, -1/7, -1/7,    0,    0,-1/7, 1/9,    0,-1/7,-1/7,    0,numeric(11)),
  "S. pennellii - S. pimpinellifolium"     =c(0, 1/9,-1/10,  1/9,  1/9,-1/10,-1/10,-1/10,    0,-1/10,  1/9,-1/10,    0,    0,    0,  1/9,    0,    0,    0,    0,-1/10,    0,  1/9,  1/9,    0,    0,    0,-1/10,    0,   0, 1/9,-1/10,   0,   0,-1/10,numeric(11))
)
colnames(CoMaSpecies) <- names(coef(mod_leaf_size2))

compSpecies <- glht(mod_leaf_size2, linfct=CoMaSpecies, df=mod_leaf_size$dims$N - mod_leaf_size$dims$p)
summary(compSpecies)


# 2.: Which Populations differ (within Species)?
################################################

coef(mod_leaf_size2)

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

colnames(CoMaPops) <- names(coef(mod_leaf_size2))

compPops <- glht(mod_leaf_size2, linfct=CoMaPops, df=mod_leaf_size$dims$N - mod_leaf_size$dims$p)
summary(compPops)


#### OUTPUT stats to table 

table_compPops <- broom::tidy(compPops)
table_comp_Species <- broom::tidy(compSpecies)

write.table(table_compPops, "C:/Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/2023_10_09_leaf_size_CompPops.csv")
write.table(table_comp_Species, "C:/Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/2023_10_09_leaf_size_CompSpecies.csv")


######################################
#### 3. Get adjusted SPECIES mean #### 
######################################
coef(mod_leaf_size2)

# WORKS

CoMaSpecies_2 <- rbind(
  "S. habrochaites"=    c(0,    0,    0,    0,    0,    0,    0,    0,  1/8,    0,    0,    0,  1/8,  1/8,  1/8,    0,   0,    0,  1/8,  1/8,    0,  1/8,    0,    0,    0,    0,    0,    0,  1/8,   0,   0,    0,   0,   0,    0,numeric(11)),
  "S. lycopersicoides"= c(0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,   0,  1/7,    0,    0,    0,    0,    0,    0,  1/7,  1/7,  1/7,    0,    0, 1/7,   0,    0, 1/7, 1/7,    0,numeric(11)),
  "S. pimpinellifolium"=c(0,    0, 1/10,    0,    0, 1/10, 1/10, 1/10,    0, 1/10,    0, 1/10,    0,    0,    0,    0,   0,    0,    0,    0, 1/10,    0,    0,    0,    0,    0,    0, 1/10,    0,   0,   0, 1/10,   0,   0, 1/10,numeric(11)),
  "S. pennellii"       =c(0,  1/9,     0, 1/9,  1/9,    0,    0,    0,    0,    0,  1/9,    0,    0,    0,    0,  1/9, 1/9,    0,    0,    0,    0,    0,  1/9,  1/9,    0,    0,    0,    0,    0,   0, 1/9,    0,   0,   0,    0,numeric(11))
)
colnames(CoMaSpecies_2) <- names(coef(mod_leaf_size2))

compSpecies_2 <- glht(mod_leaf_size2, linfct=CoMaSpecies_2, df=mod_leaf_size$dims$N - mod_leaf_size$dims$p)
compSpecies_2
summary(compSpecies_2)

data_id <- data_full %>%
  filter(genotype!="LA2727") %>% 
  dplyr::select(species, genotype) %>% 
  unique() %>% 
  mutate(contrast = species) %>% 
  left_join(broom::tidy(compSpecies_2), by="contrast") %>% 
  dplyr::select(species, estimate) %>% 
  unique()

write.csv(data_id,"C:/Users/suaph281/Nextcloud/ResiDEvo/2023_ScleroSolanumAuswertung/2023_PaperVizStats/stats_tables/2023_10_10_leaf_size_Species_Means.csv")


## 
# Get pairwise comparisons of s. pennellii accessions
##


# get all possible combinations of S. pennellii accessions 

spen <- data_filter %>% 
  ungroup() %>% 
  dplyr::select(genotype, species) %>% 
  filter(species=="S. pennellii") %>% 
  unique() 

comb <- expand.grid(spen$genotype,spen$genotype)

combinations <- comb %>% 
  unite(pairs, Var1, Var2, sep=" - ")

# n= sample size per group
# = sample size per accession 

nn = data_filter %>% 
  dplyr::select(genotype, SampleID) %>%
  group_by(genotype) %>% 
  summarise(n = n(), 
            genotype=genotype) %>% 
  unique()

nnn = nn$n
names(nnn) = nn$genotype

KoMaPFsub <- contrMat(n=nnn, type="Tukey")
combina <- combinations$pairs

KoMaPF <- KoMaPFsub[row.names(KoMaPFsub) %in% combina,]

new_columns <- matrix(0, nrow = nrow(KoMaPF), ncol = 11)
result_matrix <- cbind(KoMaPF, new_columns)
colnames(result_matrix) = names(coef(mod_leaf_size2))
# Get pairwise stats
comp_spen_leaf_size = glht(mod_leaf_size2, linfct=result_matrix)
summary(comp_spen_leaf_size)

#get CLD

analysis <-   comp_spen_leaf_size %>% 
  broom::tidy() %>% 
  mutate(contrast = str_replace_all(contrast, "\\s", "")) %>% 
  pull(adj.p.value, contrast)  %>% 
  multcompView::multcompLetters()

letters_tbl <- analysis$Letters %>% 
  enframe("trt", "label") %>% 
  mutate(rsp = 0)

# add population mean 
# wont do this 
estimates <- data.frame(mod_leaf_size2$coefficients[0:35]) %>% 
  rownames_to_column("genotype") %>% 
  mutate(genotype = str_replace(genotype, "genotype", ""))


write.table(letters_tbl, "C://Users/suaph281/Nextcloud/ResiDEvo/data/2023_ScleroPhenotypes/stats_tables/2023_10_26_spen_leafsize_pairwise.txt")