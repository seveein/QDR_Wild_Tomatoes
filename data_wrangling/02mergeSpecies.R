# Combine all annotated slope-files into one big file
# CONTAMINANTS are INcluded
rm(list=ls())

#### SLOPE pure ####

library(tidyverse)
library(ggpubr)

shabro <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/habrochaites/2023_06_shabro_sclero_ann_ALL.csv")
shabro <- shabro %>% 
  select(!X.1 & !X)
spimpi <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/pimpinellifolium/2023_06_spimp_sclero_ann_ALL.csv")
spimpi <- spimpi%>% 
  select(!X.1 & !X)
spenne <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/pennellii/2023_06_19_spen_sclero_ann_ALL.csv")
spenne <- spenne%>% 
  select(!X.1 & !X)
slycop <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/lycopersicoides/2023_06_19_slycop_sclero_ann_ALL.csv")
slycop <- slycop%>% 
  select(!X.1 & !X)

merge <- shabro %>% 
  bind_rows(slycop, spenne, spimpi)

merge$species[merge$species == "pennellii "]="pennellii"

merge$species[merge$species == "pennellii"]="S. pennellii"
merge$species[merge$species == "pimpinellifolium"]="S. pimpinellifolium"
merge$species[merge$species == "habrochaites"]="S. habrochaites"
merge$species[merge$species == "lycopersicoides"]="S. lycopersicoides"
merge$species[merge$species == "lycopersicum"]="S. lycopersicum"

write.csv(merge, "/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16AllSlopes_Ann.csv")


#### COUNTS ####

shabro <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/habrochaites/2023_08_16_shabro_sclero_COUNTS_ann_ALL_wthMOCK.csv")
spimpi <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/pimpinellifolium/2023_08_16_spimp_sclero_COUNTS_ann_wthMOCK_ALL.csv")
spimpi <- spimpi %>% 
  select(!X.1 & !X)
spenne <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/pennellii/2023_08_16_spen_sclero_COUNTS_ann_ALL_edited_wthMOCK.csv")
spenne <- spenne %>% 
  select(!X.1 & !X) 
slycop <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/lycopersicoides/2023_08_16_slycop_sclero_COUNTS_wthMOCK_ALL.csv")
slycop <- slycop %>% 
  select(!X.1 & !X)

head(shabro)

merge <- shabro %>% 
  select(!X.1 & !X) %>% 
  bind_rows(slycop, spenne, spimpi)

merge$species[merge$species == "pennellii "]="pennellii"

merge$species[merge$species == "pennellii"]="S. pennellii"
merge$species[merge$species == "pimpinellifolium"]="S. pimpinellifolium"
merge$species[merge$species == "habrochaites"]="S. habrochaites"
merge$species[merge$species == "lycopersicoides"]="S. lycopersicoides"
merge$species[merge$species == "lycopersicum"]="S. lycopersicum"

#write.csv(merge, "/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16AllCounts_Ann_wthMOCK.csv")

#### SLOPES and COUNTS ####
merge <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16AllCounts_Ann_wthMOCK.csv", header=T)
slopes <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16AllSlopes_Ann.csv", header=T)
slopes <- slopes %>% 
  select(!X)

merge_counts_clopes <- merge %>% 
  left_join(slopes, by=c("start_date", "Box", "ID"))

write.csv(merge_counts_clopes, "/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16AllSlopesCounts_Ann_wthMOCK.csv")
