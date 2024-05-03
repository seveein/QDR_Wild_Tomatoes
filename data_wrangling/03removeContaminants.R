## REMOVE CONTAMINANTS

# load file that contains contaminants-id

exclude <- read.csv2("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/2023_SolSclerSkip.csv", header=T)
exclude$species[exclude$species=="S.pimpinellifolium"]="S. pimpinellifolium"
exclude$species[exclude$species=="S.lycopersicoides"]="S. lycopersicoides"

#### SLOPE ####
# load slope file and edit

data_slope <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16AllSlopes_Ann.csv", header=T)
#data_slope$species[data_slope$species=="habrochaites"]="S. habrochaites"

data_cleaned_slope <- data_slope %>% 
  filter(!(species%in%exclude$species & start_date %in% exclude$experiment & genotype%in%exclude$exclude))  

write.csv(x= data_cleaned_slope, file= "/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16AllSlopes_Ann_noContaminants.csv")


#### COUNTS ####

data_count <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16AllCounts_Ann_wthMOCK.csv", header=T)

data_cleaned_counts <- data_count %>% 
  filter(!(species%in%exclude$species & start_date %in% exclude$experiment & genotype%in%exclude$exclude))

write.csv(x= data_cleaned_counts, file= "/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16AllCounts_Ann_wthMOCK_noContaminants.csv")


#### Counts Slopes ####
data_count_slopes <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16AllSlopesCounts_Ann_wthMOCK.csv", header=T)

data_cleaned_counts_slopes <- data_count_slopes %>% 
  select(!species.y & !inoculum.y & !indv.y & !X & !X.1 & !genotype.y) %>% 
  filter(!(species.x%in%exclude$species & start_date %in% exclude$experiment & genotype.x%in%exclude$exclude))
head(data_cleaned_counts_slopes)

colnames(data_cleaned_counts_slopes) = c("ID", "time", "lesion_area", "leaf_area", "box", "start_date", "species", "genotype", "indv", "inoculum", "lag", "slope", "comm", "LDT")

write.csv(x= data_cleaned_counts, file= "/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16AllSlopesCounts_Ann_wthMock_coContaminants.csv")
