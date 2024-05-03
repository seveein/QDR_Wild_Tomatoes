#### Focus on the three most plausible repetitions per species
rm(list=ls())

library(tidyverse)

### slopes ###
data_slope <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16AllSlopes_Ann_noContaminants.csv", header=T)

data_slope_threeReps <-data_slope %>% 
  filter(!start_date %in% c("2023_03_24", "2023_03_28", 
                            "2023_05_04", "2023_05_30", 
                            "2023_03_31"))

write.csv(data_slope_threeReps,"/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_28AllSlopes_Ann_noContaminants_threeReps.csv" )

#### counts ####

data_counts <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16AllCounts_Ann_wthMOCK_noContaminants.csv", header=T)

data_counts_threeReps <-data_counts %>% 
  filter(!start_date %in% c("2023_03_24", "2023_03_28",  
                            "2023_05_04", "2023_05_30", 
                            "2023_03_31"))

write.csv(data_counts_threeReps,"/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_28AllCounts_Ann_wthMOCK_noContaminants_threeReps.csv" )

#### counts & slopes ####

data_counts_slopes <- read.csv("/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_16AllSlopesCounts_Ann_wthMock_coContaminants.csv", header=T)

data_counts_slopes_threeReps <-data_counts_slopes %>% 
  filter(!start_date %in% c("2023_03_24", "2023_03_28", 
                            "2023_05_04", "2023_05_30", 
                            "2023_03_31"))

write.csv(data_counts_slopes_threeReps,"/home/se_residevo/Desktop/SEVERIN LAB/ResiDEvo/data/2023_ScleroPhenotypes/2023_08_28AllSlopesCounts_Ann_wthMock_coContaminants_threeReps.csv" )
