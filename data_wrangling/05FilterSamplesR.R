## FILTER for samples with variable leaf size
# stand: 2023_09_08 

rm(list=ls())

library(tidyverse)

setwd("C://Users/suaph281/Nextcloud/ResiDEvo/")
data_raw <- read.csv("data/2023_ScleroPhenotypes/2023_08_28AllSlopesCounts_Ann_wthMock_coContaminants_threeReps.csv", header=T) 

colnames(data_raw)=c("X.1", "X", "id", "time", "lesion_area", "leaf_area", "box", 
                      "start_date", "species", "genotype", "indv.", "inoculum", "lag", "slope",
                      "comm", "LDT")

data_raw <- data_raw %>% 
  select(!X.1 & !X)

# remove excessive data points (1.5 % IQR) 
# and remove odd samples with variable leaf size (variability of more than 15 % of leaf size)
# THIS IS FINAL FILTER:

get_filter <- function(data){
  out <- data %>%
    filter(time < 990 & time > 20, 
           comm!="s"|is.na(comm)==T) %>%
    # get individual sample ID
    mutate(SampleID = paste(start_date, genotype, id, inoculum, sep="_")) %>%
    group_by(SampleID) %>% 
    mutate(
      mean = mean(leaf_area),
      p2.5 = quantile(leaf_area, 0.015),
      p97.5 = quantile(leaf_area, 0.985),
      IQR = p97.5 - p2.5,
      IQR_mean = IQR/mean) %>%
    filter(leaf_area < p97.5 & leaf_area >p2.5) %>% 
    group_by(SampleID) %>%
    mutate(mean_leaf =  mean(leaf_area),
           sd_leaf = sd(leaf_area),
           ratio = sd_leaf/mean_leaf) %>%
    filter(ratio < .15) %>%
    return(out)
}

length <- dim(data_raw)[1]

##### TEST EACH FILTER STEP #####

# to trust the data, check, how much is removed by each step:

t990 <- data_raw %>% 
  filter(time > 990) %>% 
  summarize(n=n()) %>% 
  mutate(diff=length-n)
t30 <- data_raw %>% 
  filter(time < 20) %>% 
  summarize(n=n()) %>% 
  mutate(diff=length-n)
skip <- data_raw %>% 
  filter(comm=="s") %>% 
  summarize(n=n()) %>% 
  mutate(diff=length-n)

# removed just by hard filtering:

(hard_filter = (t990$n + t30$n + skip$n))



data_points_weid <- data_raw %>%
  filter(time < 990 & time > 20) %>% 
  filter(comm!="s"|is.na(comm)==T) %>%
    # get individual sample ID
    mutate(SampleID = paste(start_date, genotype, id, inoculum, sep="_")) %>%
    group_by(SampleID) %>% 
    mutate(
      mean = mean(leaf_area),
      p2.5 = quantile(leaf_area, 0.02),
      p97.5 = quantile(leaf_area, 0.98),
      IQR = p97.5 - p2.5,
      IQR_mean = IQR/mean) %>%
    filter(leaf_area > p97.5 | leaf_area < p2.5)  %>% 
  ungroup() %>% 
  summarize(n=n()) %>% 
  select(n) 

(soft_filter = data_points_weid$n)
# removes 323.000 data points
hard_filter + soft_filter

###### 'ON FINAL DATA' #####

data_filtered <- get_filter(data_raw)


# get the SampleIDs that are still in the panel
SampleIDs_kept <- data_filtered %>% 
  select(SampleID) %>% 
  unique()

# ANTI join so you have the samples that are excluded form the panel!
SampleIDs_filtered <- data_raw %>%
  mutate(SampleID = paste(start_date, genotype, id, inoculum, sep="_")) %>%
  select(SampleID) %>% 
  unique() %>% 
  anti_join(SampleIDs_kept)
  
# skipped samples are the most frequent reason for removal! not filtering.

#####
## Edit the actual data files 
#####

slope <- data_filtered %>% 
  select(SampleID,start_date, species, genotype, inoculum, comm, slope, lag, LDT) %>% 
  unique()

write.csv(slope,"data/2023_ScleroPhenotypes/2023_09_08AllSlopes_Ann_wthMock_coContaminants_threeReps_filtered.csv" )


count_slope <- data_filtered %>% 
  select(SampleID,start_date, species, genotype, inoculum, comm, time, leaf_area, lesion_area, slope, lag, LDT) 

write.csv(count_slope,"data/2023_ScleroPhenotypes/2023_09_08AllSlopesCounts_Ann_wthMock_coContaminants_threeReps_filtered.csv" )

