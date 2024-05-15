#Dan Watson, danwatson@vt.edu
#05/10/2024
#working on process to analyze birdnet data collected with wildlife acoustic recorders from 2023 field season

#clear workspace
rm(list=ls())

#load packages 
library(tidyverse)
library(tidyr)
library(dplyr)

### start here to upload unprocessed birdnet data from txt files downloaded from the ARC, NAS, or other hardrive ###

#unzip folder from working directory 
#this is just the data for July 2023, but should work for an entire dataset
unzip("data_2023_07_output.zip")  

#string together and load the unzipped txt files and name the dataframe "acoustic_dat" that can be downloaded as a single large dataset in .csv format
#this will take awhile
acoustic_dat <- list.files(pattern = ".txt") %>% 
  map_dfr(~read_tsv(.x, col_names = FALSE), .id = "file") %>% #read each file and combine them
  rename_all(tolower) %>% #convert all column names to lowercase
  rename_all(~str_replace_all(., "V", "col")) %>% #replace V with col in column names 
  rename_all(~str_replace_all(., "X", "row")) %>% #replace X with row in column names
  mutate(file = str_remove(file, ".txt")) %>% #remove .txt from file column
  mutate(file = as.numeric(file)) #convert file column to numeric

#change column names x1 = "selection", x2 = "view", x3 = "channel", x4 = "begin_file", x5 = "begin_time", x6 = "end_time", x7 = "low_freq", x8 = "high_freq", x9 = "species_code", x10 = "species", x11 = "confidence", x12 = "rank" 
colnames(acoustic_dat) <- c("file", "selection", "view", "channel", "begin_file", "begin_time", "end_time", "low_freq", "high_freq", "species_code", "species", "confidence", "rank") #file column remains "file" 

#remove old headers left over from the txt files by filtering out rows that do not contain "Spectrogram 1" in the 'view' column
acoustic_dat <- acoustic_dat %>% 
  filter(view == "Spectrogram 1")

#remove un-needed columns for analysis if you'd like
acoustic_dat <- acoustic_dat %>% 
  select(-c(file, selection, view, channel, low_freq, high_freq, species_code))

#save the dataframe as a .csv file to avoid re-running previous code
#this example .csv is already uploaded to our github repo
write_csv(acoustic_dat, "acoustic_dat_07.csv")

### start here if you already have some birdnet data processed to .csv format, perhaps from our github repo ###

#read in csv file
acoustic_dat <- read_csv("acoustic_dat_07.csv") 

#create new dataframe with rows only beginning with "SMA10458' in the 'begin_file' column for analysis of this specific site 'R-3' and its associated ARU name 'SMA10458'
r3_sma10458 <- acoustic_dat %>% 
  filter(str_detect(begin_file, "SMA10458"))

#confidence score filtering 
#filter out all observations with a confidence score less than 0.5 for now
r3_sma10458 <- r3_sma10458 %>% 
  filter(confidence >= 0.5)



