#Dan Watson, danwatson@vt.edu
#latest update: 05/15/2024
#working on process to analyze birdnet data collected with wildlife acoustic recorders from 2023 field season

#clear workspace
rm(list=ls())

#load packages 
library(tidyverse)
library(tidyr)
library(dplyr)

### start here to upload unprocessed birdnet data from txt files downloaded from the ARC, NAS, or other hardrive ###

#unzip folder from local_data folder already within working directory 
unzip("data/data_2023_07_output.zip")

#string together and load the unzipped txt files and name the dataframe "acoustic_dat" that can be downloaded as a single large dataset in .csv format
#this will take awhile
acoustic_dat <- list.files(pattern = ".txt") %>% 
  map_dfr(~read_tsv(.x, col_names = FALSE), .id = "file") %>% #read each file and combine them
  rename_all(tolower) %>% #convert all column names to lowercase
  rename_all(~str_replace_all(., "V", "col")) %>% #replace V with col in column names 
  rename_all(~str_replace_all(., "X", "row")) %>% #replace X with row in column names
  mutate(file = str_remove(file, ".txt")) %>% #remove .txt from file column
  mutate(file = as.numeric(file)) #convert file column to numeric

#change column names x1 = "selection", x2 = "view", x3 = "channel", x4 = "begin_file", x5 = "begin_time", x6 = "end_time", x7 = "low_freq", x8 = "high_freq", x9 = "species_code", x10 = "common_name", x11 = "confidence", x12 = "rank" 
colnames(acoustic_dat) <- c("file", "selection", "view", "channel", "begin_file", "begin_time", "end_time", "low_freq", "high_freq", "species_code", "common_name", "confidence", "rank") #file column remains "file" 

#remove old headers left over from the txt files by filtering out rows that do not contain "Spectrogram 1" in the 'view' column
acoustic_dat <- acoustic_dat %>% 
  filter(view == "Spectrogram 1")

#remove un-needed columns for analysis if you'd like
acoustic_dat <- acoustic_dat %>% 
  select(-c(file, selection, view, channel, low_freq, high_freq, species_code))

#save the dataframe as a .csv file to avoid re-running previous code
#this example .csv is already uploaded to our github repo
write_csv(acoustic_dat, "data/acoustic_dat_07.csv")

### start here if you already have some birdnet data processed to .csv format, perhaps our example for the github repo ###

#read in csv files
acoustic_dat <- read_csv("data/acoustic_dat_07.csv") 

#change 'species' to 'common_name' for matching up with other .csv with alpha codes-- my mistake
acoustic_dat <- acoustic_dat %>% 
  rename(common_name = species)

sp_codes <- read_csv("data/species.list.csv") #for alpha codes and scientific names

#add column for sp_codes and scienfitic_name to acoustic_dat that contains the 'sp_code' (aka alpha codes) for each corresponding 'common_name'
acoustic_dat <- acoustic_dat %>% 
  left_join(sp_codes, by = "common_name")

#view the first few rows of the dataframe
head(acoustic_dat)

#extract date and time from the 'begin_file' column
acoustic_dat <- acoustic_dat %>%
  mutate(
    date = sub(".*_(\\d{8})_.*", "\\1", begin_file),
    time = sub(".*_(\\d{6})\\.wav", "\\1", begin_file)
  )

#convert date to MMDDYYYY structure and then format as date instead of chr
acoustic_dat$date <- format(as.Date(acoustic_dat$date, format="%Y%m%d"), "%m%d%Y")
acoustic_dat$date <- as.Date(acoustic_dat$date, format="%m%d%Y") 

#convert time to HH:MM:SS format
acoustic_dat$time <- format(strptime(acoustic_dat$time, "%H%M%S"), "%H:%M:%S")

str(acoustic_dat) #check structure of dataframe

#create new dataframe with rows only beginning with "SMA10458' in the 'begin_file' column for analysis of this specific site 'R-3' and its associated ARU name 'SMA10458'
r3 <- acoustic_dat %>% 
  filter(str_detect(begin_file, "SMA10458")) %>% #add column 'site' with value 'R-3' for every row
  mutate(site = "R-3") 

#view the dataframe
View(r3)

## confidence score filtering ##
#filter out all observations with a confidence score less than 0.5 for now
r3 <- r3 %>% 
  filter(confidence >= 0.5)

## recording time selection ##
#diurnal/sunrise recordings only
#remove rows with a time of 13:00:00 or greater just for sunrise and diurnal recordings 
r3 <- r3 %>% 
  filter(time < "13:00:00")

## check out the dataset a bit now ##

#look for data gaps in the 'date' column for missing days in the month 
#extract the month and year from the first date
start_date <- min(r3$date)
end_date <- max(r3$date)

#generate a sequence of all dates for the month
all_dates <- seq(from = floor_date(start_date, "month"), to = ceiling_date(end_date, "month") - days(1), by = "day")

#find missing dates
missing_dates <- setdiff(all_dates, unique(acoustic_dat$date))
print(missing_dates)

#look at the number of all unique species detected and display their names with total counts
unique_species <- r3 %>% 
  group_by(common_name) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n))
View(unique_species)

#write the dataframe to a .csv file for further analysis
write_csv(r3, "data/r3_07_23.csv")
