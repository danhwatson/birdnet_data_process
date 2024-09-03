# Dan Watson, danwatson@vt.edu
# Latest update: 05/15/2024
# Working on process to analyze birdnet data collected with Wildlife Acoustic recorders (ARUs) from 2023 field season
# Script relies heavily on Wildlife Acoustic's file naming system for their ARUs

# Clear workspace
rm(list=ls())

# Load packages 
library(tidyverse)

### Start here to upload unprocessed birdnet data from txt files downloaded from the ARC, NAS, or other hardrive ###

# Unzip birdnet output folders from each month within working directory 
unzip("data/data.zip")

# Concatenate birdnet files and name the dataframe "acoustic_dat" that can be downloaded as a single large dataset in .csv format
# This may take awhile
acoustic_dat <- list.files(pattern = ".txt") %>% 
  map_dfr(~read_tsv(.x, col_names = FALSE), .id = "file") %>% #read each file and combine them
  rename_all(tolower) %>% #convert all column names to lowercase
  rename_all(~str_replace_all(., "V", "col")) %>% #replace V with col in column names 
  rename_all(~str_replace_all(., "X", "row")) %>% #replace X with row in column names
  mutate(file = str_remove(file, ".txt")) %>% #remove .txt from file column
  mutate(file = as.numeric(file)) #convert file column to numeric

#remove all txt files out birdnet_data_process directory 
file.remove(list.files(pattern = ".txt"))

# Change column names x1 = "selection", x2 = "view", x3 = "channel", x4 = "begin_file", x5 = "begin_time", x6 = "end_time", x7 = "low_freq", x8 = "high_freq", x9 = "species_code", x10 = "common_name", x11 = "confidence", x12 = "rank"
colnames(acoustic_dat) <- c("file", "selection", "view", "channel", "begin_time(s)", "end_time(s)", "low_freq", "high_freq", "common_name", "species_code", "confidence", "file_name", "file_offset(s)") #file column remains "file" 

# Remove old headers left over from the txt files by filtering out rows that do not contain "Spectrogram 1" in the 'view' column
acoustic_dat <- acoustic_dat %>% 
  filter(view == "Spectrogram 1")

#remove "/projects/birdnet/chemours/Data_2024/" from the beginning of 'file_name' rows
acoustic_dat$file_name <- str_remove(acoustic_dat$file_name, "/projects/birdnet/chemours/Data_2024_pc_ac/")

# Remove un-needed columns for analysis if you'd like
acoustic_dat <- acoustic_dat %>% 
  select(-c(file, selection, view, channel, low_freq, high_freq, species_code))

# Extract date and time from the 'begin_file' column
acoustic_dat <- acoustic_dat %>%
  mutate(
    date = sub(".*_(\\d{8})_.*", "\\1", file_name),
    time = sub(".*_(\\d{6})\\.wav", "\\1", file_name)
  )

# Convert date to MMDDYYYY structure and then format as date instead of chr
acoustic_dat$date <- format(as.Date(acoustic_dat$date, format="%Y%m%d"), "%m%d%Y")
acoustic_dat$date <- as.Date(acoustic_dat$date, format="%m%d%Y") 

# Convert time to HH:MM:SS format
acoustic_dat$time <- format(strptime(acoustic_dat$time, "%H%M%S"), "%H:%M:%S")

# Add an "aru" column to the dataframe that contains the ARU from the first 8 characters of the 'begin_file' column
acoustic_dat <- acoustic_dat %>% 
  mutate(aru = substr(file_name, 1, 8))
         
### Start here if you already have some birdnet data processed to .csv format

sp_codes <- read_csv("data/species_list.csv") #for alpha codes and scientific names

# Add column for sp_codes and scientific_name to acoustic_dat that contains the 'sp_code' (aka alpha codes) for each corresponding 'common_name'
acoustic_dat <- acoustic_dat %>% 
  left_join(sp_codes, by = "common_name")

# Rename "sp_code.x" and "scientific_name.x" to "sp_code" and "scientific_name"
acoustic_dat <- acoustic_dat %>% 
  rename(sp_code = sp_code.y, scientific_name = scientific_name.y)

# Read in aru_timeline
aru_timeline <- read_csv("data/aru_timeline.csv")

# Convert the date column in aru_timeline to Date type
aru_timeline$date <- as.Date(aru_timeline$date, format = "%m/%d/%y")

# Transform aru_timeline to long format
aru_long <- aru_timeline %>%
  pivot_longer(cols = -date, names_to = "site", values_to = "aru")

# Merge the datasets based on aru and date
acoustic_dat <- acoustic_dat %>%
  left_join(aru_long, by = c("aru", "date"))

#rearrange columns
acoustic_dat <- acoustic_dat %>% 
  select("site", "aru", "common_name", "scientific_name", "sp_code", "confidence", "date", "time", "begin_time(s)", "end_time(s)", "file_offset(s)", "file_name",)

# Write the dataframe to a .csv file for further analysis
write_csv(acoustic_dat, "data/acoustic_dat_24_pc.csv")

# Load the dataframe
acoustic_dat <- read_csv("data/2024/acoustic_dat_24_pc.csv")


