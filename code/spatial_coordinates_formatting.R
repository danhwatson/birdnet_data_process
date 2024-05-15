



# Clear workspace 
rm(list=ls())

# Load libraries
library(tidyverse)

# Load data
aru_data <- read_csv("data/acoustic_dat_24_07.csv")

cords <- read_csv("data/aru_site_coordinates.csv") #For spOccupancy, you must use projected coordinates (UTM zone 3 for this example)
# It does not use lat/long

# Merge X and Y columns based on site 
aru_data <- aru_data %>% 
  left_join(cords, by = c("site" = "site")) %>% 
  select(-lat, -long)

# Write 
