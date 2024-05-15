
## Dan Watson 
## Last edit: 2024-08-02
## Optimized script for formatting weather station data and merging with ARU data

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(pbapply)
library(lubridate)


# Load weather data
weather <- read_csv("data/brunswick_weather_clean.csv")
# Load ARU data
aru_data <- read_csv("data/acoustic_dat_24_07.csv")

# Convert date and time to POSIXct datetime format in weather
weather <- weather %>%
  mutate(date_time = as.POSIXct(paste(date, time), format = "%Y-%m-%d %H:%M:%S"))

# Convert date and time to POSIXct datetime format in aru_data
aru_data <- aru_data %>%
  mutate(date_time = as.POSIXct(paste(date, time), format = "%Y-%m-%d %H:%M:%S"))

# Filter weather data by the desired date range
weather <- weather %>%
  filter(date >= as.Date("2024-03-10") & date <= as.Date("2024-06-28"))

# Filter aru_data by the desired date range
aru_data <- aru_data %>%
  filter(date >= as.Date("2024-03-10") & date <= as.Date("2024-06-28"))

# Remove rows with NAs in the 'site' column
aru_data <- aru_data %>%
  filter(!is.na(site))

# Find the closest date_time in weather for each date_time in aru_data
closest_weather <- aru_data %>%
  mutate(
    closest_weather_index = map_dbl(
      date_time, 
      ~ which.min(abs(difftime(weather$date_time, ., units = "secs")))
    )
  ) %>%
  left_join(weather, by = c("closest_weather_index" = "row_number"))

# Merge with aru_data
merged_data <- aru_data %>%
  left_join(closest_weather, by = "date_time")

# Function to find the closest weather time for each ARU datetime
find_closest_weather <- function(aru_time) {
  diffs <- abs(difftime(weather$date_time, aru_time, units = "secs"))
  closest_index <- which.min(diffs)
  return(weather[closest_index, ])
}

# Apply the function across all rows in ARU data with a progress bar
closest_weather <- pblapply(aru_data$date_time, find_closest_weather)

# Combine the closest weather data with ARU data
closest_weather <- bind_rows(closest_weather)
merged_data <- bind_cols(aru_data, closest_weather)



## Old version that works but is not fast enough, but neither is the above 

## Dan Watson 
## Last edit: 2024-08-02
## Formatting weather station data from Brunswick, GA (KBQK) for use in the ARU analysis
## Downloaded from MesoWest (https://mesowest.utah.edu/)

# Clear workspace 
rm(list=ls())

# Load packages 
library(tidyverse)

# Load weather data
weather <- read_csv("data/brunswick_weather_clean.csv")

# Create new columns for date and time from Date_Time column by extracting the first 10 characters for date and characters 12-16 for time 
weather <- weather %>% 
  mutate(date = substr(Date_Time, 1, 10),
         time = substr(Date_Time, 12, 16)) %>%
  select(-Date_Time)

# NAs in precipitation column are replaced with 0
weather$precip_accum_inches[is.na(weather$precip_accum_inches)] <- 0

#check structure of weather data
str(weather)
#make sure date structure is correct 
weather$date <- as.Date(weather$date, format = "%Y-%m-%d")

# Load aru timeline
aru_timeline <- read_csv("data/aru_timeline.csv")

# Load aru data
aru_data <- read_csv("data/acoustic_dat_24_07.csv")

# Remove weather rows before 2024-03-10 and after 2024-06-29
weather <- weather %>% 
  filter(date >= "2024-03-10") %>%
  filter(date <= "2024-06-29")
aru_data <- aru_data %>% 
  filter(date >= "2024-03-10") %>%
  filter(date <= "2024-06-29")

# Convert date and time to POSIXct for accurate time comparison
aru_data <- aru_data %>%
  mutate(datetime = as.POSIXct(paste(date, time), format = "%Y-%m-%d %H:%M:%S"))

weather <- weather %>%
  mutate(datetime = as.POSIXct(paste(date, time), format = "%Y-%m-%d %H:%M:%S"))

# Function to find the closest time match
find_closest_time <- function(datetime) {
  diffs <- abs(difftime(weather$datetime, datetime, units = "secs"))
  closest_index <- which.min(diffs)
  return(weather[closest_index, ])
}

# Apply the function to each row of aru_data
merged_df <- aru_data %>%
  rowwise() %>%
  mutate(
    closest_weather = list(find_closest_time(datetime))
  ) %>%
  unnest_wider(closest_weather, names_sep = "_")  # Expanding the nested list column



 