
## Dan Watson 
## Last edit: 2024-08-02
## Optimized script for formatting weather station data and merging with ARU data

# Load packages 
library(data.table) # Using data.table for faster data processing

# Load weather data (only necessary columns)
weather <- fread("data/brunswick_weather_clean.csv", select = c("Date_Time", "precip_accum_inches"))

# Create new columns for date, time, and handle missing precipitation data
weather[, `:=`(date = as.Date(substr(Date_Time, 1, 10)),
               time = substr(Date_Time, 12, 16),
               precip_accum_inches = fifelse(is.na(precip_accum_inches), 0, precip_accum_inches))]

# Convert Date_Time to POSIXct datetime format
weather[, datetime := as.POSIXct(Date_Time, format = "%Y-%m-%d %H:%M:%S")]

# Filter weather data by the desired date range
weather <- weather[date >= as.Date("2024-03-10") & date <= as.Date("2024-06-29")]

# Load aru data (only necessary columns)
aru_data <- fread("data/acoustic_dat_24_07.csv", select = c("date", "time", "other_needed_columns"))

# Convert date and time to POSIXct datetime format in aru_data
aru_data[, datetime := as.POSIXct(paste(date, time), format = "%Y-%m-%d %H:%M:%S")]

# Filter aru_data by the desired date range
aru_data <- aru_data[date >= as.Date("2024-03-10") & date <= as.Date("2024-06-29")]

# Perform an efficient join based on nearest datetime using data.table's rolling join
setkey(weather, datetime)
setkey(aru_data, datetime)

# Rolling join to find the closest match for datetime
merged_data <- aru_data[weather, roll = "nearest"]

# Save the merged data (optional step if you want to save the processed data)
# fwrite(merged_data, "data/merged_aru_weather_data.csv")

