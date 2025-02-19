#Formatting count data and covariates to be used as a unmarkedFrameOccu object
#Last edited: 2025-02-12

# Clear workspace
rm(list=ls()) 

# Load packages
library(tidyverse)
library(corrplot)
library(car) 

# Load the acoustic data from 2024
acoustic_dat <- read_csv("data/count_data/count_data_24_bacs.csv")
# load ARU timeline
aru_timeline <- read_csv("data/aru_timeline.csv")
# Load the vegetation data from 2024
veg <- read_csv("data/line_intercept_summary.csv")
# Load weather data 
weather <- read_csv("data/weather_station_daily_2024.csv")

# Create a new column for the site treatment grouped by the first letter of the site name
acoustic_dat <- acoustic_dat %>%
  mutate(treatment = case_when(
    str_detect(site, "^M") ~ "mine",
    str_detect(site, "^T") ~ "timber",
    str_detect(site, "^R") ~ "young_savanna",
    str_detect(site, "^O") ~ "mature_savanna",
    TRUE ~ "Unknown"
  ))

#filter dates out of acoustic_dat before 2024-03-10 - 2024-06-28
acoustic_dat <- acoustic_dat %>%
  filter(date >= '2024-03-10' & date <= '2024-06-28')

# Add a column for sp_code and make it bacs
acoustic_dat$sp_code <- "bacs"

# Convert the date column in aru_timeline to Date type
aru_timeline$date <- as.Date(aru_timeline$date, format = "%m/%d/%y")

# filter aru_timeline for dates on between 2024-03-10 and 2024-06-28
aru_timeline <- aru_timeline %>%
  filter(date >= '2024-03-10' & date <= '2024-06-28')

# Transform aru_timeline to long format
aru_long <- aru_timeline %>%
  pivot_longer(cols = -date, names_to = "site", values_to = "aru")

# If a site has an 'NA' value for an ARU on a given date, acoustic_dat should have an NA instead of 0
acoustic_dat <- acoustic_dat %>%
  left_join(aru_long, by = c("date", "site"))

# Update the 'count' to NA where 'aru' is NA
acoustic_dat$count <- ifelse(is.na(acoustic_dat$aru), NA, acoustic_dat$count)

# Drop unnecessary columns
acoustic_dat <- acoustic_dat %>% select(-aru)

# Convert 'date' column to date format in acoustic_dat
acoustic_dat$date <- as.Date(acoustic_dat$date)

# Generate a sequence of all dates present in the df
all_dates <- seq(min(acoustic_dat$date), max(acoustic_dat$date), by="day")

# Create a replicate column for each ARU and day in chronological order
acoustic_dat <- acoustic_dat %>%
  arrange(site, date) %>%
  group_by(site) %>%
  mutate(replicate = as.numeric(factor(date))) %>%
  ungroup()

y.long <- acoustic_dat %>% 
  group_by(site, date, replicate, sp_code) %>% 
  summarize(count = ifelse(all(is.na(count)), NA, max(count, na.rm = TRUE))) %>%  # Handle all NA cases
  ungroup() %>%
  glimpse()

## Beginning to set these up for loops ##

# Values to create a 3D array 
# Species codes
sp.codes <- sort(unique(y.long$sp_code))
# Plot (site) codes
plot.codes <- sort(unique(y.long$site))
# Number of species
N <- length(sp.codes)
# Maximum number of replicates at a site
K <- length(unique(y.long$replicate))
# Number of sites
J <- length(unique(y.long$site))
# Array for detection-nondetection data
y <- array(NA, dim = c(N, J, K))
# Label the dimensions of y (not necessary, but helpful)
dimnames(y)[[1]] <- sp.codes
dimnames(y)[[2]] <- plot.codes
dimnames(y)[[3]] <- 1:K

# Fill the array with detection data
for (i in 1:N) {
  for (j in 1:J) {
    site_data <- y.long %>% 
      filter(sp_code == sp.codes[i], site == plot.codes[j]) %>% 
      arrange(date, replicate)
    
    y[i, j, 1:nrow(site_data)] <- site_data$count
  }
}

# Look at the structure of the array y
str(y)

# Total number of observations for species
apply(y, 1, sum, na.rm = TRUE)

# Detection Covariates 
day.2024 <- acoustic_dat %>%
  group_by(site, replicate) %>%
  summarize(date = unique(date)) %>%
  ungroup() %>%
  glimpse()

# Loop through sites and replicates to extract Julian day
day <- matrix(NA, nrow = J, ncol = K)
for (j in 1:J) { # Loop through sites
  for (k in 1:K) { # Loop through replicate surveys
    # Get current date for each survey
    curr.vals <- day.2024 %>%
      filter(site == plot.codes[j], replicate == k) %>%
      mutate(date = yday(date)) %>%
      select(date) %>%
      arrange(date)
    
    # If the site was surveyed for the given replicate, 
    # extract the first date value. 
    if (nrow(curr.vals) > 0) {
      day[j, k] <- curr.vals$date[1]
    }
  } # k (replicates)
}  # j (sites)

# Rename rows
site_names_day <- unique(plot.codes)
rownames(day) <- site_names_day

# Convert weather dates to julian day
weather$julian_day <- yday(as.Date(weather$date))

# Remove any rows with dates before 2024-03-10 and after 2024-06-28
weather <- weather %>%
  filter(date >= "2024-03-10" & date <= "2024-06-28")

# Initialize matrices for each detection covariate with appropriate dimensions
temp <- matrix(NA, nrow = J, ncol = K)  # Temperature
wind_speed <- matrix(NA, nrow = J, ncol = K)  # Wind speed
precipitation <- matrix(NA, nrow = J, ncol = K)  # Precipitation

# Set row names for all matrices to match the 'site' column from 'weather'
site_names_weather <- unique(weather$site)
rownames(temp) <- site_names_weather
rownames(wind_speed) <- site_names_weather
rownames(precipitation) <- site_names_weather

# Loop through sites and replicates to populate the detection covariates
for (j in 1:J) {
  for (k in 1:K) {
    # Find the corresponding date for this site/replicate
    curr_date <- day[j, k]
    if (!is.na(curr_date)) {
      # Match the site and date to extract the corresponding weather data
      weather_row <- weather[weather$site == plot.codes[j] & weather$julian_day == curr_date, ]
      
      if (nrow(weather_row) > 0) {
        # Populate the matrices with the corresponding weather data
        temp[j, k] <- weather_row$avg_temp_f
        wind_speed[j, k] <- weather_row$avg_wind_speed_mph
        precipitation[j, k] <- weather_row$total_precipitation_in
      }
    }
  }
}

# Combine the detection covariates into a list
det.covs <- list(
  day = day,  # Julian day matrix
  temp = temp,  # Temperature matrix
  wind_speed = wind_speed,  # Wind speed matrix
  precipitation = precipitation  # Precipitation matrix
)

# Checking for multicollinearity between detection covariates
det.covs.df <- data.frame(
  day = as.vector(det.covs$day),
  temp = as.vector(det.covs$temp),
  wind_speed = as.vector(det.covs$wind_speed),
  precipitation = as.vector(det.covs$precipitation)
)
# Plot correlation matrix
cor_matrix_det <- cor(det.covs.df, use = "complete.obs")

# Change column and row names for clarity
colnames(cor_matrix_det) <- c("Julian Day", "Temperature Avg. (°F)", "Wind Speed Avg. (kph)", "Precipitation Avg. (in)")
rownames(cor_matrix_det) <- c("Julian Day", "Temperature Avg. (°F)", "Wind Speed Avg. (kph)", "Precipitation Avg. (in)")

# Visualize the correlation matrix with corrplot
corrplot(cor_matrix_det, method = "color", order = "alphabet", type = 'lower', 
         addCoef.col = "black", tl.col = "black", tl.srt = 45, tl.cex = 1,  number.cex = 1, cl.cex = 1)

# Select and arrange vegetation covariates
veg.covs <- veg %>% 
  select(site, treatment, shrub_cover, grass_cover, forb_cover, shrub_height, grass_height) %>%
  arrange(site)

# Match the order of sites with plot.codes
veg.covs <- veg.covs[match(plot.codes, veg.covs$site), ]

# Select only the covariate columns (without site column)
abund.covs <- veg.covs %>%
  select(treatment, shrub_cover, grass_cover, forb_cover, 
         shrub_height, grass_height)

# Convert occ.covs to a dataframe
abund.covs <- as.data.frame(abund.covs)

# Set the row names to be the site names
rownames(abund.covs) <- veg.covs$site

# View the structure to ensure it worked correctly
str(abund.covs)

# 3D y-array into 2D
y_flat <- apply(y, c(2, 3), identity)

# Check the structure 
str(y_flat)

# Create a model call to build formulas 
abund_data_bacs <- list(y = y_flat, abund.covs = abund.covs, det.covs = det.covs)

# Ensure that every NA in y has a corresponding NA in detection covariates
for(i in seq_along(abund_data_bacs$det.covs)) {
  # Assign NA where y has NA
  abund_data_bacs$det.covs[[i]][is.na(abund_data_bacs$y)] <- NA
  
  # Align row and column names with y
  rownames(abund_data_bacs$det.covs[[i]]) <- rownames(abund_data_bacs$y)
  colnames(abund_data_bacs$det.covs[[i]]) <- colnames(abund_data_bacs$y)
}

# Save as data object 
save(abund_data_bacs, file = "data/abundance_data/abundance_24_bacs.RData")

