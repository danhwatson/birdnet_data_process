# Clear
rm(list=ls())

# Load packages
library(tidyverse)
library(unmarked)

# Load the acoustic data from July 2023
acoustic_dat <- read_csv("data/acoustic_dat_07.csv")

# Load the vegetation data from 2024
veg <- read_csv("data/veg_24.csv")

# Load weather data 
#weather <- read_csv("data/weather.csv")

# Convert 'date' column to Date type in acoustic_dat
acoustic_dat$date <- as.Date(acoustic_dat$date)

# Generate a sequence of all dates present in the dataset
all_dates <- seq(min(acoustic_dat$date), max(acoustic_dat$date), by="day")

# Create a replicate column for each ARU and day in chronological order
acoustic_dat <- acoustic_dat %>%
  arrange(aru, date) %>%
  group_by(aru) %>%
  mutate(replicate = as.numeric(factor(date))) %>%
  ungroup()

y.long <- acoustic_dat %>% 
  group_by(aru, date, replicate, sp_code) %>% 
  summarize(count=n()) %>%
  ungroup() %>%
  glimpse()

# Species codes
sp.codes <- sort(unique(y.long$sp_code))
# Plot (site) codes
plot.codes <- sort(unique(y.long$aru))
# Number of species
N <- length(sp.codes)
# Maximum number of replicates at a site
K <- 31
# Number of sites
J <- length(unique(y.long$aru))
# Array for detection-nondetection data
y <- array(NA, dim = c(N, J, K))
# Label the dimensions of y (not necessary, but helpful)
dimnames(y)[[1]] <- sp.codes
dimnames(y)[[2]] <- plot.codes

# Look at the structure of our array y
str(y)

for (j in 1:J) { # Loop through sites
  for (k in 1:K) { # Loop through replicates at each site
    # Extract data for current site/replicate combination
    curr.df <- y.long %>%
      filter(aru == plot.codes[j], replicate == k)
    # Check if more than one date for a given replicate
    if (n_distinct(curr.df$date) > 1) {
      # If there is more than 1 date, only use the data from the first date
      curr.dates <- unique(sort(curr.df$date))
      curr.df <- curr.df %>% 
        filter(date == curr.dates[1])
    }
    # If plot j was sampled during replicate k, 
    # curr.df will have at least 1 row (i.e., at least 
    # one species will be observed). If not, assume it 
    # was not sampled for that replicate.
    if (nrow(curr.df) > 0) {
      # Extract the species that were observed during this site/replicate
      curr.sp <- which(sp.codes %in% curr.df$sp_code)
      # Set value to 1 for species that were observed
      y[curr.sp, j, k] <- 1
      # Set value to 0 for all other species
      y[-curr.sp, j, k] <- 0
    }
  } # k (replicates)
} # j (sites)

# Look at the structure of the array y
str(y)

# Total number of observations for each species
apply(y, 1, sum, na.rm = TRUE)

# Detection-nondetection matrix for species of interest
NOBO <- c('NOBO')
y_nobo <- y[which(sp.codes %in% NOBO), , ]
str(y_nobo)

# Manually specify the ARU names you want to keep
selected_arus <- c('SMA10463', 'SMA10416', 'SMA10460', 'SMA10448', 'SMA10426',
                   'SMA10470', 'SMA10449', 'SMA10424', 'SMA10458', 'SMA10451',
                   'SMA10465', 'SMA10428', 'SMA10469', 'SMA10471', 'SMA10417', 'SMA10395')

# Assuming ARU names are in the row names of y_nobo, filter the data
filtered_y_nobo <- y_nobo[rownames(y_nobo) %in% selected_arus, , drop = FALSE]

# Print structure of filtered data
str(filtered_y_nobo)

# Filter vegetation data for selected ARUs for siteCovs
site_occurence_covariates <- veg %>%
  filter(aru_23 %in% selected_arus) %>%
  select(aru_23, avg_cover_conifer, avg_cover_grass, avg_cover_shrub, 
         avg_cover_forb, avg_height_conifer, avg_height_grass,
         avg_height_shrub,avg_height_forb, overall_avg_cover,overall_avg_height, site_treat)

# Remove rows containing "SMA10471" and "SMA10416" from site_covariates 
filtered_site_occurence_covariates <- site_occurence_covariates %>%
  filter(aru_23 != "SMA10471" & aru_23 != "SMA10416")

rownames(filtered_site_occurence_covariates) <- filtered_site_occurence_covariates$aru_23

# Create the unmarkedFrame object
umf <- unmarkedFrameOccu(y = filtered_y_nobo, siteCovs = filtered_site_occurence_covariates)

# Fit the occupancy model
model <- occu(~1 ~ scale(avg_cover_shrub) + scale(avg_cover_grass) + scale(avg_cover_forb), data = umf)

# Print summary of the model
summary(model)
