# Clear 
rm(list = ls())

# Load packages
library(tidyverse)
library(unmarked)

# Load acoustic data from July 2023
acoustic_dat <- read.csv("data/acoustic_dat_07.csv", header = TRUE)

#load veg data from 2023
veg <- read.csv("data/veg_23.csv", header = TRUE)

# Convert 'date' column to Date type in acoustic_dat
acoustic_dat$date <- as.Date(acoustic_dat$date)

# Generate a sequence of all dates present in the dataset
all_dates <- seq(min(acoustic_dat$date), max(acoustic_dat$date), by = "day")

# Create a replicate column for each ARU and day in chronological order
acoustic_dat <- acoustic_dat %>%
  arrange(aru, date) %>%
  group_by(aru) %>%
  mutate(replicate = as.numeric(factor(date))) %>%
  ungroup()

y.long <- acoustic_dat %>%
  group_by(aru, date, replicate, sp_code) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  glimpse()

# Remove ARUs not found in veg df
y.long <- y.long %>%
  filter(aru %in% veg$aru)

# Remove matches from veg that are not in y.long
veg <- veg %>%
  filter(aru %in% y.long$aru)

# Species codes
sp.codes <- sort(unique(y.long$sp_code))
# Plot (site) codes.
plot.codes <- sort(unique(y.long$aru))
# Number of species
N <- length(sp.codes)
# Maximum number of replicates at a site
K <- 31
# Number of sites
J <- length(unique(y.long$aru))
# Array for detection-nondetection data.
y <- array(NA, dim = c(N, J, K))
# Label the dimensions of y (not necessary, but helpful)
dimnames(y)[[1]] <- sp.codes
dimnames(y)[[2]] <- plot.codes
# Look at the structure of our array y
str(y)

for (j in 1:J) { # Loop through sites.
  for (k in 1:K) { # Loop through replicates at each site.
    # Extract data for current site/replicate combination.
    curr.df <- y.long %>%
      filter(aru == plot.codes[j], replicate == k)
    # Check if more than one date for a given replicate
    if (n_distinct(curr.df$date) > 1) {
      # If there is more than 1 date, only use the data
      # from the first date.
      curr.dates <- unique(sort(curr.df$date))
      curr.df <- curr.df %>%
        filter(date == curr.dates[1])
    }
    # If plot j was sampled during replicate k,
    # curr.df will have at least 1 row (i.e., at least
    # one species will be observed). If not, assume it
    # was not sampled for that replicate.
    if (nrow(curr.df) > 0) {
      # Extract the species that were observed during
      # this site/replicate.
      curr.sp <- which(sp.codes %in% curr.df$sp_code)
      # Set value to 1 for species that were observed.
      y[curr.sp, j, k] <- 1
      # Set value to 0 for all other species.
      y[-curr.sp, j, k] <- 0
    }
  } # k (replicates)
} # j (sites)
str(y)

# Total number of observations for each species
apply(y, 1, sum, na.rm = TRUE)

#Detection-nondetection matrix for species of interest
NOBO <- c('NOBO')
y_nobo <- y[which(sp.codes %in% NOBO), , ]
str(y_nobo)


occurence_covariates <- veg %>%
  select(aru, average_shrubs, average_grass, average_forb, average_bareground, average_veg_height, site_treat) %>%
  distinct()

rownames(occurence_covariates)<- occurence_covariates$aru

# Create unmarkedFrameOccu object
umf <- unmarkedFrameOccu(y = y_nobo, siteCovs = occurence_covariates)

# Fit the occupancy model including vegetation covariates
model <- occu(~1 ~ scale(average_shrubs) + scale(average_grass) + scale(average_forb), data = umf)

# Summarize the model results
summary(model)

