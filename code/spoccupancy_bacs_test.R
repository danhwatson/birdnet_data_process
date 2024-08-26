# Clear
rm(list=ls())

# Load packages
library(tidyverse)
library(spOccupancy)
library(corrplot)
library(car) 
library(pander)
library(knitr)

### Lots of housekeeping needed ###
# Load the acoustic data from July 2023
acoustic_dat <- read_csv("data/bacs_count_data_24.csv")
# load ARU timeline
aru_timeline <- read_csv("data/aru_timeline.csv")
# Load the vegetation data from 2024
veg <- read_csv("data/line_intercept_summary.csv")
# Load weather data 
weather <- read_csv("data/weather_daily.csv")
# Load coordinates for each site
coords <- read_csv("data/aru_site_coordinates.csv")

## Check for multicollinearity in the veg data ## 
cor_matrix_veg <- cor(veg %>% select(shrub_cover, forb_cover, grass_cover, conifer_cover, shrub_height, forb_height, grass_height, conifer_height), use = "complete.obs")

corrplot(cor_matrix_veg, method = "color", addCoef.col = "black", tl.col = "black", tl.srt = 45)

veg_treatment_cor_mod <- lm(shrub_cover ~ treatment + forb_cover + grass_cover + conifer_cover + 
                            shrub_height + forb_height + grass_height + conifer_height, 
                            data = veg)

# Calculate VIF
vif_treatment <- vif(veg_treatment_cor_mod)

# Create a data frame for VIF results, including row names as a 'Variable' column
vif_table <- data.frame(
  Variable = rownames(vif_treatment),
  GVIF = vif_treatment[, "GVIF"],
  Df = vif_treatment[, "Df"],
  `GVIF^(1/(2*Df))` = vif_treatment[, "GVIF^(1/(2*Df))"],
  stringsAsFactors = FALSE
)

# Remove row names from the data frame to prevent pander from adding them as a separate column
rownames(vif_table) <- NULL

# Use pander to create a clean table without the extra &nbsp; column
pander(
  vif_table, 
  caption = "VIF Results for Veg Treatment Model", 
  row.names = FALSE
)


## Some detection data wrangling ##
# Change counts greater than 0 to 1 for detection / non-detection data
acoustic_dat$count[acoustic_dat$count > 0] <- 1

# Add a column for sp_code and make it BACS
acoustic_dat$sp_code <- "BACS"

# Remove all 2023 dates from date column 
aru_timeline <- aru_timeline %>%
  filter(!grepl("2023", date)) %>%
  filter(date >= '2024-03-10' & date <= '2024-06-28')

# Convert the date column in aru_timeline to Date type
aru_timeline$date <- as.Date(aru_timeline$date, format = "%m/%d/%y")

# Transform aru_timeline to long format
aru_long <- aru_timeline %>%
  pivot_longer(cols = -date, names_to = "site", values_to = "aru")

# IF a site has an 'NA' value for an ARU on a given date, acoustic_dat should have an NA instead of 0
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

# write_csv(y.long, "data/bacs_detection_data.csv")
### You can just read this in after the above formatting ###
# y.long <- read_csv("data/bacs_detection_data.csv")

## Beginning to set these up for loops ##

# Values to create a 3D array for spOccupancy
# Species codes
sp.codes <- sort(unique(y.long$sp_code))
# Plot (site) codes
plot.codes <- sort(unique(y.long$site))
# Number of species
N <- length(sp.codes)
# Maximum number of replicates at a site
K <- 111
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
corrplot(cor_matrix_det, method = "color", addCoef.col = "black", tl.col = "black", tl.srt = 45)

# Select and arrange vegetation covariates
veg.covs <- veg %>% 
  select(site, treatment, shrub_cover, grass_cover, forb_cover) %>%
  arrange(site)

# Match the order of sites with plot.codes
veg.covs <- veg.covs[match(plot.codes, veg.covs$site), ]

# Select only the covariate columns (without site column)
occ.covs <- veg.covs %>%
  select(treatment, shrub_cover, grass_cover, forb_cover)

# Convert occ.covs to a dataframe
occ.covs <- as.data.frame(occ.covs)

# Set the row names to be the site names
rownames(occ.covs) <- veg.covs$site

# View the structure to ensure it worked correctly
str(occ.covs)

# Format UTM coordinate data for use in the spatial models 
coords <- coords %>%
  select(site, X, Y) %>%
  arrange(site)

# Convert to a matrix and assign row names as site names
coords_matrix <- as.matrix(coords[, c("X", "Y")])
rownames(coords_matrix) <- coords$site

# Check coordinates are right
plot(bacs_data$coords, pch = 19) 

# Not a MSOM, so turn 3D y-array into 2D
y_flat <- apply(y, c(2, 3), identity)

# Check the structure 
str(y_flat)

# Create a model call to build formulas 
bacs_data <- list(y = y_flat, occ.covs = occ.covs, det.covs = det.covs, coords = coords_matrix)

# Save as data object 
# save(bacs_data, file = "data/bacs_data.RData")

# Load data object
load("data/bacs_data.RData")

# Create covariate formulas for clearer argument building 

# Null
null.occ.formula <- ~ 1
null.det.formula <- ~ scale(day) + scale(wind_speed) + scale(precipitation)

# Without treatment
occ.formula <- ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover)
det.formula <- ~ scale(day) + scale(wind_speed) + scale(precipitation)

# With treatment
t_occ.formula <- ~ (treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover)
t_det.formula <- ~ scale(day) +  scale(wind_speed) + scale(precipitation)

# Null model
null_out <- PGOcc(
  occ.formula = ~1,
  det.formula = det.formula,
  data = bacs_data,
  inits =,
  n.samples = 10000,
  priors =,
  n.omp.threads = 1,
  verbose = TRUE,
  n.report = 1000,
  n.burn = 1000,
  n.thin = 1,
  n.chains = 5
)
summary(null_out)

# Null model with NNGP
sp_null_out <- spPGOcc(
  occ.formula = ~1,
  det.formula = det.formula,
  data = bacs_data,
  inits =,
  priors =,
  n.batch = 400,
  batch.length = 25,
  NNGP = TRUE,
  n.neighbors = 5,
  n.burn = 5000,
  n.thin = 10,
  n.chains = 5
)
summary(sp_null_out)

# Run the non-spatial model
out <- PGOcc(
  occ.formula = occ.formula,
  det.formula = det.formula,
  data = bacs_data,
  inits =,
  n.samples = 10000,
  priors =,
  n.omp.threads = 1,
  verbose = TRUE,
  n.report = 1000,
  n.burn = 1000,
  n.thin = 1,
  n.chains = 5
)
summary(out)


# Run the spatial model
sp_out <- spPGOcc(
  occ.formula = occ.formula,
  det.formula = det.formula,
  data = bacs_data,
  inits =,
  priors =,
  n.batch = 400,
  batch.length = 25,
  NNGP = TRUE,
  n.neighbors = 5,
  n.burn = 5000,
  n.thin = 10,
  n.chains = 5,
  n.report = 100,
  verbose = TRUE
)
summary(sp_out)
plot(sp_out, 'beta', density = FALSE) # Good way to visualize convergence 
plot(sp_out, 'alpha', density = FALSE)

# Run the non-spatial model with treatment
t_out <- PGOcc(
  occ.formula = t_occ.formula,
  det.formula = t_det.formula,
  data = bacs_data,
  inits =,
  n.samples = 10000,
  priors =,
  n.omp.threads = 1,
  verbose = TRUE,
  n.report = 1000,
  n.burn = 1000,
  n.thin = 1,
  n.chains = 5
)

summary(t_out)

# Run the spatial model with treatment
t_sp_out <- spPGOcc(
  occ.formula = t_occ.formula,
  det.formula = t_det.formula,
  data = bacs_data,
  inits =,
  priors =,
  n.batch = 400,
  batch.length = 25,
  NNGP = TRUE,
  n.neighbors = 5,
  n.burn = 5000,
  n.thin = 10,
  n.chains = 5,
  n.report = 100,
  verbose = TRUE
)
summary(t_sp_out)

# Calculate WAIC
null_out_waic <- waicOcc(null_out)
sp_null_out_waic <- waicOcc(sp_null_out)
out_waic0cc <- waicOcc(out)
sp_out_waic <- waicOcc(sp_out) 
t_out_waic <- waicOcc(t_out)
t_sp_out_waic <- waicOcc(t_sp_out)


# Create a data frame of WAIC values
waic_values <- data.frame(
  Model = c("null_out", "sp_null_out", "out", "sp_out", "t_out", "t_sp_out"),
  elpd = c(null_out_waic["elpd"],
           sp_null_out_waic["elpd"],
           out_waic0cc["elpd"],
           sp_out_waic["elpd"],
           t_out_waic["elpd"],
           t_sp_out_waic["elpd"]),
  pD = c(null_out_waic["pD"],
         sp_null_out_waic["pD"],
         out_waic0cc["pD"],
         sp_out_waic["pD"],
         t_out_waic["pD"],
         t_sp_out_waic["pD"]),
  WAIC = c(null_out_waic["WAIC"],
           sp_null_out_waic["WAIC"],
           out_waic0cc["WAIC"],
           sp_out_waic["WAIC"],
           t_out_waic["WAIC"],
           t_sp_out_waic["WAIC"])
)

# WAIC table 
knitr::kable(waic_values, caption = "WAIC Results for Models")

# Posterior predictive checks for the best model (sp_out)
ppc_sp_out_group_1 <- ppcOcc(object = sp_out, fit.stat = "freeman-tukey", group = 1)

ppc_sp_out_group_2 <- ppcOcc(object = sp_out, fit.stat = "freeman-tukey", group = 2)
summary(ppc_sp_out_group_1)
summary(ppc_sp_out_group_2)


ppc_df_group_1 <- data.frame(fit = ppc_sp_out_group_1$fit., 
                     fit.rep = ppc_sp_out_group_1$fit.y.rep, 
                     color = 'lightskyblue1')
ppc_df_group_1$color[ppc_df_group_1$fit.rep > ppc_df_group_1$fit] <- 'lightsalmon'
plot(ppc_df_group_1$fit, ppc_df_group_1$fit.rep, bg = ppc_df_group_1$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc_df_group_1$fit, ppc_df_group_1$fit, col = 'black')

ppc_df_group_2 <- data.frame(fit = ppc_sp_out_group_2$fit., 
                     fit.rep = ppc_sp_out_group_2$fit.y.rep, 
                     color = 'lightskyblue1')
ppc_df_group_2$color[ppc_df_group_2$fit.rep > ppc_df_group_2$fit] <- 'lightsalmon'
plot(ppc_df_group_2$fit, ppc_df_group_2$fit.rep, bg = ppc_df_group_2$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc_df_group_2$fit, ppc_df_group_2$fit, col = 'black')



diff.fit <- ppc_sp_out_group_1$fit.y.rep.group.quants[3, ] - ppc_sp_out_group_1$fit.y.group.quants[3, ] 
plot(diff.fit, pch = 19, xlab = 'Sites (group 1 PPC)', ylab = 'Replicate - True Discrepancy')

diff.fit <- ppc_sp_out_group_2$fit.y.rep.group.quants[3, ] - ppc_sp_out_group_2$fit.y.group.quants[3, ]
plot(diff.fit, pch = 19, xlab = 'Replicate Surveys (group 2 PPC)', ylab = 'Replicate - True Discrepancy')


ppc_df_group_1 <- data.frame(fit = ppc_sp_out_group_1$fit., 
                             fit.rep = ppc_sp_out_group_1$fit.y.rep, 
                             color = 'lightskyblue1')
ppc_df_group_1$color[ppc_df_group_1$fit.rep > ppc_df_group_1$fit] <- 'lightsalmon'
plot(ppc_df_group_1$fit, ppc_df_group_1$fit.rep, bg = ppc_df_group_1$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc_df_group_1$fit, ppc_df_group_1$fit, col = 'black')

ppc_df_group_2 <- data.frame(fit = ppc_sp_out_group_2$fit., 
                             fit.rep = ppc_sp_out_group_2$fit.y.rep, 
                             color = 'lightskyblue1')
ppc_df_group_2$color[ppc_df_group_2$fit.rep > ppc_df_group_2$fit] <- 'lightsalmon'
plot(ppc_df_group_2$fit, ppc_df_group_2$fit.rep, bg = ppc_df_group_2$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc_df_group_2$fit, ppc_df_group_2$fit, col = 'black')





