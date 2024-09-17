# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  

# Load praw abundance data
load("data/abundance_data/abundance_24_praw.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_praw$y <- ifelse(abund_data_praw$y > 0 | is.na(abund_data_praw$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as RN model can't handle NAs
abund_data_praw$y <- abund_data_praw$y[, -c(1:55)]
abund_data_praw$det.covs$day <- abund_data_praw$det.covs$day[, -c(1:55)]
abund_data_praw$det.covs$temp <- abund_data_praw$det.covs$temp[, -c(1:55)]
abund_data_praw$det.covs$wind_speed <- abund_data_praw$det.covs$wind_speed[, -c(1:55)]
abund_data_praw$det.covs$precipitation <- abund_data_praw$det.covs$precipitation[, -c(1:55)]


# Create the unmarkedFrameOccu object
umf_praw <- unmarkedFrameOccu(
  y = as.matrix(abund_data_praw$y),
  siteCovs = abund_data_praw$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_praw$det.covs$day),
    temp = as.matrix(abund_data_praw$det.covs$temp),
    wind_speed = as.matrix(abund_data_praw$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_praw$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model
rn_model_t_praw <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_praw,
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)
summary(rn_model_t_praw)

# Fit the second Royle-Nichols model
rn_model_praw <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_praw,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_praw)

# Fit the null model
rn_model_null_praw <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ 1, 
  data = umf_praw,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_praw)

# AICc comparison
models_list_praw <- list(rn_model_t_praw = rn_model_t_praw, rn_model_praw = rn_model_praw, rn_model_null_praw = rn_model_null_praw)
model_names_praw <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_praw <- aictab(cand.set = models_list_praw, modnames = model_names_praw)
print(aicc_table_praw)

# Function returning fit-statistics
fitstats_praw <- function(rn_model_t_praw) {
  observed <- getY(rn_model_t_praw@data)
  expected <- fitted(rn_model_t_praw)
  resids <- residuals(rn_model_t_praw)
  sse <- sum(resids^2, na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
  fit_out <- c(cSSE = sse, Chisq = chisq, freemanTukey = freeTuke)
  return(fit_out)
}

#parametric bootstrapping 
# Assuming 'mod' is your fitted model
pb_praw <- parboot(rn_model_t_praw, fitstats_praw, nsim=100)
print(pb_praw)

