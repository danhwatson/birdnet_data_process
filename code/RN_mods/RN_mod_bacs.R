# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  

# Load bacs abundance data
load("data/abundance_data/abundance_24_bacs.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_bacs$y <- ifelse(abund_data_bacs$y > 0 | is.na(abund_data_bacs$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as RN model can't handle NAs
abund_data_bacs$y <- abund_data_bacs$y[, -c(1:55)]
abund_data_bacs$det.covs$day <- abund_data_bacs$det.covs$day[, -c(1:55)]
abund_data_bacs$det.covs$temp <- abund_data_bacs$det.covs$temp[, -c(1:55)]
abund_data_bacs$det.covs$wind_speed <- abund_data_bacs$det.covs$wind_speed[, -c(1:55)]
abund_data_bacs$det.covs$precipitation <- abund_data_bacs$det.covs$precipitation[, -c(1:55)]


# Create the unmarkedFrameOccu object
umf_bacs <- unmarkedFrameOccu(
  y = as.matrix(abund_data_bacs$y),
  siteCovs = abund_data_bacs$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_bacs$det.covs$day),
    temp = as.matrix(abund_data_bacs$det.covs$temp),
    wind_speed = as.matrix(abund_data_bacs$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_bacs$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model
rn_model_t_bacs <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_bacs,
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)
summary(rn_model_t_bacs)

# Fit the second Royle-Nichols model
rn_model_bacs <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_bacs,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_bacs)

# Fit the null model
rn_model_null_bacs <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ 1, 
  data = umf_bacs,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_bacs)

# AICc comparison
models_list_bacs <- list(rn_model_t_bacs = rn_model_t_bacs, rn_model_bacs = rn_model_bacs, rn_model_null_bacs = rn_model_null_bacs)
model_names_bacs <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_bacs <- aictab(cand.set = models_list_bacs, modnames = model_names_bacs)
print(aicc_table_bacs)

# Function returning fit-statistics
fitstats_bacs <- function(rn_model_t_bacs) {
  observed <- getY(rn_model_t_bacs@data)
  expected <- fitted(rn_model_t_bacs)
  resids <- residuals(rn_model_t_bacs)
  sse <- sum(resids^2, na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
  fit_out <- c(cSSE = sse, Chisq = chisq, freemanTukey = freeTuke)
  return(fit_out)
}

#parametric bootstrapping 
# Assuming 'mod' is your fitted model
pb_bacs <- parboot(rn_model_t_bacs, fitstats_bacs, nsim=100)
print(pb_bacs)