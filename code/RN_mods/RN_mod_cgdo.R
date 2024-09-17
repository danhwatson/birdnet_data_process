# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  

# Load cgdo abundance data
load("data/abundance_data/abundance_24_cgdo.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_cgdo$y <- ifelse(abund_data_cgdo$y > 0 | is.na(abund_data_cgdo$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as RN model can't handle NAs
abund_data_cgdo$y <- abund_data_cgdo$y[, -c(1:55)]
abund_data_cgdo$det.covs$day <- abund_data_cgdo$det.covs$day[, -c(1:55)]
abund_data_cgdo$det.covs$temp <- abund_data_cgdo$det.covs$temp[, -c(1:55)]
abund_data_cgdo$det.covs$wind_speed <- abund_data_cgdo$det.covs$wind_speed[, -c(1:55)]
abund_data_cgdo$det.covs$precipitation <- abund_data_cgdo$det.covs$precipitation[, -c(1:55)]


# Create the unmarkedFrameOccu object
umf_cgdo <- unmarkedFrameOccu(
  y = as.matrix(abund_data_cgdo$y),
  siteCovs = abund_data_cgdo$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_cgdo$det.covs$day),
    temp = as.matrix(abund_data_cgdo$det.covs$temp),
    wind_speed = as.matrix(abund_data_cgdo$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_cgdo$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model
rn_model_t_cgdo <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_cgdo,
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)
summary(rn_model_t_cgdo)

# Fit the second Royle-Nichols model
rn_model_cgdo <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_cgdo,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_cgdo)

# Fit the null model
rn_model_null_cgdo <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ 1, 
  data = umf_cgdo,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_cgdo)

# AICc comparison
models_list_cgdo <- list(rn_model_t_cgdo = rn_model_t_cgdo, rn_model_cgdo = rn_model_cgdo, rn_model_null_cgdo = rn_model_null_cgdo)
model_names_cgdo <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_cgdo <- aictab(cand.set = models_list_cgdo, modnames = model_names_cgdo)
print(aicc_table_cgdo)

# Function returning fit-statistics
fitstats_cgdo <- function(rn_model_t_cgdo) {
  observed <- getY(rn_model_t_cgdo@data)
  expected <- fitted(rn_model_t_cgdo)
  resids <- residuals(rn_model_t_cgdo)
  sse <- sum(resids^2, na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
  fit_out <- c(cSSE = sse, Chisq = chisq, freemanTukey = freeTuke)
  return(fit_out)
}

#parametric bootstrapping 
# Assuming 'mod' is your fitted model
pb_cgdo <- parboot(rn_model_t_cgdo, fitstats_cgdo, nsim=100)
print(pb_cgdo)