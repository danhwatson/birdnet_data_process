# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  

# Load wevi abundance data
load("data/abundance_data/abundance_24_wevi.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_wevi$y <- ifelse(abund_data_wevi$y > 0 | is.na(abund_data_wevi$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as RN model can't handle NAs
abund_data_wevi$y <- abund_data_wevi$y[, -c(1:55)]
abund_data_wevi$det.covs$day <- abund_data_wevi$det.covs$day[, -c(1:55)]
abund_data_wevi$det.covs$temp <- abund_data_wevi$det.covs$temp[, -c(1:55)]
abund_data_wevi$det.covs$wind_speed <- abund_data_wevi$det.covs$wind_speed[, -c(1:55)]
abund_data_wevi$det.covs$precipitation <- abund_data_wevi$det.covs$precipitation[, -c(1:55)]


# Create the unmarkedFrameOccu object
umf_wevi <- unmarkedFrameOccu(
  y = as.matrix(abund_data_wevi$y),
  siteCovs = abund_data_wevi$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_wevi$det.covs$day),
    temp = as.matrix(abund_data_wevi$det.covs$temp),
    wind_speed = as.matrix(abund_data_wevi$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_wevi$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model
rn_model_t_wevi <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_wevi,
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)
summary(rn_model_t_wevi)

# Fit the second Royle-Nichols model
rn_model_wevi <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_wevi,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_wevi)

# Fit the null model
rn_model_null_wevi <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ 1, 
  data = umf_wevi,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_wevi)

# AICc comparison
models_list_wevi <- list(rn_model_t_wevi = rn_model_t_wevi, rn_model_wevi = rn_model_wevi, rn_model_null_wevi = rn_model_null_wevi)
model_names_wevi <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_wevi <- aictab(cand.set = models_list_wevi, modnames = model_names_wevi)
print(aicc_table_wevi)

# Function returning fit-statistics
fitstats_wevi <- function(rn_model_t_wevi) {
  observed <- getY(rn_model_t_wevi@data)
  expected <- fitted(rn_model_t_wevi)
  resids <- residuals(rn_model_t_wevi)
  sse <- sum(resids^2, na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
  fit_out <- c(cSSE = sse, Chisq = chisq, freemanTukey = freeTuke)
  return(fit_out)
}

#parametric bootstrapping 
# Assuming 'mod' is your fitted model
pb_wevi <- parboot(rn_model_t_wevi, fitstats_wevi, nsim=100)
print(pb_wevi)