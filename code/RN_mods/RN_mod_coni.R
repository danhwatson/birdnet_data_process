# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  

# Load coni abundance data
load("data/abundance_data/abundance_24_coni.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_coni$y <- ifelse(abund_data_coni$y > 0 | is.na(abund_data_coni$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as RN model can't handle NAs
abund_data_coni$y <- abund_data_coni$y[, -c(1:55)]
abund_data_coni$det.covs$day <- abund_data_coni$det.covs$day[, -c(1:55)]
abund_data_coni$det.covs$temp <- abund_data_coni$det.covs$temp[, -c(1:55)]
abund_data_coni$det.covs$wind_speed <- abund_data_coni$det.covs$wind_speed[, -c(1:55)]
abund_data_coni$det.covs$precipitation <- abund_data_coni$det.covs$precipitation[, -c(1:55)]


# Create the unmarkedFrameOccu object
umf_coni <- unmarkedFrameOccu(
  y = as.matrix(abund_data_coni$y),
  siteCovs = abund_data_coni$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_coni$det.covs$day),
    temp = as.matrix(abund_data_coni$det.covs$temp),
    wind_speed = as.matrix(abund_data_coni$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_coni$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model
rn_model_t_coni <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_coni,
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)
summary(rn_model_t_coni)

# Fit the second Royle-Nichols model
rn_model_coni <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_coni,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_coni)

# Fit the null model
rn_model_null_coni <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ 1, 
  data = umf_coni,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_coni)

# AICc comparison
models_list_coni <- list(rn_model_t_coni = rn_model_t_coni, rn_model_coni = rn_model_coni, rn_model_null_coni = rn_model_null_coni)
model_names_coni <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_coni <- aictab(cand.set = models_list_coni, modnames = model_names_coni)
print(aicc_table_coni)

# Function returning fit-statistics
fitstats_coni <- function(rn_model_t_coni) {
  observed <- getY(rn_model_t_coni@data)
  expected <- fitted(rn_model_t_coni)
  resids <- residuals(rn_model_t_coni)
  sse <- sum(resids^2, na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
 fit_out <- c(cSSE = sse, Chisq = chisq, freemanTukey = freeTuke)
 return(fit_out)
}

#parametric bootstrapping 
# Assuming 'mod' is your fitted model
pb_coni <- parboot(rn_model_t_coni, fitstats_coni, nsim=100)
print(pb_coni)
