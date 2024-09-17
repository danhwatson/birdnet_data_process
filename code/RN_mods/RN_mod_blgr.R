# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  

# Load blgr abundance data
load("data/abundance_data/abundance_24_blgr.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_blgr$y <- ifelse(abund_data_blgr$y > 0 | is.na(abund_data_blgr$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as RN model can't handle NAs
abund_data_blgr$y <- abund_data_blgr$y[, -c(1:55)]
abund_data_blgr$det.covs$day <- abund_data_blgr$det.covs$day[, -c(1:55)]
abund_data_blgr$det.covs$temp <- abund_data_blgr$det.covs$temp[, -c(1:55)]
abund_data_blgr$det.covs$wind_speed <- abund_data_blgr$det.covs$wind_speed[, -c(1:55)]
abund_data_blgr$det.covs$precipitation <- abund_data_blgr$det.covs$precipitation[, -c(1:55)]


# Create the unmarkedFrameOccu object
umf_blgr <- unmarkedFrameOccu(
  y = as.matrix(abund_data_blgr$y),
  siteCovs = abund_data_blgr$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_blgr$det.covs$day),
    temp = as.matrix(abund_data_blgr$det.covs$temp),
    wind_speed = as.matrix(abund_data_blgr$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_blgr$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model
rn_model_t_blgr <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_blgr,
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)
summary(rn_model_t_blgr)

# Fit the second Royle-Nichols model
rn_model_blgr <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_blgr,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_blgr)

# Fit the null model
rn_model_null_blgr <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ 1, 
  data = umf_blgr,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_blgr)

# AICc comparison
models_list_blgr <- list(rn_model_t_blgr = rn_model_t_blgr, rn_model_blgr = rn_model_blgr, rn_model_null_blgr = rn_model_null_blgr)
model_names_blgr <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_blgr <- aictab(cand.set = models_list_blgr, modnames = model_names_blgr)
print(aicc_table_blgr)

# Function returning fit-statistics
fitstats_blgr <- function(rn_model_t_blgr) {
  observed <- getY(rn_model_t_blgr@data)
  expected <- fitted(rn_model_t_blgr)
  resids <- residuals(rn_model_t_blgr)
  sse <- sum(resids^2, na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
  fit_out <- c(cSSE = sse, Chisq = chisq, freemanTukey = freeTuke)
  return(fit_out)
}

#parametric bootstrapping 
# Assuming 'mod' is your fitted model
pb_blgr <- parboot(rn_model_t_blgr, fitstats_blgr, nsim=100)
print(pb_blgr)