# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  

# Load nobo abundance data
load("data/abundance_data/abundance_24_nobo.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_nobo$y <- ifelse(abund_data_nobo$y > 0 | is.na(abund_data_nobo$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as RN model can't handle NAs
abund_data_nobo$y <- abund_data_nobo$y[, -c(1:55)]
abund_data_nobo$det.covs$day <- abund_data_nobo$det.covs$day[, -c(1:55)]
abund_data_nobo$det.covs$temp <- abund_data_nobo$det.covs$temp[, -c(1:55)]
abund_data_nobo$det.covs$wind_speed <- abund_data_nobo$det.covs$wind_speed[, -c(1:55)]
abund_data_nobo$det.covs$precipitation <- abund_data_nobo$det.covs$precipitation[, -c(1:55)]


# Create the unmarkedFrameOccu object
umf_nobo <- unmarkedFrameOccu(
  y = as.matrix(abund_data_nobo$y),
  siteCovs = abund_data_nobo$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_nobo$det.covs$day),
    temp = as.matrix(abund_data_nobo$det.covs$temp),
    wind_speed = as.matrix(abund_data_nobo$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_nobo$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model
rn_model_t_nobo <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_nobo,
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)
summary(rn_model_t_nobo)

# Fit the second Royle-Nichols model
rn_model_nobo <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_nobo,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_nobo)

# Fit the null model
rn_model_null_nobo <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ 1, 
  data = umf_nobo,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_nobo)

# AICc comparison
models_list_nobo <- list(rn_model_t_nobo = rn_model_t_nobo, rn_model_nobo = rn_model_nobo, rn_model_null_nobo = rn_model_null_nobo)
model_names_nobo <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_nobo <- aictab(cand.set = models_list_nobo, modnames = model_names_nobo)
print(aicc_table_nobo)

# Function returning fit-statistics
fitstats_nobo <- function(rn_model_t_nobo) {
  observed <- getY(rn_model_t_nobo@data)
  expected <- fitted(rn_model_t_nobo)
  resids <- residuals(rn_model_t_nobo)
  sse <- sum(resids^2, na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
 fit_out <- c(cSSE = sse, Chisq = chisq, freemanTukey = freeTuke)
 return(fit_out)
}

#parametric bootstrapping 
# Assuming 'mod' is your fitted model
pb_nobo <- parboot(rn_model_t_nobo, fitstats_nobo, nsim=100)
print(pb_nobo)

