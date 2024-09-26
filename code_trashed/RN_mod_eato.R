# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  

# Load eato abundance data
load("data/abundance_data/abundance_24_eato.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_eato$y <- ifelse(abund_data_eato$y > 0 | is.na(abund_data_eato$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as RN model can't handle NAs
abund_data_eato$y <- abund_data_eato$y[, -c(1:55)]
abund_data_eato$det.covs$day <- abund_data_eato$det.covs$day[, -c(1:55)]
abund_data_eato$det.covs$temp <- abund_data_eato$det.covs$temp[, -c(1:55)]
abund_data_eato$det.covs$wind_speed <- abund_data_eato$det.covs$wind_speed[, -c(1:55)]
abund_data_eato$det.covs$precipitation <- abund_data_eato$det.covs$precipitation[, -c(1:55)]


# Create the unmarkedFrameOccu object
umf_eato <- unmarkedFrameOccu(
  y = as.matrix(abund_data_eato$y),
  siteCovs = abund_data_eato$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_eato$det.covs$day),
    temp = as.matrix(abund_data_eato$det.covs$temp),
    wind_speed = as.matrix(abund_data_eato$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_eato$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model
rn_model_t_eato <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_eato,
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)
summary(rn_model_t_eato)

# Fit the second Royle-Nichols model
rn_model_eato <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_eato,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_eato)

# Fit the null model
rn_model_null_eato <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ 1, 
  data = umf_eato,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_eato)

# AICc comparison
models_list_eato <- list(rn_model_t_eato = rn_model_t_eato, rn_model_eato = rn_model_eato, rn_model_null_eato = rn_model_null_eato)
model_names_eato <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_eato <- aictab(cand.set = models_list_eato, modnames = model_names_eato)
print(aicc_table_eato)

# Function returning fit-statistics
fitstats_eato <- function(rn_model_t_eato) {
  observed <- getY(rn_model_t_eato@data)
  expected <- fitted(rn_model_t_eato)
  resids <- residuals(rn_model_t_eato)
  sse <- sum(resids^2, na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
  fit_out <- c(cSSE = sse, Chisq = chisq, freemanTukey = freeTuke)
  return(fit_out)
}

#parametric bootstrapping 
# Assuming 'mod' is your fitted model
pb_eato <- parboot(rn_model_t_eato, fitstats_eato, nsim=100)
print(pb_eato)