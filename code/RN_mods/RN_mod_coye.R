# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  

# Load coye abundance data
load("data/abundance_data/abundance_24_coye.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_coye$y <- ifelse(abund_data_coye$y > 0 | is.na(abund_data_coye$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as RN model can't handle NAs
abund_data_coye$y <- abund_data_coye$y[, -c(1:55)]
abund_data_coye$det.covs$day <- abund_data_coye$det.covs$day[, -c(1:55)]
abund_data_coye$det.covs$temp <- abund_data_coye$det.covs$temp[, -c(1:55)]
abund_data_coye$det.covs$wind_speed <- abund_data_coye$det.covs$wind_speed[, -c(1:55)]
abund_data_coye$det.covs$precipitation <- abund_data_coye$det.covs$precipitation[, -c(1:55)]


# Create the unmarkedFrameOccu object
umf_coye <- unmarkedFrameOccu(
  y = as.matrix(abund_data_coye$y),
  siteCovs = abund_data_coye$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_coye$det.covs$day),
    temp = as.matrix(abund_data_coye$det.covs$temp),
    wind_speed = as.matrix(abund_data_coye$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_coye$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model
rn_model_t_coye <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_coye,
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)
summary(rn_model_t_coye)

# Fit the second Royle-Nichols model
rn_model_coye <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_coye,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_coye)

# Fit the null model
rn_model_null_coye <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ 1, 
  data = umf_coye,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_coye)

# AICc comparison
models_list_coye <- list(rn_model_t_coye = rn_model_t_coye, rn_model_coye = rn_model_coye, rn_model_null_coye = rn_model_null_coye)
model_names_coye <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_coye <- aictab(cand.set = models_list_coye, modnames = model_names_coye)
print(aicc_table_coye)

# Function returning fit-statistics
fitstats_coye <- function(rn_model_t_coye) {
  observed <- getY(rn_model_t_coye@data)
  expected <- fitted(rn_model_t_coye)
  resids <- residuals(rn_model_t_coye)
  sse <- sum(resids^2, na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
  fit_out <- c(cSSE = sse, Chisq = chisq, freemanTukey = freeTuke)
  return(fit_out)
}

#parametric bootstrapping 
# Assuming 'mod' is your fitted model
pb_coye <- parboot(rn_model_t_coye, fitstats_coye, nsim=100)
print(pb_coye)