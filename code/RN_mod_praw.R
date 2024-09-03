# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  # For AICc comparison
library(cvTools)     # For cross-validation

# Load nobo abundance data
load("data/nobo_abundance_24.RData")

# Convert values greater than 0 or NA to 1 in "y"
nobo_abund_data_24$y <- ifelse(nobo_abund_data_24$y > 0 | is.na(nobo_abund_data_24$y), 1, 0)

# Create the unmarkedFrameOccu object
umf <- unmarkedFrameOccu(
  y = as.matrix(nobo_abund_data_24$y),
  siteCovs = nobo_abund_data_24$abund.covs,
  obsCovs = list(
    day = as.matrix(nobo_abund_data_24$det.covs$day),
    temp = as.matrix(nobo_abund_data_24$det.covs$temp),
    wind_speed = as.matrix(nobo_abund_data_24$det.covs$wind_speed),
    precipitation = as.matrix(nobo_abund_data_24$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model
rn_model_t <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf,
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)
summary(rn_model_t)

# Fit the second Royle-Nichols model
rn_model <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model)

# Fit the null model
rn_model_null <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ 1, 
  data = umf,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null)

# AICc comparison
models_list <- list(rn_model_t = rn_model_t, rn_model = rn_model, rn_model_null = rn_model_null)
model_names <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table <- aictab(cand.set = models_list, modnames = model_names)
print(aicc_table)



