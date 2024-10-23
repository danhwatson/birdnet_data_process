# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  
library(pROC)
library(fastDummies)

# Load bacs abundance data
load("data/abundance_data/abundance_24_bacs.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_bacs$y <- ifelse(abund_data_bacs$y > 0 | is.na(abund_data_bacs$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as these RN models can't handle NAs we have in our dataset 
abund_data_bacs$y <- abund_data_bacs$y[, -c(1:55)]
abund_data_bacs$det.covs$day <- abund_data_bacs$det.covs$day[, -c(1:55)]
abund_data_bacs$det.covs$temp <- abund_data_bacs$det.covs$temp[, -c(1:55)]
abund_data_bacs$det.covs$wind_speed <- abund_data_bacs$det.covs$wind_speed[, -c(1:55)]
abund_data_bacs$det.covs$precipitation <- abund_data_bacs$det.covs$precipitation[, -c(1:55)]

# Create dummy variables for the 'treatment' factor
abund_data_bacs$abund.covs <- dummy_cols(abund_data_bacs$abund.covs, select_columns = "treatment", remove_first_dummy = TRUE)

# Remove the original 'treatment' column
abund_data_bacs$abund.covs <- abund_data_bacs$abund.covs %>% select(-treatment)

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

# Number of folds
K <- 5

# Create folds
set.seed(123)  # For reproducibility
folds <- sample(rep(1:K, length.out = nrow(abund_data_bacs$abund.covs)))

# Initialize vectors to store performance metrics for each model
performance_t <- numeric(K)
performance_b <- numeric(K)
performance_null <- numeric(K)

for (i in 1:K) {
  # Split data into training and testing sets
  train_indices <- which(folds != i)
  test_indices <- which(folds == i)
  
  train_data <- umf_bacs[train_indices, ]
  test_data <- umf_bacs[test_indices, ]
  
  # Fit the Royle-Nichols model with treatment covariates to the training data
  rn_model_t_bacs <- occuRN(
    formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
    ~ treatment_rx_fire_sec_growth + treatment_rx_fire_young + treatment_timber + scale(shrub_cover) + scale(grass_cover),
    data = train_data, 
    K = 25,  
    method = "BFGS"
  )
  
  # Predict on the test data for rn_mod_t_bacs
  predictions_t <- predict(rn_model_t_bacs, newdata = test_data, type = "state")
  
  # Evaluate the model (e.g., using MSE) for rn_mod_t_bacs
  performance_t[i] <- mean((test_data@y - predictions_t$Predicted)^2)
  
  # Fit the Royle-Nichols model without treatment covariates to the training data
  rn_model_bacs <- occuRN(
    formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
    ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
    data = train_data, 
    K = 25,  
    method = "BFGS"
  )
  
  # Predict on the test data for rn_mod_bacs
  predictions_b <- predict(rn_model_bacs, newdata = test_data, type = "state")
  
  # Evaluate the model (e.g., using MSE) for rn_mod_bacs
  performance_b[i] <- mean((test_data@y - predictions_b$Predicted)^2)
  
  # Fit the null Royle-Nichols model to the training data
  rn_model_null_bacs <- occuRN(
    formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
    ~ 1, 
    data = train_data, 
    K = 25,  
    method = "BFGS"
  )

  # Predict on the test data for rn_null_bacs
  predictions_null <- predict(rn_null_bacs, newdata = test_data, type = "state")
  
  # Evaluate the model (e.g., using MSE) for rn_null_bacs
  performance_null[i] <- mean((test_data@y - predictions_null$Predicted)^2)
}

# Average performance metric for each model
mean_performance_t <- mean(performance_t)
mean_performance_b <- mean(performance_b)
mean_performance_null <- mean(performance_null)

print(mean_performance_t)
print(mean_performance_b)
print(mean_performance_null)
  

# MB tests
gof_rn_model_t_bacs <- mb.gof.test(rn_model_t_bacs, nsim = 100, c.hat.est = TRUE, model.type = "royle-nichols")
gof_rn_model_bacs <- mb.gof.test(rn_model_bacs, nsim = 100, c.hat.est = TRUE, model.type = "royle-nichols")
gof_rn_model_null_bacs <- mb.gof.test(rn_model_null_bacs, nsim = 100, c.hat.est = TRUE, model.type = "royle-nichols")

print(gof_rn_model_t_bacs)
print(gof_rn_model_bacs)
print(gof_rn_model_null_bacs)

