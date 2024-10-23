# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  

# Load BACS abundance data
load("data/abundance_data/abundance_24_bacs.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_bacs$y <- ifelse(abund_data_bacs$y > 0 | is.na(abund_data_bacs$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as these RN models can't handle NAs we have in our dataset 
abund_data_bacs$y <- abund_data_bacs$y[, -c(1:55)]
abund_data_bacs$det.covs$day <- abund_data_bacs$det.covs$day[, -c(1:55)]
abund_data_bacs$det.covs$temp <- abund_data_bacs$det.covs$temp[, -c(1:55)]
abund_data_bacs$det.covs$wind_speed <- abund_data_bacs$det.covs$wind_speed[, -c(1:55)]
abund_data_bacs$det.covs$precipitation <- abund_data_bacs$det.covs$precipitation[, -c(1:55)]


abund_data_bacs$abund.covs$treatment <- factor(abund_data_bacs$abund.covs$treatment, 
                                               levels = c("mine", "rx_fire_sec_growth", 
                                                          "rx_fire_young", "timber"))


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

# Fit the models 
rn_mod_t_bacs <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
  ~ (treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + 
    scale(shrub_height) + scale(grass_height), 
  data = umf_bacs, 
  method = "BFGS"  
)

summary(rn_mod_t_bacs) 


rn_mod_bacs <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_bacs,
  K = 25,  
  method = "BFGS" 
)
summary(rn_mod_bacs)

# Fit the null model
rn_null_bacs <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
  ~ 1, 
  data = umf_bacs,
  K = 25,  
  method = "BFGS"  
)
summary(rn_null_bacs)

#AICc
AICc(rn_mod_t_bacs)
AICc(rn_mod_bacs)
AICc(rn_null_bacs)


# MB tests
gof_rn_model_t_bacs <- mb.gof.test(rn_mod_t_bacs, nsim = 500, c.hat.est = TRUE, model.type = "royle-nichols")

print(gof_rn_model_t_bacs)


## Cross Validation 
K <- 5
set.seed(123)  
folds <- sample(rep(1:K, length.out = nrow(abund_data_bacs$abund.covs)))
performance <- numeric(K)



######################### Keeping factors #########################


for (i in 1:K) {
  train_indices <- which(folds != i)
  test_indices <- which(folds == i)
  
  train_data <- umf_bacs[train_indices, ]
  test_data <- umf_bacs[test_indices, ]
  
  # Fit model to the training data 
  model <- occuRN(
    formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
    ~ factor(treatment) - 1 + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
    data = train_data, 
    K = 25,  
    method = "BFGS"
  )
  
  # Predict on the test data
  predictions <- predict(model, newdata = test_data, type = "state")
  
  # Evaluate the model using MSE
  performance[i] <- mean((test_data@y - predictions$Predicted)^2)
}

# Performance metric
mean_performance <- mean(performance)
print(mean_performance)
