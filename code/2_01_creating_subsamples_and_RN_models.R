#Testing how sub-sampling data for Blue Grosbeak and Bachman's Sparrow affects RN model estimates
#Sub-sampling a full data set
#Fitting RN models to them
#Comparing AICc values
#Checking for multi-collinearity with VIF 
#Creating predicted relatives abundance estimates for each treatment
#Last updated 2025/02/12

# Clear workspace
rm(list = ls()) 

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)

#How rows correspond to dates, just to orient sub-sampling
date_range <- read.csv("data/date_range_for_subsampling.csv")
#56 = 5/3/24, 105 = 6/21/24

# Load blgr abundance data
load("data/abundance_data/abundance_24_blgr.RData")

# Restrict to columns 56 to 105
abund_data_blgr$y <- abund_data_blgr$y[, 56:105]
abund_data_blgr$det.covs$day <- abund_data_blgr$det.covs$day[, 56:105]
abund_data_blgr$det.covs$temp <- abund_data_blgr$det.covs$temp[, 56:105]
abund_data_blgr$det.covs$wind_speed <- abund_data_blgr$det.covs$wind_speed[, 56:105]
abund_data_blgr$det.covs$precipitation <- abund_data_blgr$det.covs$precipitation[, 56:105]

# Convert values greater than 0 or NA to 1 in "y"
abund_data_blgr$y <- ifelse(abund_data_blgr$y > 0 | is.na(abund_data_blgr$y), 1, 0)

# Extract detection covariates
det.covs <- list(
  day = abund_data_blgr$det.covs$day,
  temp = abund_data_blgr$det.covs$temp,
  wind_speed = abund_data_blgr$det.covs$wind_speed,
  precipitation = abund_data_blgr$det.covs$precipitation
)

## Get the actual starting column index
start_col_label <- "56"
start_col <- which(colnames(abund_data_blgr$y) == start_col_label)

# Retain only every nth day
retain_every_nth <- function(data, start_col, n) {
  cols_to_keep <- seq(start_col, ncol(data), by = n)
  return(data[, cols_to_keep])
}

# Define the intervals
n_values <- c(2, 3, 4, 7, 10)

# Create datasets for each interval
datasets <- list()

for (n in n_values) {
  # Retain only the desired columns for y and det.covs
  reduced_y <- retain_every_nth(abund_data_blgr$y, start_col, n)
  reduced_det.covs <- lapply(det.covs, function(cov) retain_every_nth(cov, start_col, n))
  
  # Store the reduced dataset in the list
  datasets[[paste0("every_", n, "_days","_blgr")]] <- list(
    y = reduced_y,
    det.covs = reduced_det.covs,
    abund.covs = abund_data_blgr$abund.covs,
    coords = abund_data_blgr$coords
  )
}

# Save each dataset in the list as an RData file with the correct name
for (name in names(datasets)) {
  # Extract the dataset from the list
  assign(name, datasets[[name]])  # Assign dataset to a variable with the desired name
  
  # Save the dataset with its unique name
  save(list = name, file = paste0("data/abundance_data/blgr_subsamples/", name, ".RData"))
}

# Load RData files for the new intervals
load("data/abundance_data/blgr_subsamples/every_2_days_blgr.RData")
load("data/abundance_data/blgr_subsamples/every_3_days_blgr.RData")
load("data/abundance_data/blgr_subsamples/every_4_days_blgr.RData")
load("data/abundance_data/blgr_subsamples/every_7_days_blgr.RData")
load("data/abundance_data/blgr_subsamples/every_10_days_blgr.RData")

#Rename to every_day_blgr
every_day_blgr <- abund_data_blgr

#umf for every_day
umf_every_day_blgr <- unmarkedFrameOccu(
  y = as.matrix(every_day_blgr$y),
  siteCovs = every_day_blgr$abund.covs,
  obsCovs = list(
    day = as.matrix(every_day_blgr$det.covs$day),
    temp = as.matrix(every_day_blgr$det.covs$temp),
    wind_speed = as.matrix(every_day_blgr$det.covs$wind_speed),
    precipitation = as.matrix(every_day_blgr$det.covs$precipitation)
  )
)

#Global model
every_day_blgr <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_every_day_blgr, 
  K = 25,  
  method = "BFGS"  
)

#Summary
summary(every_day_blgr)

#Just treatment
every_day_t_blgr <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ factor(treatment) -1,
  data = umf_every_day_blgr, 
  K = 25,  
  method = "BFGS"  
)

#Just vegetation
every_day_v_blgr <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_every_day_blgr, 
  K = 25,  
  method = "BFGS"  
)

#Null Model
every_day_null_blgr <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ 1,
  data = umf_every_day_blgr, 
  K = 25,  
  method = "BFGS"  
)

########AICc comparison and AICc table formatting###########

models_list_blgr <- list(every_day_blgr = every_day_blgr, 
                         every_day_t_blgr = every_day_t_blgr, 
                         every_day_v_blgr = every_day_v_blgr, 
                         every_day_null_blgr = every_day_null_blgr)

model_names_blgr <- c("TV", "T", "V", "Null")
aicc_table_blgr <- aictab(cand.set = models_list_blgr, modnames = model_names_blgr)
print(aicc_table_blgr)


# Every 2 days
umf_every_2_days_blgr <- unmarkedFrameOccu(
  y = as.matrix(every_2_days_blgr$y),
  siteCovs = every_2_days_blgr$abund.covs,
  obsCovs = list(
    day = as.matrix(every_2_days_blgr$det.covs$day),
    temp = as.matrix(every_2_days_blgr$det.covs$temp),
    wind_speed = as.matrix(every_2_days_blgr$det.covs$wind_speed),
    precipitation = as.matrix(every_2_days_blgr$det.covs$precipitation)
  )
)

every_2_days_blgr <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ factor(treatment) - 1 + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_every_2_days_blgr, 
  K = 25,  
  method = "BFGS"  
)

# Every 3 days
umf_every_3_days_blgr <- unmarkedFrameOccu(
  y = as.matrix(every_3_days_blgr$y),
  siteCovs = every_3_days_blgr$abund.covs,
  obsCovs = list(
    day = as.matrix(every_3_days_blgr$det.covs$day),
    temp = as.matrix(every_3_days_blgr$det.covs$temp),
    wind_speed = as.matrix(every_3_days_blgr$det.covs$wind_speed),
    precipitation = as.matrix(every_3_days_blgr$det.covs$precipitation)
  )
)

every_3_days_blgr <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ factor(treatment) - 1 + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_every_3_days_blgr, 
  K = 25,  
  method = "BFGS"  
)

# Every 4 days
umf_every_4_days_blgr <- unmarkedFrameOccu(
  y = as.matrix(every_4_days_blgr$y),
  siteCovs = every_4_days_blgr$abund.covs,
  obsCovs = list(
    day = as.matrix(every_4_days_blgr$det.covs$day),
    temp = as.matrix(every_4_days_blgr$det.covs$temp),
    wind_speed = as.matrix(every_4_days_blgr$det.covs$wind_speed),
    precipitation = as.matrix(every_4_days_blgr$det.covs$precipitation)
  )
)

every_4_days_blgr <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ factor(treatment) - 1 + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_every_4_days_blgr, 
  K = 25,  
  method = "BFGS"  
)

# Every 7 days
umf_every_7_days_blgr <- unmarkedFrameOccu(
  y = as.matrix(every_7_days_blgr$y),
  siteCovs = every_7_days_blgr$abund.covs,
  obsCovs = list(
    day = as.matrix(every_7_days_blgr$det.covs$day),
    temp = as.matrix(every_7_days_blgr$det.covs$temp),
    wind_speed = as.matrix(every_7_days_blgr$det.covs$wind_speed),
    precipitation = as.matrix(every_7_days_blgr$det.covs$precipitation)
  )
)

every_7_days_blgr <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ factor(treatment) - 1 + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_every_7_days_blgr, 
  K = 25,  
  method = "BFGS"  
)

# Every 10 days
umf_every_10_days_blgr <- unmarkedFrameOccu(
  y = as.matrix(every_10_days_blgr$y),
  siteCovs = every_10_days_blgr$abund.covs,
  obsCovs = list(
    day = as.matrix(every_10_days_blgr$det.covs$day),
    temp = as.matrix(every_10_days_blgr$det.covs$temp),
    wind_speed = as.matrix(every_10_days_blgr$det.covs$wind_speed),
    precipitation = as.matrix(every_10_days_blgr$det.covs$precipitation)
  )
)

every_10_days_blgr <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ factor(treatment) - 1 + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_every_10_days_blgr, 
  K = 25,  
  method = "BFGS"  
)


#Save all models 
save(every_day_blgr, every_2_days_blgr, every_3_days_blgr, every_4_days_blgr, every_7_days_blgr, every_10_days_blgr, file = "data/rn_models/rn_subsample_models_blgr.RData")


##GOF test 
gof_blgr <- mb.gof.test(every_day_blgr, nsim=200, c.hat.est=TRUE, model.type="royle-nichols")
print(gof_blgr)

##########Predictions for effect of treatment for every day##########

# Create a new data frame with the mean values of the covariates
newdata_blgr <- data.frame(treatment=levels(umf_every_day_blgr@siteCovs$treatment), shrub_cover=mean(umf_every_day_blgr@siteCovs$shrub_cover),grass_cover=mean(umf_every_day_blgr@siteCovs$grass_cover), forb_cover=mean(umf_every_day_blgr@siteCovs$forb_cover),shrub_height=mean(umf_every_day_blgr@siteCovs$shrub_height), grass_height=mean(umf_every_day_blgr@siteCovs$grass_height))

# Generate predictions for the new data with the predict function
predictions_blgr <- predict(every_day_blgr, newdata_blgr, type = "state", se.fit = TRUE)

# Add predictions and SE to the newdata_blgr dataframe
newdata_blgr$predicted_state <- predictions_blgr$Predicted
newdata_blgr$SE <- predictions_blgr$SE

# Calculate 95% confidence intervals from SE
newdata_blgr$lower_CI <- newdata_blgr$predicted_state - 1.96 * newdata_blgr$SE
newdata_blgr$upper_CI <- newdata_blgr$predicted_state + 1.96 * newdata_blgr$SE

# View the results
print(newdata_blgr)

# Save results as csv for plotting in another R session
write.csv(newdata_blgr, file = "data/means_abund_parameters/means_treatment_parameters_blgr_1.csv", row.names = FALSE)

##########Predictions for effect of treatment for every 2 days##########
newdata_2_blgr <- data.frame(treatment=levels(umf_every_2_days_blgr@siteCovs$treatment), shrub_cover=mean(umf_every_2_days_blgr@siteCovs$shrub_cover),grass_cover=mean(umf_every_2_days_blgr@siteCovs$grass_cover), forb_cover=mean(umf_every_2_days_blgr@siteCovs$forb_cover),shrub_height=mean(umf_every_2_days_blgr@siteCovs$shrub_height), grass_height=mean(umf_every_2_days_blgr@siteCovs$grass_heigh))

predictions_2_blgr <- predict(every_2_days_blgr, newdata_2_blgr, type = "state", se.fit = TRUE)

newdata_2_blgr$predicted_state <- predictions_2_blgr$Predicted
newdata_2_blgr$SE <- predictions_2_blgr$SE

newdata_2_blgr$lower_CI <- newdata_2_blgr$predicted_state - 1.96 * newdata_2_blgr$SE
newdata_2_blgr$upper_CI <- newdata_2_blgr$predicted_state + 1.96 * newdata_2_blgr$SE

print(newdata_2_blgr)

write.csv(newdata_2_blgr, file = "data/means_abund_parameters/means_treatment_parameters_blgr_2.csv", row.names = FALSE)

##########Predictions for effect of treatment for every 3 days##########
newdata_3_blgr <- data.frame(treatment=levels(umf_every_3_days_blgr@siteCovs$treatment), shrub_cover=mean(umf_every_3_days_blgr@siteCovs$shrub_cover),grass_cover=mean(umf_every_3_days_blgr@siteCovs$grass_cover), forb_cover=mean(umf_every_3_days_blgr@siteCovs$forb_cover),shrub_height=mean(umf_every_3_days_blgr@siteCovs$shrub_height), grass_height=mean(umf_every_3_days_blgr@siteCovs$grass_height))

predictions_3_blgr <- predict(every_3_days_blgr, newdata_3_blgr, type = "state", se.fit = TRUE)

newdata_3_blgr$predicted_state <- predictions_3_blgr$Predicted
newdata_3_blgr$SE <- predictions_3_blgr$SE

newdata_3_blgr$lower_CI <- newdata_3_blgr$predicted_state - 1.96 * newdata_3_blgr$SE
newdata_3_blgr$upper_CI <- newdata_3_blgr$predicted_state + 1.96 * newdata_3_blgr$SE

print(newdata_3_blgr)

write.csv(newdata_3_blgr, file = "data/means_abund_parameters/means_treatment_parameters_blgr_3.csv", row.names = FALSE)

##########Predictions for effect of treatment for every 4 days##########
newdata_4_blgr <- data.frame(treatment=levels(umf_every_4_days_blgr@siteCovs$treatment), shrub_cover=mean(umf_every_4_days_blgr@siteCovs$shrub_cover),grass_cover=mean(umf_every_4_days_blgr@siteCovs$grass_cover), forb_cover=mean(umf_every_4_days_blgr@siteCovs$forb_cover),shrub_height=mean(umf_every_4_days_blgr@siteCovs$shrub_height), grass_height=mean(umf_every_4_days_blgr@siteCovs$grass_height))

predictions_4_blgr <- predict(every_4_days_blgr, newdata_4_blgr, type = "state", se.fit = TRUE)

newdata_4_blgr$predicted_state <- predictions_4_blgr$Predicted
newdata_4_blgr$SE <- predictions_4_blgr$SE

newdata_4_blgr$lower_CI <- newdata_4_blgr$predicted_state - 1.96 * newdata_4_blgr$SE
newdata_4_blgr$upper_CI <- newdata_4_blgr$predicted_state + 1.96 * newdata_4_blgr$SE

print(newdata_4_blgr)

write.csv(newdata_4_blgr, file = "data/means_abund_parameters/means_treatment_parameters_blgr_4.csv", row.names = FALSE)

##########Predictions for effect of treatment for every 7 days##########
newdata_7_blgr <- data.frame(treatment=levels(umf_every_7_days_blgr@siteCovs$treatment), shrub_cover=mean(umf_every_7_days_blgr@siteCovs$shrub_cover),grass_cover=mean(umf_every_7_days_blgr@siteCovs$grass_cover), forb_cover=mean(umf_every_7_days_blgr@siteCovs$forb_cover),shrub_height=mean(umf_every_7_days_blgr@siteCovs$shrub_height), grass_height=mean(umf_every_7_days_blgr@siteCovs$grass_height))

predictions_7_blgr <- predict(every_7_days_blgr, newdata_7_blgr, type = "state", se.fit = TRUE)

newdata_7_blgr$predicted_state <- predictions_7_blgr$Predicted
newdata_7_blgr$SE <- predictions_7_blgr$SE

newdata_7_blgr$lower_CI <- newdata_7_blgr$predicted_state - 1.96 * newdata_7_blgr$SE
newdata_7_blgr$upper_CI <- newdata_7_blgr$predicted_state + 1.96 * newdata_7_blgr$SE

print(newdata_7_blgr)

write.csv(newdata_7_blgr, file = "data/means_abund_parameters/means_treatment_parameters_blgr_7.csv", row.names = FALSE)

##########Predictions for effect of treatment for every 10 days##########
newdata_10_blgr <- data.frame(treatment=levels(umf_every_10_days_blgr@siteCovs$treatment), shrub_cover=mean(umf_every_10_days_blgr@siteCovs$shrub_cover),grass_cover=mean(umf_every_10_days_blgr@siteCovs$grass_cover), forb_cover=mean(umf_every_10_days_blgr@siteCovs$forb_cover),shrub_height=mean(umf_every_10_days_blgr@siteCovs$shrub_height), grass_height=mean(umf_every_10_days_blgr@siteCovs$grass_height))

predictions_10_blgr <- predict(every_10_days_blgr, newdata_10_blgr, type = "state", se.fit = TRUE)

newdata_10_blgr$predicted_state <- predictions_10_blgr$Predicted
newdata_10_blgr$SE <- predictions_10_blgr$SE

newdata_10_blgr$lower_CI <- newdata_10_blgr$predicted_state - 1.96 * newdata_10_blgr$SE
newdata_10_blgr$upper_CI <- newdata_10_blgr$predicted_state + 1.96 * newdata_10_blgr$SE

print(newdata_10_blgr)

write.csv(newdata_10_blgr, file = "data/means_abund_parameters/means_treatment_parameters_blgr_10.csv", row.names = FALSE)


#############################Repeat for Bachman's Sparrow######################################

# Load bacs abundance data
load("data/abundance_data/abundance_24_bacs.RData")

# Restrict to columns 56 to 105
abund_data_bacs$y <- abund_data_bacs$y[, 56:105]
abund_data_bacs$det.covs$day <- abund_data_bacs$det.covs$day[, 56:105]
abund_data_bacs$det.covs$temp <- abund_data_bacs$det.covs$temp[, 56:105]
abund_data_bacs$det.covs$wind_speed <- abund_data_bacs$det.covs$wind_speed[, 56:105]
abund_data_bacs$det.covs$precipitation <- abund_data_bacs$det.covs$precipitation[, 56:105]

# Convert values greater than 0 or NA to 1 in "y"
abund_data_bacs$y <- ifelse(abund_data_bacs$y > 0 | is.na(abund_data_bacs$y), 1, 0)

# Extract detection covariates
det.covs <- list(
  day = abund_data_bacs$det.covs$day,
  temp = abund_data_bacs$det.covs$temp,
  wind_speed = abund_data_bacs$det.covs$wind_speed,
  precipitation = abund_data_bacs$det.covs$precipitation
)

## Get the actual starting column index
start_col_label <- "56"
start_col <- which(colnames(abund_data_bacs$y) == start_col_label)

# Retain only every nth day
retain_every_nth <- function(data, start_col, n) {
  cols_to_keep <- seq(start_col, ncol(data), by = n)
  return(data[, cols_to_keep])
}

# Define the intervals
n_values <- c(2, 3, 4, 7, 10)

# Create datasets for each interval
datasets <- list()

for (n in n_values) {
  # Retain only the desired columns for y and det.covs
  reduced_y <- retain_every_nth(abund_data_bacs$y, start_col, n)
  reduced_det.covs <- lapply(det.covs, function(cov) retain_every_nth(cov, start_col, n))
  
  # Store the reduced dataset in the list
  datasets[[paste0("every_", n, "_days","_bacs")]] <- list(
    y = reduced_y,
    det.covs = reduced_det.covs,
    abund.covs = abund_data_bacs$abund.covs,
    coords = abund_data_bacs$coords
  )
}

# Save each dataset in the list as an RData file with the correct name
for (name in names(datasets)) {
  # Extract the dataset from the list
  assign(name, datasets[[name]])  # Assign dataset to a variable with the desired name
  
  # Save the dataset with its unique name
  save(list = name, file = paste0("data/abundance_data/bacs_subsamples/", name, ".RData"))
}

# Load RData files for the new intervals
load("data/abundance_data/bacs_subsamples/every_2_days_bacs.RData")
load("data/abundance_data/bacs_subsamples/every_3_days_bacs.RData")
load("data/abundance_data/bacs_subsamples/every_4_days_bacs.RData")
load("data/abundance_data/bacs_subsamples/every_7_days_bacs.RData")
load("data/abundance_data/bacs_subsamples/every_10_days_bacs.RData")

#Rename to every_day_bacs
every_day_bacs <- abund_data_bacs

#umf for every_day
umf_every_day_bacs <- unmarkedFrameOccu(
  y = as.matrix(every_day_bacs$y),
  siteCovs = every_day_bacs$abund.covs,
  obsCovs = list(
    day = as.matrix(every_day_bacs$det.covs$day),
    temp = as.matrix(every_day_bacs$det.covs$temp),
    wind_speed = as.matrix(every_day_bacs$det.covs$wind_speed),
    precipitation = as.matrix(every_day_bacs$det.covs$precipitation)
  )
)

#Global model
every_day_bacs <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_every_day_bacs, 
  K = 25,  
  method = "BFGS"  
)

#Summary
summary(every_day_bacs)

#Just treatment
every_day_t_bacs <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ factor(treatment) -1,
  data = umf_every_day_bacs, 
  K = 25,  
  method = "BFGS"  
)

#Just vegetation
every_day_v_bacs <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_every_day_bacs, 
  K = 25,  
  method = "BFGS"  
)

#Null Model
every_day_null_bacs <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ 1,
  data = umf_every_day_bacs, 
  K = 25,  
  method = "BFGS"  
)

########AICc comparison and AICc table formatting###########

models_list_bacs <- list(every_day_bacs = every_day_bacs, 
                         every_day_t_bacs = every_day_t_bacs, 
                         every_day_v_bacs = every_day_v_bacs, 
                         every_day_null_bacs = every_day_null_bacs)

model_names_bacs <- c("TV", "T", "V", "Null")
aicc_table_bacs <- aictab(cand.set = models_list_bacs, modnames = model_names_bacs)
print(aicc_table_bacs)


# Every 2 days
umf_every_2_days_bacs <- unmarkedFrameOccu(
  y = as.matrix(every_2_days_bacs$y),
  siteCovs = every_2_days_bacs$abund.covs,
  obsCovs = list(
    day = as.matrix(every_2_days_bacs$det.covs$day),
    temp = as.matrix(every_2_days_bacs$det.covs$temp),
    wind_speed = as.matrix(every_2_days_bacs$det.covs$wind_speed),
    precipitation = as.matrix(every_2_days_bacs$det.covs$precipitation)
  )
)

every_2_days_bacs <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ factor(treatment) - 1 + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_every_2_days_bacs, 
  K = 25,  
  method = "BFGS"  
)

# Every 3 days
umf_every_3_days_bacs <- unmarkedFrameOccu(
  y = as.matrix(every_3_days_bacs$y),
  siteCovs = every_3_days_bacs$abund.covs,
  obsCovs = list(
    day = as.matrix(every_3_days_bacs$det.covs$day),
    temp = as.matrix(every_3_days_bacs$det.covs$temp),
    wind_speed = as.matrix(every_3_days_bacs$det.covs$wind_speed),
    precipitation = as.matrix(every_3_days_bacs$det.covs$precipitation)
  )
)

every_3_days_bacs <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ factor(treatment) - 1 + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_every_3_days_bacs, 
  K = 25,  
  method = "BFGS"  
)

# Every 4 days
umf_every_4_days_bacs <- unmarkedFrameOccu(
  y = as.matrix(every_4_days_bacs$y),
  siteCovs = every_4_days_bacs$abund.covs,
  obsCovs = list(
    day = as.matrix(every_4_days_bacs$det.covs$day),
    temp = as.matrix(every_4_days_bacs$det.covs$temp),
    wind_speed = as.matrix(every_4_days_bacs$det.covs$wind_speed),
    precipitation = as.matrix(every_4_days_bacs$det.covs$precipitation)
  )
)

every_4_days_bacs <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ factor(treatment) - 1 + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_every_4_days_bacs, 
  K = 25,  
  method = "BFGS"  
)

# Every 7 days
umf_every_7_days_bacs <- unmarkedFrameOccu(
  y = as.matrix(every_7_days_bacs$y),
  siteCovs = every_7_days_bacs$abund.covs,
  obsCovs = list(
    day = as.matrix(every_7_days_bacs$det.covs$day),
    temp = as.matrix(every_7_days_bacs$det.covs$temp),
    wind_speed = as.matrix(every_7_days_bacs$det.covs$wind_speed),
    precipitation = as.matrix(every_7_days_bacs$det.covs$precipitation)
  )
)

every_7_days_bacs <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ factor(treatment) - 1 + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_every_7_days_bacs, 
  K = 25,  
  method = "BFGS"  
)

# Every 10 days
umf_every_10_days_bacs <- unmarkedFrameOccu(
  y = as.matrix(every_10_days_bacs$y),
  siteCovs = every_10_days_bacs$abund.covs,
  obsCovs = list(
    day = as.matrix(every_10_days_bacs$det.covs$day),
    temp = as.matrix(every_10_days_bacs$det.covs$temp),
    wind_speed = as.matrix(every_10_days_bacs$det.covs$wind_speed),
    precipitation = as.matrix(every_10_days_bacs$det.covs$precipitation)
  )
)

every_10_days_bacs <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed)
  ~ factor(treatment) - 1 + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_every_10_days_bacs, 
  K = 25,  
  method = "BFGS"  
)

#Save all models 
save(every_day_bacs, every_2_days_bacs, every_3_days_bacs, every_4_days_bacs, every_7_days_bacs, every_10_days_bacs, file = "data/rn_models/rn_subsample_models_bacs.RData")

##GOF test 
gof_bacs <- mb.gof.test(every_day_bacs, nsim=200, c.hat.est=TRUE, model.type="royle-nichols")
print(gof_bacs)

##########Predictions for effect of treatment for every day##########

# Create a new data frame with the mean values of the covariates
newdata_bacs <- data.frame(treatment=levels(umf_every_day_bacs@siteCovs$treatment), shrub_cover=mean(umf_every_day_bacs@siteCovs$shrub_cover),grass_cover=mean(umf_every_day_bacs@siteCovs$grass_cover), forb_cover=mean(umf_every_day_bacs@siteCovs$forb_cover),shrub_height=mean(umf_every_day_bacs@siteCovs$shrub_height), grass_height=mean(umf_every_day_bacs@siteCovs$grass_height))

# Generate predictions for the new data with the predict function
predictions_bacs <- predict(every_day_bacs, newdata_bacs, type = "state", se.fit = TRUE)

# Add predictions and SE to the newdata_bacs dataframe
newdata_bacs$predicted_state <- predictions_bacs$Predicted
newdata_bacs$SE <- predictions_bacs$SE

# Calculate 95% confidence intervals from SE
newdata_bacs$lower_CI <- newdata_bacs$predicted_state - 1.96 * newdata_bacs$SE
newdata_bacs$upper_CI <- newdata_bacs$predicted_state + 1.96 * newdata_bacs$SE

# View the results
print(newdata_bacs)

# Save results as csv for plotting in another R session
write.csv(newdata_bacs, file = "data/means_abund_parameters/means_treatment_parameters_bacs_1.csv", row.names = FALSE)

##########Predictions for effect of treatment for every 2 days##########
newdata_2_bacs <- data.frame(treatment=levels(umf_every_2_days_bacs@siteCovs$treatment), shrub_cover=mean(umf_every_2_days_bacs@siteCovs$shrub_cover),grass_cover=mean(umf_every_2_days_bacs@siteCovs$grass_cover), forb_cover=mean(umf_every_2_days_bacs@siteCovs$forb_cover),shrub_height=mean(umf_every_2_days_bacs@siteCovs$shrub_height), grass_height=mean(umf_every_2_days_bacs@siteCovs$grass_heigh))

predictions_2_bacs <- predict(every_2_days_bacs, newdata_2_bacs, type = "state", se.fit = TRUE)

newdata_2_bacs$predicted_state <- predictions_2_bacs$Predicted
newdata_2_bacs$SE <- predictions_2_bacs$SE

newdata_2_bacs$lower_CI <- newdata_2_bacs$predicted_state - 1.96 * newdata_2_bacs$SE
newdata_2_bacs$upper_CI <- newdata_2_bacs$predicted_state + 1.96 * newdata_2_bacs$SE

print(newdata_2_bacs)

write.csv(newdata_2_bacs, file = "data/means_abund_parameters/means_treatment_parameters_bacs_2.csv", row.names = FALSE)

##########Predictions for effect of treatment for every 3 days##########
newdata_3_bacs <- data.frame(treatment=levels(umf_every_3_days_bacs@siteCovs$treatment), shrub_cover=mean(umf_every_3_days_bacs@siteCovs$shrub_cover),grass_cover=mean(umf_every_3_days_bacs@siteCovs$grass_cover), forb_cover=mean(umf_every_3_days_bacs@siteCovs$forb_cover),shrub_height=mean(umf_every_3_days_bacs@siteCovs$shrub_height), grass_height=mean(umf_every_3_days_bacs@siteCovs$grass_height))

predictions_3_bacs <- predict(every_3_days_bacs, newdata_3_bacs, type = "state", se.fit = TRUE)

newdata_3_bacs$predicted_state <- predictions_3_bacs$Predicted
newdata_3_bacs$SE <- predictions_3_bacs$SE

newdata_3_bacs$lower_CI <- newdata_3_bacs$predicted_state - 1.96 * newdata_3_bacs$SE
newdata_3_bacs$upper_CI <- newdata_3_bacs$predicted_state + 1.96 * newdata_3_bacs$SE

print(newdata_3_bacs)

write.csv(newdata_3_bacs, file = "data/means_abund_parameters/means_treatment_parameters_bacs_3.csv", row.names = FALSE)

##########Predictions for effect of treatment for every 4 days##########
newdata_4_bacs <- data.frame(treatment=levels(umf_every_4_days_bacs@siteCovs$treatment), shrub_cover=mean(umf_every_4_days_bacs@siteCovs$shrub_cover),grass_cover=mean(umf_every_4_days_bacs@siteCovs$grass_cover), forb_cover=mean(umf_every_4_days_bacs@siteCovs$forb_cover),shrub_height=mean(umf_every_4_days_bacs@siteCovs$shrub_height), grass_height=mean(umf_every_4_days_bacs@siteCovs$grass_height))

predictions_4_bacs <- predict(every_4_days_bacs, newdata_4_bacs, type = "state", se.fit = TRUE)

newdata_4_bacs$predicted_state <- predictions_4_bacs$Predicted
newdata_4_bacs$SE <- predictions_4_bacs$SE

newdata_4_bacs$lower_CI <- newdata_4_bacs$predicted_state - 1.96 * newdata_4_bacs$SE
newdata_4_bacs$upper_CI <- newdata_4_bacs$predicted_state + 1.96 * newdata_4_bacs$SE

print(newdata_4_bacs)

write.csv(newdata_4_bacs, file = "data/means_abund_parameters/means_treatment_parameters_bacs_4.csv", row.names = FALSE)

##########Predictions for effect of treatment for every 7 days##########
newdata_7_bacs <- data.frame(treatment=levels(umf_every_7_days_bacs@siteCovs$treatment), shrub_cover=mean(umf_every_7_days_bacs@siteCovs$shrub_cover),grass_cover=mean(umf_every_7_days_bacs@siteCovs$grass_cover), forb_cover=mean(umf_every_7_days_bacs@siteCovs$forb_cover),shrub_height=mean(umf_every_7_days_bacs@siteCovs$shrub_height), grass_height=mean(umf_every_7_days_bacs@siteCovs$grass_height))

predictions_7_bacs <- predict(every_7_days_bacs, newdata_7_bacs, type = "state", se.fit = TRUE)

newdata_7_bacs$predicted_state <- predictions_7_bacs$Predicted
newdata_7_bacs$SE <- predictions_7_bacs$SE

newdata_7_bacs$lower_CI <- newdata_7_bacs$predicted_state - 1.96 * newdata_7_bacs$SE
newdata_7_bacs$upper_CI <- newdata_7_bacs$predicted_state + 1.96 * newdata_7_bacs$SE

print(newdata_7_bacs)

write.csv(newdata_7_bacs, file = "data/means_abund_parameters/means_treatment_parameters_bacs_7.csv", row.names = FALSE)

##########Predictions for effect of treatment for every 10 days##########
newdata_10_bacs <- data.frame(treatment=levels(umf_every_10_days_bacs@siteCovs$treatment), shrub_cover=mean(umf_every_10_days_bacs@siteCovs$shrub_cover),grass_cover=mean(umf_every_10_days_bacs@siteCovs$grass_cover), forb_cover=mean(umf_every_10_days_bacs@siteCovs$forb_cover),shrub_height=mean(umf_every_10_days_bacs@siteCovs$shrub_height), grass_height=mean(umf_every_10_days_bacs@siteCovs$grass_height))

predictions_10_bacs <- predict(every_10_days_bacs, newdata_10_bacs, type = "state", se.fit = TRUE)

newdata_10_bacs$predicted_state <- predictions_10_bacs$Predicted
newdata_10_bacs$SE <- predictions_10_bacs$SE

newdata_10_bacs$lower_CI <- newdata_10_bacs$predicted_state - 1.96 * newdata_10_bacs$SE
newdata_10_bacs$upper_CI <- newdata_10_bacs$predicted_state + 1.96 * newdata_10_bacs$SE

print(newdata_10_bacs)

write.csv(newdata_10_bacs, file = "data/means_abund_parameters/means_treatment_parameters_bacs_10.csv", row.names = FALSE)