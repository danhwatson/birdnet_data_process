# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  

# Load cgdo abundance data
load("data/abundance_data/abundance_24_cgdo.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_cgdo$y <- ifelse(abund_data_cgdo$y > 0 | is.na(abund_data_cgdo$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as RN model can't handle NAs
abund_data_cgdo$y <- abund_data_cgdo$y[, -c(1:55)]
abund_data_cgdo$det.covs$day <- abund_data_cgdo$det.covs$day[, -c(1:55)]
abund_data_cgdo$det.covs$temp <- abund_data_cgdo$det.covs$temp[, -c(1:55)]
abund_data_cgdo$det.covs$wind_speed <- abund_data_cgdo$det.covs$wind_speed[, -c(1:55)]
abund_data_cgdo$det.covs$precipitation <- abund_data_cgdo$det.covs$precipitation[, -c(1:55)]


# Create the unmarkedFrameOccu object
umf_cgdo <- unmarkedFrameOccu(
  y = as.matrix(abund_data_cgdo$y),
  siteCovs = abund_data_cgdo$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_cgdo$det.covs$day),
    temp = as.matrix(abund_data_cgdo$det.covs$temp),
    wind_speed = as.matrix(abund_data_cgdo$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_cgdo$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model
rn_model_t_cgdo <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ factor(treatment) -1 + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_cgdo,
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)
summary(rn_model_t_cgdo)
saveRDS(rn_model_t_cgdo, file = "data/rn_model_t_cgdo.rds")

# Fit the second Royle-Nichols model
rn_model_cgdo <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_cgdo,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_cgdo)
saveRDS(rn_model_cgdo, file = "data/rn_model_cgdo.rds")

# Fit the null model
rn_model_null_cgdo <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ 1, 
  data = umf_cgdo,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_cgdo)

# AICc comparison
models_list_cgdo <- list(rn_model_t_cgdo = rn_model_t_cgdo, rn_model_cgdo = rn_model_cgdo, rn_model_null_cgdo = rn_model_null_cgdo)
model_names_cgdo <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_cgdo <- aictab(cand.set = models_list_cgdo, modnames = model_names_cgdo)
print(aicc_table_cgdo)

# Function returning fit-statistics
#fitstats_cgdo <- function(rn_model_t_cgdo) {
#  observed <- getY(rn_model_t_cgdo@data)
#  expected <- fitted(rn_model_t_cgdo)
#  resids <- residuals(rn_model_t_cgdo)
#  sse <- sum(resids^2, na.rm=TRUE)
#  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
#  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
#  fit_out <- c(cSSE = sse, Chisq = chisq, freemanTukey = freeTuke)
#  return(fit_out)
#}

gof <- mb.gof.test(rn_model_t_cgdo, nsim=100, c.hat.est=TRUE, model.type="RN")

#parametric bootstrapping 
# Assuming 'mod' is your fitted model
#pb_cgdo <- parboot(rn_model_t_cgdo, fitstats_cgdo, nsim=100)
#print(pb_cgdo)

# Predictions
newdata_cgdo <- data.frame(treatment=levels(umf_cgdo@siteCovs$treatment), shrub_cover=mean(umf_cgdo@siteCovs$shrub_cover),grass_cover=mean(umf_cgdo@siteCovs$grass_cover), forb_cover=mean(umf_cgdo@siteCovs$forb_cover),shrub_height=mean(umf_cgdo@siteCovs$shrub_height), grass_height=mean(umf_cgdo@siteCovs$grass_height))
#

# Predictions for the new data
predictions_cgdo <- predict(rn_model_t_cgdo, newdata_cgdo, type = "state", se.fit = TRUE)

# Add predictions and SE to the newdata_cgdo dataframe
newdata_cgdo$predicted_state <- predictions_cgdo$Predicted
newdata_cgdo$SE <- predictions_cgdo$SE

# Calculate 95% confidence intervals
newdata_cgdo$lower_CI <- newdata_cgdo$predicted_state - 1.96 * newdata_cgdo$SE
newdata_cgdo$upper_CI <- newdata_cgdo$predicted_state + 1.96 * newdata_cgdo$SE

# View the results
print(newdata_cgdo)

# Save results as csv
write.csv(newdata_cgdo, file = "data/mean_parameters_cgdo.csv", row.names = FALSE)



# Optionally, save the dataframe
saveRDS(newdata_cgdo, file = "data/newdata_cgdo_with_CI.rds")

# Predictions
# Generate new data for shrub cover while keeping other covariates constant
shrub_range <- seq(min(umf_cgdo@siteCovs$shrub_cover), max(umf_cgdo@siteCovs$shrub_cover), length.out = 100)
newdata_shrub <- data.frame(
  treatment = factor(1, levels = c(1, 2, 3, 4)),  # Choose one treatment level, adjust as needed
  shrub_cover = shrub_range,
  grass_cover = mean(umf_cgdo@siteCovs$grass_cover),
  forb_cover = mean(umf_cgdo@siteCovs$forb_cover),
  shrub_height = mean(umf_cgdo@siteCovs$shrub_height),
  grass_height = mean(umf_cgdo@siteCovs$grass_height)
)

# Check the levels of the treatment factor in the original data
treatment_levels <- levels(umf_cgdo@siteCovs$treatment)

# Generate new data for shrub cover while keeping other covariates constant
shrub_range <- seq(min(umf_cgdo@siteCovs$shrub_cover), max(umf_cgdo@siteCovs$shrub_cover), length.out = 100)

# Create a new data frame with shrub cover and other covariates constant for each treatment level
newdata_shrub <- expand.grid(
  treatment = factor(treatment_levels, levels = treatment_levels),
  shrub_cover = shrub_range,
  grass_cover = mean(umf_cgdo@siteCovs$grass_cover),
  forb_cover = mean(umf_cgdo@siteCovs$forb_cover),
  shrub_height = mean(umf_cgdo@siteCovs$shrub_height),
  grass_height = mean(umf_cgdo@siteCovs$grass_height)
)

# Predict state (occupancy or abundance) for shrub cover range across all treatment levels
predictions_shrub <- predict(rn_model_t_cgdo, newdata_shrub, type = "state")

# Create a dataframe for plotting
plot_data_shrub <- data.frame(
  shrub_cover = newdata_shrub$shrub_cover,
  treatment = newdata_shrub$treatment,
  predicted_state = predictions_shrub$Predicted,
  lower_CI = predictions_shrub$lower,
  upper_CI = predictions_shrub$upper
)

# Plot the predictive response to shrub cover for each treatmentx
ggplot(plot_data_shrub, aes(x = shrub_cover, y = predicted_state, color = treatment)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill = treatment), alpha = 0.2) +
  labs(x = "Shrub Cover", y = "Predicted State", title = "Predictive Response to Shrub Cover by Treatment") +
  theme_minimal()


# Predict state (occupancy or abundance) for shrub cover range regardless of treatment 
predictions_shrub_all <- predict(rn_model_t_cgdo, newdata_shrub, type = "state")

# Create a dataframe for plotting
plot_data_shrub_all <- data.frame(
  shrub_cover = newdata_shrub$shrub_cover,
  predicted_state = predictions_shrub_all$Predicted,
  lower_CI = predictions_shrub_all$lower,
  upper_CI = predictions_shrub_all$upper
)

# Plot the predictive response to shrub cover for all treatments
ggplot(plot_data_shrub_all, aes(x = shrub_cover, y = predicted_state)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill = predicted), alpha = 0.2) +
  labs(x = "Shrub Cover", y = "Predicted State", title = "Predictive Response to Shrub Cover") +
  theme_minimal()



