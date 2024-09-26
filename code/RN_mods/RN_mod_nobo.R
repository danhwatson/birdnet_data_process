# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  


##########Loading data and formatting for RN mods###########

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


##########Fitting Royle-Nichols models###########

# Create the unmarkedFrameOccu object needed to be used for all models made with unmarked
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

# Fit the first Royle-Nichols model (with site treatment as an abundance covariate, t for treatment)
rn_model_t_nobo <- occuRN( # This is the unmarked function to fit Royle-Nichols (RN) models
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) # Detection model covariates
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), # Abundance model covariates 
  data = umf_nobo, # Specifying the data from the unmarkedFrameOccu object we created
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)

summary(rn_model_t_nobo) # Look at the summary output of the model
# Models have '-1' term in the abundance model formula to remove the intercept, therefore providing mean parameter estimates for each level of the factor (treatment)

#saveRDS(rn_model_t_nobo, file = "data/rn_model_t_nobo.rds") 
# Optional step to save model object to load in another R session


# Fit the second Royle-Nichols model (without site treatment as an abundance covariate)
rn_model_nobo <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_nobo,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_nobo)

# Fit the null model
rn_model_null_nobo <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
  ~ 1, 
  data = umf_nobo,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_nobo)

########AICc comparison###########

# Load the master AICc table
aicc_table <- read.csv("data/aicc_table.csv")

models_list_nobo <- list(rn_model_t_nobo = rn_model_t_nobo, rn_model_nobo = rn_model_nobo, rn_model_null_nobo = rn_model_null_nobo)
model_names_nobo <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_nobo <- aictab(cand.set = models_list_nobo, modnames = model_names_nobo)
print(aicc_table_nobo)

# Add species name to the table
aicc_table_nobo$species_code <- "NOBO"

# Add the AICc table to the master table
aicc_table_master <- rbind(aicc_table, aicc_table_nobo)

# Save the AICc table to csv
write.csv(aicc_table_master, "data/aicc_table.csv", row.names = FALSE)


##########Goodness of fit##########

# Function returning fit-statistics
#fitstats_nobo <- function(rn_model_t_nobo) {
#  observed <- getY(rn_model_t_nobo@data)
#  expected <- fitted(rn_model_t_nobo)
#  resids <- residuals(rn_model_t_nobo)
#  sse <- sum(resids^2, na.rm=TRUE)
#  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
#  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
#  fit_out <- c(cSSE = sse, Chisq = chisq, freemanTukey = freeTuke)
#  return(fit_out)
#}


# Need to look more into GoF tests
#gof <- mb.gof.test(rn_model_t_nobo, nsim=100, c.hat.est=TRUE, model.type="RN")

#parametric bootstrapping 
# Assuming 'mod' is your fitted model
#pb_nobo <- parboot(rn_model_t_nobo, fitstats_nobo, nsim=100)
#print(pb_nobo)


##########Predictions for effect of treatment##########

# Create a new data frame with the mean values of the covariates
newdata_nobo <- data.frame(treatment=levels(umf_nobo@siteCovs$treatment), shrub_cover=mean(umf_nobo@siteCovs$shrub_cover),grass_cover=mean(umf_nobo@siteCovs$grass_cover), forb_cover=mean(umf_nobo@siteCovs$forb_cover),shrub_height=mean(umf_nobo@siteCovs$shrub_height), grass_height=mean(umf_nobo@siteCovs$grass_height))

# Generate predictions for the new data with the predict function
predictions_nobo <- predict(rn_model_t_nobo, newdata_nobo, type = "state", se.fit = TRUE)

# Add predictions and SE to the newdata_nobo dataframe
newdata_nobo$predicted_state <- predictions_nobo$Predicted
newdata_nobo$SE <- predictions_nobo$SE

# Calculate 95% confidence intervals from SE
newdata_nobo$lower_CI <- newdata_nobo$predicted_state - 1.96 * newdata_nobo$SE
newdata_nobo$upper_CI <- newdata_nobo$predicted_state + 1.96 * newdata_nobo$SE

# View the results
print(newdata_nobo)

# Optionally, save results as csv for plotting in another R session
write.csv(newdata_nobo, file = "data/means_treatment_parameters_nobo.csv", row.names = FALSE)


##########Predictions for effect of grass cover##########

# Generate new data for grass cover while keeping other covariates constant
grass_range <- seq(min(umf_nobo@siteCovs$grass_cover), max(umf_nobo@siteCovs$grass_cover), length.out = 100)

# Check the levels of the treatment factor in the original data
treatment_levels <- levels(umf_nobo@siteCovs$treatment)

newdata_grass <- data.frame(
  treatment = factor(treatment_levels, levels = treatment_levels),
  grass_cover = grass_range,
  shrub_cover = mean(umf_nobo@siteCovs$shrub_cover),
  forb_cover = mean(umf_nobo@siteCovs$forb_cover),
  shrub_height = mean(umf_nobo@siteCovs$shrub_height),
  grass_height = mean(umf_nobo@siteCovs$grass_height)
)

# Predict state (occupancy or abundance) for grass  cover range across all treatment levels
predictions_grass <- predict(rn_model_t_nobo, newdata_grass, type = "state")

# Create a dataframe for plotting
plot_data_grass <- data.frame(
  grass_cover = newdata_grass$grass_cover,
  treatment = newdata_grass$treatment,
  predicted_state = predictions_grass$Predicted,
  lower_CI = predictions_grass$lower,
  upper_CI = predictions_grass$upper
)

# Plot the predictive response to grass cover for each treatment
ggplot(plot_data_grass, aes(x = grass_cover, y = predicted_state, color = treatment)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill = treatment), alpha = 0.2) +
  labs(x = "Grass Cover", y = "Predicted Abundance", title = "Northern Bobwhite (nobo) -- Predicted Response to Grass Cover across Treatments") +
  theme_minimal()

# Generate new data for grass cover without the treatment variable
grass_range <- seq(min(umf_nobo@siteCovs$grass_cover), max(umf_nobo@siteCovs$grass_cover), length.out = 100)

newdata_grass <- data.frame(
  grass_cover = grass_range,
  shrub_cover = mean(umf_nobo@siteCovs$shrub_cover),
  forb_cover = mean(umf_nobo@siteCovs$forb_cover),
  shrub_height = mean(umf_nobo@siteCovs$shrub_height),
  grass_height = mean(umf_nobo@siteCovs$grass_height)
)

# Predict state (occupancy or abundance) for grass cover range without treatment
predictions_grass <- predict(rn_model_nobo, newdata_grass, type = "state", se.fit = TRUE)

# Add predictions and SE to the new data frame
newdata_grass$predicted_state <- predictions_grass$Predicted
newdata_grass$SE <- predictions_grass$SE

# Calculate 95% confidence intervals
newdata_grass$lower_CI <- newdata_grass$predicted_state - 1.96 * newdata_grass$SE
newdata_grass$upper_CI <- newdata_grass$predicted_state + 1.96 * newdata_grass$SE

# Plot the predicted response to grass cover 
ggplot(newdata_grass, aes(x = grass_cover, y = predicted_state)) +
  geom_line(size = 1, color = "black") +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.2, fill = "darkgray") +
  labs(x = "Percent Grass Cover", y = "Relative Abundance per Site", title = "Northern Bobwhite (nobo) -- Predicted Response to Grass Cover") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14),
    title = element_text(size = 14, face = "bold", hjust = .5),
    plot.title = element_text(hjust = 0.4, vjust = -5)
    
    
  )
