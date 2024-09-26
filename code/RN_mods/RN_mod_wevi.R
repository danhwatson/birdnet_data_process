# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  


##########Loading data and formatting for RN mods###########

# Load wevi abundance data
load("data/abundance_data/abundance_24_wevi.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_wevi$y <- ifelse(abund_data_wevi$y > 0 | is.na(abund_data_wevi$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as RN model can't handle NAs
abund_data_wevi$y <- abund_data_wevi$y[, -c(1:55)]
abund_data_wevi$det.covs$day <- abund_data_wevi$det.covs$day[, -c(1:55)]
abund_data_wevi$det.covs$temp <- abund_data_wevi$det.covs$temp[, -c(1:55)]
abund_data_wevi$det.covs$wind_speed <- abund_data_wevi$det.covs$wind_speed[, -c(1:55)]
abund_data_wevi$det.covs$precipitation <- abund_data_wevi$det.covs$precipitation[, -c(1:55)]


##########Fitting Royle-Nichols models###########

# Create the unmarkedFrameOccu object needed to be used for all models made with unmarked
umf_wevi <- unmarkedFrameOccu(
  y = as.matrix(abund_data_wevi$y),
  siteCovs = abund_data_wevi$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_wevi$det.covs$day),
    temp = as.matrix(abund_data_wevi$det.covs$temp),
    wind_speed = as.matrix(abund_data_wevi$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_wevi$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model (with site treatment as an abundance covariate, t for treatment)
rn_model_t_wevi <- occuRN( # This is the unmarked function to fit Royle-Nichols (RN) models
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) # Detection model covariates
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), # Abundance model covariates 
  data = umf_wevi, # Specifying the data from the unmarkedFrameOccu object we created
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)

summary(rn_model_t_wevi) # Look at the summary output of the model
# Models have '-1' term in the abundance model formula to remove the intercept, therefore providing mean parameter estimates for each level of the factor (treatment)

#saveRDS(rn_model_t_wevi, file = "data/rn_model_t_wevi.rds") 
# Optional step to save model object to load in another R session


# Fit the second Royle-Nichols model (without site treatment as an abundance covariate)
rn_model_wevi <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_wevi,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_wevi)

# Fit the null model
rn_model_null_wevi <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
  ~ 1, 
  data = umf_wevi,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_wevi)

########AICc comparison###########

# Load the master AICc table
aicc_table <- read.csv("data/aicc_table.csv")

models_list_wevi <- list(rn_model_t_wevi = rn_model_t_wevi, rn_model_wevi = rn_model_wevi, rn_model_null_wevi = rn_model_null_wevi)
model_names_wevi <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_wevi <- aictab(cand.set = models_list_wevi, modnames = model_names_wevi)
print(aicc_table_wevi)

# Add species name to the table
aicc_table_wevi$species_code <- "WEVI"

# Add the AICc table to the master table
aicc_table_master <- rbind(aicc_table, aicc_table_wevi)

# Save the AICc table to csv
write.csv(aicc_table_master, "data/aicc_table.csv", row.names = FALSE)


##########Goodness of fit##########

# Function returning fit-statistics
#fitstats_wevi <- function(rn_model_t_wevi) {
#  observed <- getY(rn_model_t_wevi@data)
#  expected <- fitted(rn_model_t_wevi)
#  resids <- residuals(rn_model_t_wevi)
#  sse <- sum(resids^2, na.rm=TRUE)
#  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
#  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
#  fit_out <- c(cSSE = sse, Chisq = chisq, freemanTukey = freeTuke)
#  return(fit_out)
#}


# Need to look more into GoF tests
#gof <- mb.gof.test(rn_model_t_wevi, nsim=100, c.hat.est=TRUE, model.type="RN")

#parametric bootstrapping 
# Assuming 'mod' is your fitted model
#pb_wevi <- parboot(rn_model_t_wevi, fitstats_wevi, nsim=100)
#print(pb_wevi)


##########Predictions for effect of treatment##########

# Create a new data frame with the mean values of the covariates
newdata_wevi <- data.frame(treatment=levels(umf_wevi@siteCovs$treatment), shrub_cover=mean(umf_wevi@siteCovs$shrub_cover),grass_cover=mean(umf_wevi@siteCovs$grass_cover), forb_cover=mean(umf_wevi@siteCovs$forb_cover),shrub_height=mean(umf_wevi@siteCovs$shrub_height), grass_height=mean(umf_wevi@siteCovs$grass_height))

# Generate predictions for the new data with the predict function
predictions_wevi <- predict(rn_model_t_wevi, newdata_wevi, type = "state", se.fit = TRUE)

# Add predictions and SE to the newdata_wevi dataframe
newdata_wevi$predicted_state <- predictions_wevi$Predicted
newdata_wevi$SE <- predictions_wevi$SE

# Calculate 95% confidence intervals from SE
newdata_wevi$lower_CI <- newdata_wevi$predicted_state - 1.96 * newdata_wevi$SE
newdata_wevi$upper_CI <- newdata_wevi$predicted_state + 1.96 * newdata_wevi$SE

# View the results
print(newdata_wevi)

# Optionally, save results as csv for plotting in another R session
write.csv(newdata_wevi, file = "data/means_treatment_parameters_wevi.csv", row.names = FALSE)


##########Predictions for effect of grass cover##########

# Generate new data for grass cover while keeping other covariates constant
grass_range <- seq(min(umf_wevi@siteCovs$grass_cover), max(umf_wevi@siteCovs$grass_cover), length.out = 100)

# Check the levels of the treatment factor in the original data
treatment_levels <- levels(umf_wevi@siteCovs$treatment)

newdata_grass <- data.frame(
  treatment = factor(treatment_levels, levels = treatment_levels),
  grass_cover = grass_range,
  shrub_cover = mean(umf_wevi@siteCovs$shrub_cover),
  forb_cover = mean(umf_wevi@siteCovs$forb_cover),
  shrub_height = mean(umf_wevi@siteCovs$shrub_height),
  grass_height = mean(umf_wevi@siteCovs$grass_height)
)

# Predict state (occupancy or abundance) for grass  cover range across all treatment levels
predictions_grass <- predict(rn_model_t_wevi, newdata_grass, type = "state")

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
  labs(x = "Grass Cover", y = "Predicted Abundance", title = "Northern Bobwhite (wevi) -- Predicted Response to Grass Cover across Treatments") +
  theme_minimal()

# Generate new data for grass cover without the treatment variable
grass_range <- seq(min(umf_wevi@siteCovs$grass_cover), max(umf_wevi@siteCovs$grass_cover), length.out = 100)

newdata_grass <- data.frame(
  grass_cover = grass_range,
  shrub_cover = mean(umf_wevi@siteCovs$shrub_cover),
  forb_cover = mean(umf_wevi@siteCovs$forb_cover),
  shrub_height = mean(umf_wevi@siteCovs$shrub_height),
  grass_height = mean(umf_wevi@siteCovs$grass_height)
)

# Predict state (occupancy or abundance) for grass cover range without treatment
predictions_grass <- predict(rn_model_wevi, newdata_grass, type = "state", se.fit = TRUE)

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
  labs(x = "Percent Grass Cover", y = "Relative Abundance per Site", title = "Northern Bobwhite (wevi) -- Predicted Response to Grass Cover") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14),
    title = element_text(size = 14, face = "bold", hjust = .5),
    plot.title = element_text(hjust = 0.4, vjust = -5)
    
    
  )

