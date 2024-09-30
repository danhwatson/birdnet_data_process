# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  


##########Loading data and formatting for RN mods###########

# Load coni abundance data
load("data/abundance_data/abundance_24_coni.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_coni$y <- ifelse(abund_data_coni$y > 0 | is.na(abund_data_coni$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as these RN models can't handle NAs we have in our dataset 
abund_data_coni$y <- abund_data_coni$y[, -c(1:55)]
abund_data_coni$det.covs$day <- abund_data_coni$det.covs$day[, -c(1:55)]
abund_data_coni$det.covs$temp <- abund_data_coni$det.covs$temp[, -c(1:55)]
abund_data_coni$det.covs$wind_speed <- abund_data_coni$det.covs$wind_speed[, -c(1:55)]
abund_data_coni$det.covs$precipitation <- abund_data_coni$det.covs$precipitation[, -c(1:55)]


##########Fitting Royle-Nichols models###########

# Create the unmarkedFrameOccu object needed to be used for all models made with unmarked
umf_coni <- unmarkedFrameOccu(
  y = as.matrix(abund_data_coni$y),
  siteCovs = abund_data_coni$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_coni$det.covs$day),
    temp = as.matrix(abund_data_coni$det.covs$temp),
    wind_speed = as.matrix(abund_data_coni$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_coni$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model (with site treatment as an abundance covariate, t for treatment)
rn_model_t_coni <- occuRN( # This is the unmarked function to fit Royle-Nichols (RN) models
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) # Detection model covariates
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), # Abundance model covariates 
  data = umf_coni, # Specifying the data from the unmarkedFrameOccu object we created
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)

summary(rn_model_t_coni) # Look at the summary output of the model
# Models have '-1' term in the abundance model formula to remove the intercept, therefore providing mean parameter estimates for each level of the factor (treatment)

#saveRDS(rn_model_t_coni, file = "data/rn_model_t_coni.rds") 
# Optional step to save model object to load in another R session


# Fit the second Royle-Nichols model (without site treatment as an abundance covariate)
rn_model_coni <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_coni,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_coni)

# Fit the null model
rn_model_null_coni <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
  ~ 1, 
  data = umf_coni,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_coni)

########AICc comparison and AICc table formatting###########

# Load the master AICc table
aicc_table <- read.csv("data/aicc_table.csv")

models_list_coni <- list(rn_model_t_coni = rn_model_t_coni, rn_model_coni = rn_model_coni, rn_model_null_coni = rn_model_null_coni)
model_names_coni <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_coni <- aictab(cand.set = models_list_coni, modnames = model_names_coni)
print(aicc_table_coni)

# Add species name to the table
aicc_table_coni$species_code <- "CONI"

# Add the AICc table to the master table
aicc_table_master <- rbind(aicc_table, aicc_table_coni)

# Save the AICc table to csv
write.csv(aicc_table_master, "data/aicc_table.csv", row.names = FALSE)

##########Goodness of fit##########

gof <- mb.gof.test(rn_model_t_coni, nsim=100, c.hat.est=TRUE, model.type="RN")
# Looking more into GoF tests, would like to do k-fold cross-validation too, but just using this for now 
print(gof)

##########Predictions for effect of treatment##########

# Create a new data frame with the mean values of the covariates
newdata_coni <- data.frame(treatment=levels(umf_coni@siteCovs$treatment), shrub_cover=mean(umf_coni@siteCovs$shrub_cover),grass_cover=mean(umf_coni@siteCovs$grass_cover), forb_cover=mean(umf_coni@siteCovs$forb_cover),shrub_height=mean(umf_coni@siteCovs$shrub_height), grass_height=mean(umf_coni@siteCovs$grass_height))

# Generate predictions for the new data with the predict function
predictions_coni <- predict(rn_model_t_coni, newdata_coni, type = "state", se.fit = TRUE)

# Add predictions and SE to the newdata_coni dataframe
newdata_coni$predicted_state <- predictions_coni$Predicted
newdata_coni$SE <- predictions_coni$SE

# Calculate 95% confidence intervals from SE
newdata_coni$lower_CI <- newdata_coni$predicted_state - 1.96 * newdata_coni$SE
newdata_coni$upper_CI <- newdata_coni$predicted_state + 1.96 * newdata_coni$SE

# View the results
print(newdata_coni)

# Optionally, save results as csv for plotting in another R session
#write.csv(newdata_coni, file = "data/means_treatment_parameters_coni.csv", row.names = FALSE)

##########Predictions for effect of grass cover##########

# Generate new data for grass cover while keeping other covariates constant
grass_range <- seq(min(umf_coni@siteCovs$grass_cover), max(umf_coni@siteCovs$grass_cover), length.out = 100)

# Check the levels of the treatment factor in the original data
treatment_levels <- levels(umf_coni@siteCovs$treatment)

newdata_grass <- data.frame(
  treatment = factor(treatment_levels, levels = treatment_levels),
  grass_cover = grass_range,
  shrub_cover = mean(umf_coni@siteCovs$shrub_cover),
  forb_cover = mean(umf_coni@siteCovs$forb_cover),
  shrub_height = mean(umf_coni@siteCovs$shrub_height),
  grass_height = mean(umf_coni@siteCovs$grass_height)
)

# Predict response in abundance for grass cover range across all treatment levels
predictions_grass <- predict(rn_model_t_coni, newdata_grass, type = "state")

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
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill = treatment), alpha = 0.2) +
  labs(x = " Percent Grass Cover", y = "Predicted Relative Abundance", title = "") +
  theme_classic() 

# Filter the plot_data_grass dataframe for only the 'mine' treatment
plot_data_grass_mine <- plot_data_grass %>% filter(treatment == "mine")

# Plot the predictive response to grass cover for the 'mine' treatment only
ggplot(plot_data_grass_mine, aes(x = grass_cover, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "Percent Grass Cover", y = "Predicted Relative Abundance", title = "") +
  theme_classic()


