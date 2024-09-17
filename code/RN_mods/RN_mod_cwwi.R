# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  

# Load cwwi abundance data
load("data/abundance_data/abundance_24_cwwi.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_cwwi$y <- ifelse(abund_data_cwwi$y > 0 | is.na(abund_data_cwwi$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as RN model can't handle NAs
abund_data_cwwi$y <- abund_data_cwwi$y[, -c(1:55)]
abund_data_cwwi$det.covs$day <- abund_data_cwwi$det.covs$day[, -c(1:55)]
abund_data_cwwi$det.covs$temp <- abund_data_cwwi$det.covs$temp[, -c(1:55)]
abund_data_cwwi$det.covs$wind_speed <- abund_data_cwwi$det.covs$wind_speed[, -c(1:55)]
abund_data_cwwi$det.covs$precipitation <- abund_data_cwwi$det.covs$precipitation[, -c(1:55)]


# Create the unmarkedFrameOccu object
umf_cwwi <- unmarkedFrameOccu(
  y = as.matrix(abund_data_cwwi$y),
  siteCovs = abund_data_cwwi$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_cwwi$det.covs$day),
    temp = as.matrix(abund_data_cwwi$det.covs$temp),
    wind_speed = as.matrix(abund_data_cwwi$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_cwwi$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model
rn_model_t_cwwi <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_cwwi,
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)
summary(rn_model_t_cwwi)

# Fit the second Royle-Nichols model
rn_model_cwwi <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_cwwi,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_cwwi)

# Fit the null model
rn_model_null_cwwi <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) ~ 1, 
  data = umf_cwwi,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_cwwi)

# AICc comparison
models_list_cwwi <- list(rn_model_t_cwwi = rn_model_t_cwwi, rn_model_cwwi = rn_model_cwwi, rn_model_null_cwwi = rn_model_null_cwwi)
model_names_cwwi <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_cwwi <- aictab(cand.set = models_list_cwwi, modnames = model_names_cwwi)
print(aicc_table_cwwi)

# Function returning fit-statistics
fitstats_cwwi <- function(rn_model_t_cwwi) {
  observed <- getY(rn_model_t_cwwi@data)
  expected <- fitted(rn_model_t_cwwi)
  resids <- residuals(rn_model_t_cwwi)
  sse <- sum(resids^2, na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
  fit_out <- c(cSSE = sse, Chisq = chisq, freemanTukey = freeTuke)
  return(fit_out)
}

#parametric bootstrapping 
# Assuming 'mod' is your fitted model
pb_cwwi <- parboot(rn_model_t_cwwi, fitstats_cwwi, nsim=100)
print(pb_cwwi)

# Create a sequence of shrub cover values for predictions
shrub_vals <- seq(min(cwwi_abund_data_24$abund.covs$shrub_cover), 
                  max(cwwi_abund_data_24$abund.covs$shrub_cover), 
                  length.out = 100)

# Create a new dataframe for predictions, using the valid 'treatment' level 'mine'
new_data <- data.frame(
  shrub_cover = shrub_vals,
  grass_cover = mean(cwwi_abund_data_24$abund.covs$grass_cover),
  forb_cover = mean(cwwi_abund_data_24$abund.covs$forb_cover),
  shrub_height = mean(cwwi_abund_data_24$abund.covs$shrub_height),
  grass_height = mean(cwwi_abund_data_24$abund.covs$grass_height),
  treatment = factor("mine", levels = levels(cwwi_abund_data_24$abund.covs$treatment))
)

# Use the predict function to predict abundance based on shrub cover
predictions <- predict(rn_model, type = "state", newdata = new_data)

# Combine predicted values with the shrub cover values in a new dataframe for plotting
plot_data <- data.frame(
  shrub_cover = shrub_vals,
  predicted_abundance = predictions$Predicted,
  lower_CI = predictions$lower,
  upper_CI = predictions$upper
)

# Plot the predicted abundance against shrub cover using ggplot2
ggplot(plot_data, aes(x = shrub_cover, y = predicted_abundance)) +
  geom_line(color = "#721121", size = 1.5) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.2, fill = "#FFCF99") +
  labs(
    x = "Shrub Cover Percent",
    y = "Relative Abundance per site",
    title = "Chuck-will's-widow -- Predicted Response to Shrub Cover" 
  ) +
  theme_classic() + 
  theme(
    plot.title = element_text(vjust = -5, hjust = 0.5, size = 18),  
    axis.title.x = element_text(vjust = -.5, size = 14),             
    axis.title.y = element_text(vjust = 1.5, size = 14)             
  )


