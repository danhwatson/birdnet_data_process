# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  
library(pROC)


##########Loading data and formatting for RN mods###########

# Load wevi abundance data
load("data/abundance_data/abundance_24_wevi.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_wevi$y <- ifelse(abund_data_wevi$y > 0 | is.na(abund_data_wevi$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as these RN models can't handle NAs we have in our dataset 
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

rn_model_t_wevi_forb_height <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height) + scale(forb_height), 
  data = umf_wevi, 
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_t_wevi_forb_height)


rn_model_t_wevi_temp <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation)
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_wevi, 
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_t_wevi_temp)

# Fit the first Royle-Nichols model (with site treatment as an abundance covariate, t for treatment)
rn_model_t_wevi <- occuRN( # This is the unmarked function to fit Royle-Nichols (RN) models
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) # Detection model covariates
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), # Abundance model covariates 
  data = umf_wevi, # Specifying the data from the unmarkedFrameOccu object we created
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)

AICc(rn_model_t_wevi)
AICc(rn_model_t_wevi_forb_height)
AICc(rn_model_t_wevi_temp)


summary(rn_model_t_wevi) # Look at the summary output of the model
# Models have '-1' term in the abundance model formula to remove the intercept, therefore providing mean parameter estimates for each level of the factor (treatment)

saveRDS(rn_model_t_wevi, file = "data/rn_model_files/rn_model_t_wevi.rds") 
# Optional step to save model object to load in another R session

# Treatment only model
rn_model_t_only_wevi <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) -1, 
  data = umf_wevi,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_t_only_wevi)

# Fit the second Royle-Nichols model (without site treatment as an abundance covariate)
rn_model_wevi <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_wevi,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_wevi)

# Fit the null model
rn_model_null_wevi <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ 1, 
  data = umf_wevi,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_wevi)

########AICc comparison and AICc table formatting###########

models_list_wevi <- list(rn_model_t_wevi = rn_model_t_wevi, 
                         rn_model_t_only_wevi = rn_model_t_only_wevi, 
                         rn_model_wevi = rn_model_wevi, 
                         rn_model_null_wevi = rn_model_null_wevi)
model_names_wevi <- c("TV", "T", "V", "Null")
aicc_table_wevi <- aictab(cand.set = models_list_wevi, modnames = model_names_wevi)
print(aicc_table_wevi)

# Load the master AICc table
aicc_table <- read.csv("data/aicc_table.csv")

# Add species name to the table
aicc_table_wevi$species_code <- "wevi"

# Add the AICc table to the master table
aicc_table_master <- rbind(aicc_table, aicc_table_wevi)

# Save the AICc table to csv
write.csv(aicc_table_master, "data/aicc_table.csv", row.names = FALSE)

##########VIF for covariates##########
#Fit another RN mode for rn_model_t_wevi with an intercept to calculate VIF
rn_model_t_wevi_vif <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_wevi,
  K = 25,  
  method = "BFGS" 
)

# Calculate VIF for state and det types
vif_table_state_wevi <- vif(rn_model_t_wevi_vif, type = "state")
vif_table_det_wevi <- vif(rn_model_t_wevi_vif, type = "det")

# Convert VIF tables to data frames
vif_table_state_wevi_df <- as.data.frame(vif_table_state_wevi)
vif_table_state_wevi_df$covariate <- rownames(vif_table_state_wevi_df)
vif_table_state_wevi_df$type <- "abund"
colnames(vif_table_state_wevi_df)[1] <- "vif_value"
vif_table_state_wevi_df$species <- "wevi"

vif_table_det_wevi_df <- as.data.frame(vif_table_det_wevi)
vif_table_det_wevi_df$covariate <- rownames(vif_table_det_wevi_df)
vif_table_det_wevi_df$type <- "det"
colnames(vif_table_det_wevi_df)[1] <- "vif_value"
vif_table_det_wevi_df$species <- "wevi"

# Combine the two VIF dfs
vif_table_wevi <- rbind(vif_table_state_wevi_df, vif_table_det_wevi_df)

# Reset row names to standard numbering
rownames(vif_table_wevi) <- NULL

# Load the master VIF table
vif_table <- read.csv("data/vif_table.csv")

# Add the VIF table to the master table
vif_table_master <- rbind(vif_table, vif_table_wevi)

# Save the VIF table to csv
write.csv(vif_table_master, "data/vif_table.csv", row.names = FALSE)

##########Goodness of fit##########

gof_wevi <- mb.gof.test(rn_model_t_wevi, nsim=500, c.hat.est=TRUE, model.type="royle-nichols")
# Looking more into GoF tests, would like to do k-fold cross-validation too, but just using this for now 
print(gof_wevi)

gof.

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
write.csv(newdata_wevi, file = "data/means_abund_parameters/means_treatment_parameters_wevi.csv", row.names = FALSE)

############Predictions for effect of treatment (repeated for t_only)########
newdata_wevi_t_only <- data.frame(treatment=levels(umf_wevi@siteCovs$treatment), shrub_cover=mean(umf_wevi@siteCovs$shrub_cover),grass_cover=mean(umf_wevi@siteCovs$grass_cover), forb_cover=mean(umf_wevi@siteCovs$forb_cover),shrub_height=mean(umf_wevi@siteCovs$shrub_height), grass_height=mean(umf_wevi@siteCovs$grass_height))

predictions_wevi_t_only <- predict(rn_model_t_only_wevi, newdata_wevi_t_only, type = "state", se.fit = TRUE)

newdata_wevi_t_only$predicted_state <- predictions_wevi_t_only$Predicted
newdata_wevi_t_only$SE <- predictions_wevi_t_only$SE
newdata_wevi_t_only$lower_CI <- newdata_wevi_t_only$predicted_state - 1.96 * newdata_wevi_t_only$SE
newdata_wevi_t_only$upper_CI <- newdata_wevi_t_only$predicted_state + 1.96 * newdata_wevi_t_only$SE

print(newdata_wevi_t_only)

write.csv(newdata_wevi_t_only, file = "data/means_abund_parameters/means_treatment_parameters_wevi_t.csv", row.names = FALSE)


##########Predictions for effect of grass cover##########

# Generate new data for grass cover while keeping other covariates constant 
grass_cover_range <- seq(min(umf_wevi@siteCovs$grass_cover), max(umf_wevi@siteCovs$grass_cover), length.out = 100)

# Check the levels of the treatment factor in the original data
treatment_levels <- levels(umf_wevi@siteCovs$treatment)

newdata_grass_cover <- data.frame(
  treatment = factor(treatment_levels, levels = treatment_levels),
  grass_cover = grass_cover_range,
  shrub_cover = mean(umf_wevi@siteCovs$shrub_cover),
  forb_cover = mean(umf_wevi@siteCovs$forb_cover),
  grass_height = mean(umf_wevi@siteCovs$grass_height),
  shrub_height = mean(umf_wevi@siteCovs$shrub_height)
)

predictions_grass_cover <- predict(rn_model_t_wevi, newdata_grass_cover, type = "state")

plot_data_grass_cover <- data.frame(
  grass_cover = newdata_grass_cover$grass_cover,
  treatment = newdata_grass_cover$treatment,
  predicted_state = predictions_grass_cover$Predicted,
  lower_CI = predictions_grass_cover$lower,
  upper_CI = predictions_grass_cover$upper
)

ggplot(plot_data_grass_cover, aes(x = grass_cover, y = predicted_state, color = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill = treatment), alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

plot_data_grass_cover_mine <- plot_data_grass_cover %>% filter(treatment == "mine")
plot_data_grass_cover_timber <- plot_data_grass_cover %>% filter(treatment == "timber")
plot_data_grass_cover_rx_fire_young <- plot_data_grass_cover %>% filter(treatment == "rx_fire_young")
plot_data_grass_cover_rx_fire_sec_growth <- plot_data_grass_cover %>% filter(treatment == "rx_fire_sec_growth")

ggplot(plot_data_grass_cover_mine, aes(x = grass_cover, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/mine
dir.create("figures/predicted_abundance/mine/wevi")
ggsave("figures/predicted_abundance/mine/wevi/mine_grass_cover_wevi.png", width = 4, height = 4)