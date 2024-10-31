# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  
library(pROC)


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

rn_model_t_coni_forb_height <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height) + scale(forb_height), 
  data = umf_coni, 
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_t_coni_forb_height)


rn_model_t_coni_temp <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation)
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_coni, 
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_t_coni_temp)

# Fit the first Royle-Nichols model (with site treatment as an abundance covariate, t for treatment)
rn_model_t_coni <- occuRN( # This is the unmarked function to fit Royle-Nichols (RN) models
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) # Detection model covariates
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), # Abundance model covariates 
  data = umf_coni, # Specifying the data from the unmarkedFrameOccu object we created
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)

AICc(rn_model_t_coni)
AICc(rn_model_t_coni_forb_height)
AICc(rn_model_t_coni_temp)


summary(rn_model_t_coni) # Look at the summary output of the model
# Models have '-1' term in the abundance model formula to remove the intercept, therefore providing mean parameter estimates for each level of the factor (treatment)

saveRDS(rn_model_t_coni, file = "data/rn_model_files/rn_model_t_coni.rds") 
# Optional step to save model object to load in another R session

# Treatment only model
rn_model_t_only_coni <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) -1, 
  data = umf_coni,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_t_only_coni)

# Fit the second Royle-Nichols model (without site treatment as an abundance covariate)
rn_model_coni <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_coni,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_coni)

# Fit the null model
rn_model_null_coni <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ 1, 
  data = umf_coni,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_coni)

########AICc comparison and AICc table formatting###########

models_list_coni <- list(rn_model_t_coni = rn_model_t_coni, 
                         rn_model_t_only_coni = rn_model_t_only_coni, 
                         rn_model_coni = rn_model_coni, 
                         rn_model_null_coni = rn_model_null_coni)
model_names_coni <- c("TV", "T", "V", "Null")
aicc_table_coni <- aictab(cand.set = models_list_coni, modnames = model_names_coni)
print(aicc_table_coni)

# Load the master AICc table
aicc_table <- read.csv("data/aicc_table.csv")

# Add species name to the table
aicc_table_coni$species_code <- "coni"

# Add the AICc table to the master table
aicc_table_master <- rbind(aicc_table, aicc_table_coni)

# Save the AICc table to csv
write.csv(aicc_table_master, "data/aicc_table.csv", row.names = FALSE)

##########VIF for covariates##########
#Fit another RN mode for rn_model_t_coni with an intercept to calculate VIF
rn_model_t_coni_vif <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_coni,
  K = 25,  
  method = "BFGS" 
)

# Calculate VIF for state and det types
vif_table_state_coni <- vif(rn_model_t_coni_vif, type = "state")
vif_table_det_coni <- vif(rn_model_t_coni_vif, type = "det")

# Convert VIF tables to data frames
vif_table_state_coni_df <- as.data.frame(vif_table_state_coni)
vif_table_state_coni_df$covariate <- rownames(vif_table_state_coni_df)
vif_table_state_coni_df$type <- "abund"
colnames(vif_table_state_coni_df)[1] <- "vif_value"
vif_table_state_coni_df$species <- "coni"

vif_table_det_coni_df <- as.data.frame(vif_table_det_coni)
vif_table_det_coni_df$covariate <- rownames(vif_table_det_coni_df)
vif_table_det_coni_df$type <- "det"
colnames(vif_table_det_coni_df)[1] <- "vif_value"
vif_table_det_coni_df$species <- "coni"

# Combine the two VIF dfs
vif_table_coni <- rbind(vif_table_state_coni_df, vif_table_det_coni_df)

# Reset row names to standard numbering
rownames(vif_table_coni) <- NULL

# Load the master VIF table
vif_table <- read.csv("data/vif_table.csv")

# Add the VIF table to the master table
vif_table_master <- rbind(vif_table, vif_table_coni)

# Save the VIF table to csv
write.csv(vif_table_master, "data/vif_table.csv", row.names = FALSE)

##########Goodness of fit##########

gof_coni <- mb.gof.test(rn_model_t_coni, nsim=100, c.hat.est=TRUE, model.type="royle-nichols")
# Looking more into GoF tests, would like to do k-fold cross-validation too, but just using this for now 
print(gof_coni)

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
write.csv(newdata_coni, file = "data/means_abund_parameters/means_treatment_parameters_coni.csv", row.names = FALSE)

############Predictions for effect of treatment (repeated for t_only)########
newdata_coni_t_only <- data.frame(treatment=levels(umf_coni@siteCovs$treatment), shrub_cover=mean(umf_coni@siteCovs$shrub_cover),grass_cover=mean(umf_coni@siteCovs$grass_cover), forb_cover=mean(umf_coni@siteCovs$forb_cover),shrub_height=mean(umf_coni@siteCovs$shrub_height), grass_height=mean(umf_coni@siteCovs$grass_height))

predictions_coni_t_only <- predict(rn_model_t_only_coni, newdata_coni_t_only, type = "state", se.fit = TRUE)

newdata_coni_t_only$predicted_state <- predictions_coni_t_only$Predicted
newdata_coni_t_only$SE <- predictions_coni_t_only$SE
newdata_coni_t_only$lower_CI <- newdata_coni_t_only$predicted_state - 1.96 * newdata_coni_t_only$SE
newdata_coni_t_only$upper_CI <- newdata_coni_t_only$predicted_state + 1.96 * newdata_coni_t_only$SE

print(newdata_coni_t_only)

write.csv(newdata_coni_t_only, file = "data/means_abund_parameters/means_treatment_parameters_coni_t.csv", row.names = FALSE)


##########Predictions for effect of grass cover##########

# Generate new data for grass cover while keeping other covariates constant 
grass_cover_range <- seq(min(umf_coni@siteCovs$grass_cover), max(umf_coni@siteCovs$grass_cover), length.out = 100)

# Check the levels of the treatment factor in the original data
treatment_levels <- levels(umf_coni@siteCovs$treatment)

newdata_grass_cover <- data.frame(
  treatment = factor(treatment_levels, levels = treatment_levels),
  grass_cover = grass_cover_range,
  shrub_cover = mean(umf_coni@siteCovs$shrub_cover),
  forb_cover = mean(umf_coni@siteCovs$forb_cover),
  grass_height = mean(umf_coni@siteCovs$grass_height),
  shrub_height = mean(umf_coni@siteCovs$shrub_height)
)

predictions_grass_cover <- predict(rn_model_t_coni, newdata_grass_cover, type = "state")

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
dir.create("figures/predicted_abundance/mine/coni")
ggsave("figures/predicted_abundance/mine/coni/mine_grass_cover_coni.png", width = 4, height = 4)