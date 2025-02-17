# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  

##########Loading data and formatting for RN mods###########

# Load inbu abundance data
load("data/abundance_data/abundance_24_inbu.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_inbu$y <- ifelse(abund_data_inbu$y > 0 | is.na(abund_data_inbu$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as these RN models can't handle NAs we have in our dataset 
abund_data_inbu$y <- abund_data_inbu$y[, -c(1:55)]
abund_data_inbu$det.covs$day <- abund_data_inbu$det.covs$day[, -c(1:55)]
abund_data_inbu$det.covs$temp <- abund_data_inbu$det.covs$temp[, -c(1:55)]
abund_data_inbu$det.covs$wind_speed <- abund_data_inbu$det.covs$wind_speed[, -c(1:55)]
abund_data_inbu$det.covs$precipitation <- abund_data_inbu$det.covs$precipitation[, -c(1:55)]


##########Fitting Royle-Nichols models###########

# Create the unmarkedFrameOccu object needed to be used for all models made with unmarked
umf_inbu <- unmarkedFrameOccu(
  y = as.matrix(abund_data_inbu$y),
  siteCovs = abund_data_inbu$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_inbu$det.covs$day),
    temp = as.matrix(abund_data_inbu$det.covs$temp),
    wind_speed = as.matrix(abund_data_inbu$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_inbu$det.covs$precipitation)
  )
)


# Fit the first Royle-Nichols model (with site treatment as an abundance covariate, t for treatment)
rn_model_t_inbu <- occuRN( # This is the unmarked function to fit Royle-Nichols (RN) models
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) # Detection model covariates
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), # Abundance model covariates 
  data = umf_inbu, # Specifying the data from the unmarkedFrameOccu object we created
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)


summary(rn_model_t_inbu) # Look at the summary output of the model
# Models have '-1' term in the abundance model formula to remove the intercept, therefore providing mean parameter estimates for each level of the factor (treatment)

# Treatment only model
rn_model_t_only_inbu <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) -1, 
  data = umf_inbu,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_t_only_inbu)

# Vegetation only
rn_model_inbu <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_inbu,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_inbu)

# Fit the null model
rn_model_null_inbu <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ 1, 
  data = umf_inbu,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_inbu)

########AICc comparison and AICc table formatting###########

models_list_inbu <- list(rn_model_t_inbu = rn_model_t_inbu, 
                         rn_model_t_only_inbu = rn_model_t_only_inbu, 
                         rn_model_inbu = rn_model_inbu, 
                         rn_model_null_inbu = rn_model_null_inbu)
model_names_inbu <- c("TV", "T", "V", "Null")
aicc_table_inbu <- aictab(cand.set = models_list_inbu, modnames = model_names_inbu)
print(aicc_table_inbu)

# Load the master AICc table
aicc_table <- read.csv("data/aicc_table.csv")

# Add species name to the table
aicc_table_inbu$species_code <- "inbu"

# Add the AICc table to the master table
aicc_table_master <- rbind(aicc_table, aicc_table_inbu)

# Save the AICc table to csv
write.csv(aicc_table_master, "data/aicc_table.csv", row.names = FALSE)

##########VIF for covariates##########
#Fit another RN mode for rn_model_t_inbu with an intercept to calculate VIF
rn_model_t_inbu_vif <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_inbu,
  K = 25,  
  method = "BFGS" 
)

# Calculate VIF for state and det types
vif_table_state_inbu <- vif(rn_model_t_inbu_vif, type = "state")
vif_table_det_inbu <- vif(rn_model_t_inbu_vif, type = "det")

# Convert VIF tables to data frames
vif_table_state_inbu_df <- as.data.frame(vif_table_state_inbu)
vif_table_state_inbu_df$covariate <- rownames(vif_table_state_inbu_df)
vif_table_state_inbu_df$type <- "abund"
colnames(vif_table_state_inbu_df)[1] <- "vif_value"
vif_table_state_inbu_df$species <- "inbu"

vif_table_det_inbu_df <- as.data.frame(vif_table_det_inbu)
vif_table_det_inbu_df$covariate <- rownames(vif_table_det_inbu_df)
vif_table_det_inbu_df$type <- "det"
colnames(vif_table_det_inbu_df)[1] <- "vif_value"
vif_table_det_inbu_df$species <- "inbu"

# Combine the two VIF dfs
vif_table_inbu <- rbind(vif_table_state_inbu_df, vif_table_det_inbu_df)

# Reset row names to standard numbering
rownames(vif_table_inbu) <- NULL

# Load the master VIF table
vif_table <- read.csv("data/vif_table.csv")

# Add the VIF table to the master table
vif_table_master <- rbind(vif_table, vif_table_inbu)

# Save the VIF table to csv
write.csv(vif_table_master, "data/vif_table.csv", row.names = FALSE)

##########Goodness of fit##########

gof_inbu_t_only <- mb.gof.test(rn_model_t_only_inbu, nsim=100, c.hat.est=TRUE, model.type="royle-nichols")
# Looking more into GoF tests, would like to do k-fold cross-validation too, but just using this for now 
print(gof_inbu_t_only)

gof_inbu <- mb.gof.test(rn_model_t_inbu, nsim=100, c.hat.est=TRUE, model.type="royle-nichols")
print(gof_inbu)

##########Predictions for effect of treatment##########

# Create a new data frame with the mean values of the covariates
newdata_inbu <- data.frame(treatment=levels(umf_inbu@siteCovs$treatment), shrub_cover=mean(umf_inbu@siteCovs$shrub_cover),grass_cover=mean(umf_inbu@siteCovs$grass_cover), forb_cover=mean(umf_inbu@siteCovs$forb_cover),shrub_height=mean(umf_inbu@siteCovs$shrub_height), grass_height=mean(umf_inbu@siteCovs$grass_height))

# Generate predictions for the new data with the predict function
predictions_inbu <- predict(rn_model_t_inbu, newdata_inbu, type = "state", se.fit = TRUE)

# Add predictions and SE to the newdata_inbu dataframe
newdata_inbu$predicted_state <- predictions_inbu$Predicted
newdata_inbu$SE <- predictions_inbu$SE

# Calculate 95% confidence intervals from SE
newdata_inbu$lower_CI <- newdata_inbu$predicted_state - 1.96 * newdata_inbu$SE
newdata_inbu$upper_CI <- newdata_inbu$predicted_state + 1.96 * newdata_inbu$SE

# View the results
print(newdata_inbu)

# Optionally, save results as csv for plotting in another R session
write.csv(newdata_inbu, file = "data/means_abund_parameters/means_treatment_parameters_inbu.csv", row.names = FALSE)

############Predictions for effect of treatment (repeated for t_only)########
newdata_inbu_t_only <- data.frame(treatment=levels(umf_inbu@siteCovs$treatment), shrub_cover=mean(umf_inbu@siteCovs$shrub_cover),grass_cover=mean(umf_inbu@siteCovs$grass_cover), forb_cover=mean(umf_inbu@siteCovs$forb_cover),shrub_height=mean(umf_inbu@siteCovs$shrub_height), grass_height=mean(umf_inbu@siteCovs$grass_height))

predictions_inbu_t_only <- predict(rn_model_t_only_inbu, newdata_inbu_t_only, type = "state", se.fit = TRUE)

newdata_inbu_t_only$predicted_state <- predictions_inbu_t_only$Predicted
newdata_inbu_t_only$SE <- predictions_inbu_t_only$SE
newdata_inbu_t_only$lower_CI <- newdata_inbu_t_only$predicted_state - 1.96 * newdata_inbu_t_only$SE
newdata_inbu_t_only$upper_CI <- newdata_inbu_t_only$predicted_state + 1.96 * newdata_inbu_t_only$SE

print(newdata_inbu_t_only)

write.csv(newdata_inbu_t_only, file = "data/means_abund_parameters/means_treatment_parameters_inbu_t.csv", row.names = FALSE)


##########Predictions for effect of grass cover##########

# Generate new data for grass cover while keeping other covariates constant 
grass_cover_range <- seq(min(umf_inbu@siteCovs$grass_cover), max(umf_inbu@siteCovs$grass_cover), length.out = 100)

# Check the levels of the treatment factor in the original data
treatment_levels <- levels(umf_inbu@siteCovs$treatment)

newdata_grass_cover <- data.frame(
  treatment = factor(treatment_levels, levels = treatment_levels),
  grass_cover = grass_cover_range,
  shrub_cover = mean(umf_inbu@siteCovs$shrub_cover),
  forb_cover = mean(umf_inbu@siteCovs$forb_cover),
  grass_height = mean(umf_inbu@siteCovs$grass_height),
  shrub_height = mean(umf_inbu@siteCovs$shrub_height)
)

predictions_grass_cover <- predict(rn_model_t_inbu, newdata_grass_cover, type = "state")

plot_data_grass_cover <- data.frame(
  grass_cover = newdata_grass_cover$grass_cover,
  treatment = newdata_grass_cover$treatment,
  predicted_state = predictions_grass_cover$Predicted,
  lower_CI = predictions_grass_cover$lower,
  upper_CI = predictions_grass_cover$upper
)

plot_data_grass_cover_mine <- plot_data_grass_cover %>% filter(treatment == "mine")

ggplot(plot_data_grass_cover_mine, aes(x = grass_cover, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#post-hoc models 

#Treatment + Shrub Cover + Grass Cover + Forb Cover + Shrub Height + Grass Height
rn_model_0_inbu <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_inbu, 
  K = 25,  
  method = "BFGS" 
)

#Shrub Cover + Grass Cover + Forb Cover + Shrub Height + Grass Height 
rn_model_1_inbu <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_inbu, 
  K = 25,  
  method = "BFGS" 
)

#Treatment
rn_model_2_inbu <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) -1, 
  data = umf_inbu, 
  K = 25,  
  method = "BFGS" 
)

#Treatment * Shrub Cover * Grass Cover
rn_model_3_inbu <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment)  * scale(shrub_cover) * scale(grass_cover) -1, 
  data = umf_inbu, 
  K = 25,  
  method = "BFGS" 
)

#Treatment + Shrub Cover + Grass Cover
rn_model_4_inbu <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) - 1 + scale(shrub_cover) + scale(grass_cover), 
  data = umf_inbu, 
  K = 25,  
  method = "BFGS" 
)

#Treatment * Grass Cover 
rn_model_5_inbu <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) * scale(grass_cover) -1, 
  data = umf_inbu, 
  K = 25,  
  method = "BFGS" 
)

summary(rn_model_5_inbu)

#Treatment + Grass Cover
rn_model_6_inbu <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) -1  + scale(grass_cover), 
  data = umf_inbu, 
  K = 25,  
  method = "BFGS" 
)

#Treatment * Shrub Cover
rn_model_7_inbu <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) * scale(shrub_cover) -1, 
  data = umf_inbu, 
  K = 25,  
  method = "BFGS" 
)

summary(rn_model_7_inbu)

#Treatment + Shrub Cover
rn_model_8_inbu <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) -1  + scale(shrub_cover), 
  data = umf_inbu, 
  K = 25,  
  method = "BFGS" 
)

#Shrub Cover * Grass Cover 
rn_model_9_inbu <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) * scale(grass_cover), 
  data = umf_inbu, 
  K = 25,  
  method = "BFGS" 
)

#Shrub Cover + Grass Cover 
rn_model_10_inbu <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) + scale(grass_cover), 
  data = umf_inbu, 
  K = 25,  
  method = "BFGS" 
)

#Null
rn_model_11_inbu <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ 1, 
  data = umf_inbu, 
  K = 25,  
  method = "BFGS" 
)

#AICc comparison
models_list_inbu_post_hoc <- list(rn_model_0_inbu, #Treatment + Shrub Cover + Grass Cover + Forb Cover + Shrub Height + Grass Height
                                  rn_model_1_inbu, #Shrub Cover + Grass Cover + Forb Cover + Shrub Height + Grass Height
                                  rn_model_2_inbu, #Treatment
                                  rn_model_3_inbu, #Treatment * Shrub Cover * Grass Cover
                                  rn_model_4_inbu, #Treatment + Shrub Cover + Grass Cover
                                  rn_model_5_inbu, #Treatment * Grass Cover
                                  rn_model_6_inbu, #Treatment + Grass Cover
                                  rn_model_7_inbu, #Treatment * Shrub Cover
                                  rn_model_8_inbu, #Treatment + Shrub Cover
                                  rn_model_9_inbu, #Shrub Cover * Grass Cover
                                  rn_model_10_inbu, #Shrub Cover + Grass Cover
                                  rn_model_11_inbu) #Null

#model names use descriptive 
model_names_inbu_post_hoc <- c("Treatment + Shrub Cover + Grass Cover + Forb Cover + Shrub Height + Grass Height", 
                               "Shrub Cover + Grass Cover + Forb Cover + Shrub Height + Grass Height", 
                               "Treatment", 
                               "Treatment * Shrub Cover * Grass Cover", 
                               "Treatment + Shrub Cover + Grass Cover", 
                               "Treatment * Grass Cover", 
                               "Treatment + Grass Cover", 
                               "Treatment * Shrub Cover", 
                               "Treatment + Shrub Cover", 
                               "Shrub Cover * Grass Cover", 
                               "Shrub Cover + Grass Cover", 
                               "Null")


aicc_table_inbu_post_hoc <- aictab(cand.set = models_list_inbu_post_hoc, modnames = model_names_inbu_post_hoc)

print(aicc_table_inbu_post_hoc)

# Load the master AICc table
aicc_table <- read.csv("data/aicc_table_posthoc.csv")

# Add species name to the table
aicc_table_inbu_post_hoc$species_code <- "Indigo Bunting"

# Add the AICc table to the master table
aicc_table_master <- rbind(aicc_table, aicc_table_inbu_post_hoc)

# Save the AICc table to csv
write.csv(aicc_table_master, "data/aicc_table_posthoc.csv", row.names = FALSE)

#vif for top model
vif_table_state_inbu <- vif(rn_model_6_inbu, type = "state")
vif_table_det_inbu <- vif(rn_model_6_inbu, type = "det")

# Convert VIF tables to data frames
vif_table_state_inbu_df <- as.data.frame(vif_table_state_inbu)
vif_table_state_inbu_df$covariate <- rownames(vif_table_state_inbu_df)
vif_table_state_inbu_df$type <- "abund"
colnames(vif_table_state_inbu_df)[1] <- "vif_value"
vif_table_state_inbu_df$species <- "inbu"

vif_table_det_inbu_df <- as.data.frame(vif_table_det_inbu)
vif_table_det_inbu_df$covariate <- rownames(vif_table_det_inbu_df)
vif_table_det_inbu_df$type <- "det"
colnames(vif_table_det_inbu_df)[1] <- "vif_value"
vif_table_det_inbu_df$species <- "inbu"

# Combine the two VIF dfs
vif_table_inbu <- rbind(vif_table_state_inbu_df, vif_table_det_inbu_df)

# Reset row names to standard numbering
rownames(vif_table_inbu) <- NULL