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

# Remove the first 55 visits (columns) from "y" and detection covariates as these RN models can't handle NAs we have in our dataset 
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
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) # Detection model covariates
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), # Abundance model covariates 
  data = umf_nobo, # Specifying the data from the unmarkedFrameOccu object we created
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)

summary(rn_model_t_nobo) # Look at the summary output of the model
# Models have '-1' term in the abundance model formula to remove the intercept, therefore providing mean parameter estimates for each level of the factor (treatment)

# Treatment only model
rn_model_t_only_nobo <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) -1, 
  data = umf_nobo,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_t_only_nobo)

# Vegetation only
rn_model_nobo <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_nobo,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_nobo)

# Fit the null model
rn_model_null_nobo <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ 1, 
  data = umf_nobo,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_nobo)

########AICc comparison and AICc table formatting###########

models_list_nobo <- list(rn_model_t_nobo = rn_model_t_nobo, 
                         rn_model_t_only_nobo = rn_model_t_only_nobo, 
                         rn_model_nobo = rn_model_nobo, 
                         rn_model_null_nobo = rn_model_null_nobo)
model_names_nobo <- c("TV", "T", "V", "Null")
aicc_table_nobo <- aictab(cand.set = models_list_nobo, modnames = model_names_nobo)
print(aicc_table_nobo)

# Load the master AICc table
aicc_table <- read.csv("data/aicc_table.csv")

# Add species name to the table
aicc_table_nobo$species_code <- "nobo"

# Add the AICc table to the master table
aicc_table_master <- rbind(aicc_table, aicc_table_nobo)

# Save the AICc table to csv
write.csv(aicc_table_master, "data/aicc_table.csv", row.names = FALSE)

##########VIF for covariates##########
#Fit another RN mode for rn_model_t_nobo with an intercept to calculate VIF
rn_model_t_nobo_vif <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_nobo,
  K = 25,  
  method = "BFGS" 
)

# Calculate VIF for state and det types
vif_table_state_nobo <- vif(rn_model_t_nobo_vif, type = "state")
vif_table_det_nobo <- vif(rn_model_t_nobo_vif, type = "det")

# Convert VIF tables to data frames
vif_table_state_nobo_df <- as.data.frame(vif_table_state_nobo)
vif_table_state_nobo_df$covariate <- rownames(vif_table_state_nobo_df)
vif_table_state_nobo_df$type <- "abund"
colnames(vif_table_state_nobo_df)[1] <- "vif_value"
vif_table_state_nobo_df$species <- "nobo"

vif_table_det_nobo_df <- as.data.frame(vif_table_det_nobo)
vif_table_det_nobo_df$covariate <- rownames(vif_table_det_nobo_df)
vif_table_det_nobo_df$type <- "det"
colnames(vif_table_det_nobo_df)[1] <- "vif_value"
vif_table_det_nobo_df$species <- "nobo"

# Combine the two VIF dfs
vif_table_nobo <- rbind(vif_table_state_nobo_df, vif_table_det_nobo_df)

# Reset row names to standard numbering
rownames(vif_table_nobo) <- NULL

# Load the master VIF table
vif_table <- read.csv("data/vif_table.csv")

# Add the VIF table to the master table
vif_table_master <- rbind(vif_table, vif_table_nobo)

# Save the VIF table to csv
write.csv(vif_table_master, "data/vif_table.csv", row.names = FALSE)

##########Goodness of fit##########

gof_nobo <- mb.gof.test(rn_model_t_only_nobo, nsim=500, c.hat.est=TRUE, model.type="royle-nichols")
# Looking more into GoF tests, would like to do k-fold cross-validation too, but just using this for now 

gof_nobo_t <- mb.gof.test(rn_model_t_nobo, nsim=500, c.hat.est=TRUE, model.type="royle-nichols")


print(gof_nobo)

print(gof_nobo_t)

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
write.csv(newdata_nobo, file = "data/means_abund_parameters/means_treatment_parameters_nobo.csv", row.names = FALSE)

############Predictions for effect of treatment (repeated for t_only)########
newdata_nobo_t_only <- data.frame(treatment=levels(umf_nobo@siteCovs$treatment), shrub_cover=mean(umf_nobo@siteCovs$shrub_cover),grass_cover=mean(umf_nobo@siteCovs$grass_cover), forb_cover=mean(umf_nobo@siteCovs$forb_cover),shrub_height=mean(umf_nobo@siteCovs$shrub_height), grass_height=mean(umf_nobo@siteCovs$grass_height))

predictions_nobo_t_only <- predict(rn_model_t_only_nobo, newdata_nobo_t_only, type = "state", se.fit = TRUE)

newdata_nobo_t_only$predicted_state <- predictions_nobo_t_only$Predicted
newdata_nobo_t_only$SE <- predictions_nobo_t_only$SE
newdata_nobo_t_only$lower_CI <- newdata_nobo_t_only$predicted_state - 1.96 * newdata_nobo_t_only$SE
newdata_nobo_t_only$upper_CI <- newdata_nobo_t_only$predicted_state + 1.96 * newdata_nobo_t_only$SE

print(newdata_nobo_t_only)

write.csv(newdata_nobo_t_only, file = "data/means_abund_parameters/means_treatment_parameters_nobo_t.csv", row.names = FALSE)


##########Predictions for effect of grass cover##########

# Generate new data for grass cover while keeping other covariates constant 
grass_cover_range <- seq(min(umf_nobo@siteCovs$grass_cover), max(umf_nobo@siteCovs$grass_cover), length.out = 100)

# Check the levels of the treatment factor in the original data
treatment_levels <- levels(umf_nobo@siteCovs$treatment)

newdata_grass_cover <- data.frame(
  treatment = factor(treatment_levels, levels = treatment_levels),
  grass_cover = grass_cover_range,
  shrub_cover = mean(umf_nobo@siteCovs$shrub_cover),
  forb_cover = mean(umf_nobo@siteCovs$forb_cover),
  grass_height = mean(umf_nobo@siteCovs$grass_height),
  shrub_height = mean(umf_nobo@siteCovs$shrub_height)
)

predictions_grass_cover <- predict(rn_model_t_nobo, newdata_grass_cover, type = "state")

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
rn_model_0_nobo <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_nobo, 
  K = 25,  
  method = "BFGS" 
)

#Shrub Cover + Grass Cover + Forb Cover + Shrub Height + Grass Height 
rn_model_1_nobo <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_nobo, 
  K = 25,  
  method = "BFGS" 
)

#Treatment
rn_model_2_nobo <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) -1, 
  data = umf_nobo, 
  K = 25,  
  method = "BFGS" 
)

#Treatment * Shrub Cover * Grass Cover
rn_model_3_nobo <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment)  * scale(shrub_cover) * scale(grass_cover) -1, 
  data = umf_nobo, 
  K = 25,  
  method = "BFGS" 
)

#Treatment + Shrub Cover + Grass Cover
rn_model_4_nobo <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) - 1 + scale(shrub_cover) + scale(grass_cover), 
  data = umf_nobo, 
  K = 25,  
  method = "BFGS" 
)

#Treatment * Grass Cover 
rn_model_5_nobo <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) * scale(grass_cover) -1, 
  data = umf_nobo, 
  K = 25,  
  method = "BFGS" 
)

#Treatment + Grass Cover
rn_model_6_nobo <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) -1  + scale(grass_cover), 
  data = umf_nobo, 
  K = 25,  
  method = "BFGS" 
)

#Treatment * Shrub Cover
rn_model_7_nobo <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) * scale(shrub_cover) -1, 
  data = umf_nobo, 
  K = 25,  
  method = "BFGS" 
)

summary(rn_model_7_nobo)

#Treatment + Shrub Cover
rn_model_8_nobo <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ factor(treatment) -1  + scale(shrub_cover), 
  data = umf_nobo, 
  K = 25,  
  method = "BFGS" 
)

#Shrub Cover * Grass Cover 
rn_model_9_nobo <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) * scale(grass_cover), 
  data = umf_nobo, 
  K = 25,  
  method = "BFGS" 
)

#Shrub Cover + Grass Cover 
rn_model_10_nobo <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) + scale(grass_cover), 
  data = umf_nobo, 
  K = 25,  
  method = "BFGS" 
)

#Null
rn_model_11_nobo <- occuRN( 
  formula = ~ I(scale(day)^2) + scale(day) + scale(wind_speed) + scale(precipitation) 
  ~ 1, 
  data = umf_nobo, 
  K = 25,  
  method = "BFGS" 
)

#AICc comparison
models_list_nobo_post_hoc <- list(rn_model_0_nobo, #Treatment + Shrub Cover + Grass Cover + Forb Cover + Shrub Height + Grass Height
                                  rn_model_1_nobo, #Shrub Cover + Grass Cover + Forb Cover + Shrub Height + Grass Height
                                  rn_model_2_nobo, #Treatment
                                  rn_model_3_nobo, #Treatment * Shrub Cover * Grass Cover
                                  rn_model_4_nobo, #Treatment + Shrub Cover + Grass Cover
                                  rn_model_5_nobo, #Treatment * Grass Cover
                                  rn_model_6_nobo, #Treatment + Grass Cover
                                  rn_model_7_nobo, #Treatment * Shrub Cover
                                  rn_model_8_nobo, #Treatment + Shrub Cover
                                  rn_model_9_nobo, #Shrub Cover * Grass Cover
                                  rn_model_10_nobo, #Shrub Cover + Grass Cover
                                  rn_model_11_nobo) #Null

#model names use descriptive 
model_names_nobo_post_hoc <- c("Treatment + Shrub Cover + Grass Cover + Forb Cover + Shrub Height + Grass Height", 
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


aicc_table_nobo_post_hoc <- aictab(cand.set = models_list_nobo_post_hoc, modnames = model_names_nobo_post_hoc)

print(aicc_table_nobo_post_hoc)

# Load the master AICc table
aicc_table <- read.csv("data/aicc_table_posthoc.csv")

# Add species name to the table
aicc_table_nobo_post_hoc$species_code <- "Northern Bobwhite"

# Add the AICc table to the master table
aicc_table_master <- rbind(aicc_table, aicc_table_nobo_post_hoc)

# Save the AICc table to csv
write.csv(aicc_table_master, "data/aicc_table_posthoc.csv", row.names = FALSE)

#vif for top model
vif_table_state_nobo <- vif(rn_model_6_nobo, type = "state")
vif_table_det_nobo <- vif(rn_model_6_nobo, type = "det")

# Convert VIF tables to data frames
vif_table_state_nobo_df <- as.data.frame(vif_table_state_nobo)
vif_table_state_nobo_df$covariate <- rownames(vif_table_state_nobo_df)
vif_table_state_nobo_df$type <- "abund"
colnames(vif_table_state_nobo_df)[1] <- "vif_value"
vif_table_state_nobo_df$species <- "nobo"

vif_table_det_nobo_df <- as.data.frame(vif_table_det_nobo)
vif_table_det_nobo_df$covariate <- rownames(vif_table_det_nobo_df)
vif_table_det_nobo_df$type <- "det"
colnames(vif_table_det_nobo_df)[1] <- "vif_value"
vif_table_det_nobo_df$species <- "nobo"

# Combine the two VIF dfs
vif_table_nobo <- rbind(vif_table_state_nobo_df, vif_table_det_nobo_df)

# Reset row names to standard numbering
rownames(vif_table_nobo) <- NULL