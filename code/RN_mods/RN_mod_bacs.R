# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(AICcmodavg)  
library(pROC)


##########Loading data and formatting for RN mods###########

# Load bacs abundance data
load("data/abundance_data/abundance_24_bacs.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_bacs$y <- ifelse(abund_data_bacs$y > 0 | is.na(abund_data_bacs$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as these RN models can't handle NAs we have in our dataset 
abund_data_bacs$y <- abund_data_bacs$y[, -c(1:55)]
abund_data_bacs$det.covs$day <- abund_data_bacs$det.covs$day[, -c(1:55)]
abund_data_bacs$det.covs$temp <- abund_data_bacs$det.covs$temp[, -c(1:55)]
abund_data_bacs$det.covs$wind_speed <- abund_data_bacs$det.covs$wind_speed[, -c(1:55)]
abund_data_bacs$det.covs$precipitation <- abund_data_bacs$det.covs$precipitation[, -c(1:55)]


##########Fitting Royle-Nichols models###########

# Create the unmarkedFrameOccu object needed to be used for all models made with unmarked
umf_bacs <- unmarkedFrameOccu(
  y = as.matrix(abund_data_bacs$y),
  siteCovs = abund_data_bacs$abund.covs,
  obsCovs = list(
    day = as.matrix(abund_data_bacs$det.covs$day),
    temp = as.matrix(abund_data_bacs$det.covs$temp),
    wind_speed = as.matrix(abund_data_bacs$det.covs$wind_speed),
    precipitation = as.matrix(abund_data_bacs$det.covs$precipitation)
  )
)

# Fit the first Royle-Nichols model (with site treatment as an abundance covariate, t for treatment)
rn_model_t_bacs <- occuRN( # This is the unmarked function to fit Royle-Nichols (RN) models
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) # Detection model covariates
  ~ factor(treatment) -1  + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), # Abundance model covariates 
  data = umf_bacs, # Specifying the data from the unmarkedFrameOccu object we created
  K = 25,  # Ensure this is large enough to cover possible abundances
  method = "BFGS"  # Optimization method
)

summary(rn_model_t_bacs) # Look at the summary output of the model
# Models have '-1' term in the abundance model formula to remove the intercept, therefore providing mean parameter estimates for each level of the factor (treatment)

saveRDS(rn_model_t_bacs, file = "data/rn_model_files/rn_model_t_bacs.rds") 
# Optional step to save model object to load in another R session


# Fit the second Royle-Nichols model (without site treatment as an abundance covariate)
rn_model_bacs <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
  ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height), 
  data = umf_bacs,
  K = 25,  
  method = "BFGS" 
)
summary(rn_model_bacs)

# Fit the null model
rn_model_null_bacs <- occuRN(
  formula = ~ I(scale(day)^2) + scale(day) + scale(temp) + scale(wind_speed) + scale(precipitation) 
  ~ 1, 
  data = umf_bacs,
  K = 25,  
  method = "BFGS"  
)
summary(rn_model_null_bacs)

########AICc comparison and AICc table formatting###########

# Load the master AICc table
aicc_table <- read.csv("data/aicc_table.csv")

models_list_bacs <- list(rn_model_t_bacs = rn_model_t_bacs, rn_model_bacs = rn_model_bacs, rn_model_null_bacs = rn_model_null_bacs)
model_names_bacs <- c("Treatment Model", "Non-treatment Model", "Null Model")
aicc_table_bacs <- aictab(cand.set = models_list_bacs, modnames = model_names_bacs)
print(aicc_table_bacs)

# Add species name to the table
aicc_table_bacs$species_code <- "bacs"

# Add the AICc table to the master table
aicc_table_master <- rbind(aicc_table, aicc_table_bacs)

# Save the AICc table to csv
write.csv(aicc_table_master, "data/aicc_table.csv", row.names = FALSE)

##########Goodness of fit##########

gof <- mb.gof.test(rn_model_t_bacs, nsim=100, c.hat.est=TRUE, model.type="royle-nichols")
# Looking more into GoF tests, would like to do k-fold cross-validation too, but just using this for now 
print(gof)


##########Predictions for effect of treatment##########

# Create a new data frame with the mean values of the covariates
newdata_bacs <- data.frame(treatment=levels(umf_bacs@siteCovs$treatment), shrub_cover=mean(umf_bacs@siteCovs$shrub_cover),grass_cover=mean(umf_bacs@siteCovs$grass_cover), forb_cover=mean(umf_bacs@siteCovs$forb_cover),shrub_height=mean(umf_bacs@siteCovs$shrub_height), grass_height=mean(umf_bacs@siteCovs$grass_height))

# Generate predictions for the new data with the predict function
predictions_bacs <- predict(rn_model_t_bacs, newdata_bacs, type = "state", se.fit = TRUE)

# Add predictions and SE to the newdata_bacs dataframe
newdata_bacs$predicted_state <- predictions_bacs$Predicted
newdata_bacs$SE <- predictions_bacs$SE

# Calculate 95% confidence intervals from SE
newdata_bacs$lower_CI <- newdata_bacs$predicted_state - 1.96 * newdata_bacs$SE
newdata_bacs$upper_CI <- newdata_bacs$predicted_state + 1.96 * newdata_bacs$SE

# View the results
print(newdata_bacs)

# Optionally, save results as csv for plotting in another R session
write.csv(newdata_bacs, file = "data/means_treatment_parameters_bacs.csv", row.names = FALSE)

##########Predictions for effect of grass cover##########

# Generate new data for grass cover while keeping other covariates constant 
grass_cover_range <- seq(min(umf_bacs@siteCovs$grass_cover), max(umf_bacs@siteCovs$grass_cover), length.out = 100)

# Check the levels of the treatment factor in the original data
treatment_levels <- levels(umf_bacs@siteCovs$treatment)

newdata_grass_cover <- data.frame(
  treatment = factor(treatment_levels, levels = treatment_levels),
  grass_cover = grass_cover_range,
  shrub_cover = mean(umf_bacs@siteCovs$shrub_cover),
  forb_cover = mean(umf_bacs@siteCovs$forb_cover),
  grass_height = mean(umf_bacs@siteCovs$grass_height),
  shrub_height = mean(umf_bacs@siteCovs$shrub_height)
)

predictions_grass_cover <- predict(rn_model_t_bacs, newdata_grass_cover, type = "state")

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
dir.create("figures/predicted_abundance/mine/bacs")
ggsave("figures/predicted_abundance/mine/bacs/mine_grass_cover_bacs.png", width = 4, height = 4)

ggplot(plot_data_grass_cover_timber, aes(x = grass_cover, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/timber
dir.create("figures/predicted_abundance/timber/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/timber/bacs/timber_grass_cover_bacs.png", width = 4, height = 4)

ggplot(plot_data_grass_cover_rx_fire_young, aes(x = grass_cover, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/rx_fire_young
dir.create("figures/predicted_abundance/rx_fire_young/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/rx_fire_young/bacs/rx_fire_young_grass_cover_bacs.png", width = 4, height = 5)

ggplot(plot_data_grass_cover_rx_fire_sec_growth, aes(x = grass_cover, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/rx_fire_sec_growth
dir.create("figures/predicted_abundance/rx_fire_sec_growth/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/rx_fire_sec_growth/bacs/rx_fire_sec_growth_grass_cover_bacs.png", width = 4, height = 4)


#repeat for shrub cover
shrub_cover_range <- seq(min(umf_bacs@siteCovs$shrub_cover), max(umf_bacs@siteCovs$shrub_cover), length.out = 100)

newdata_shrub_cover <- data.frame(
  treatment = factor(treatment_levels, levels = treatment_levels),
  shrub_cover = shrub_cover_range,
  grass_cover = mean(umf_bacs@siteCovs$grass_cover),
  forb_cover = mean(umf_bacs@siteCovs$forb_cover),
  grass_height = mean(umf_bacs@siteCovs$grass_height),
  shrub_height = mean(umf_bacs@siteCovs$shrub_height)
)

predictions_shrub_cover <- predict(rn_model_t_bacs, newdata_shrub_cover, type = "state")

plot_data_shrub_cover <- data.frame(
  shrub_cover = newdata_shrub_cover$shrub_cover,
  treatment = newdata_shrub_cover$treatment,
  predicted_state = predictions_shrub_cover$Predicted,
  lower_CI = predictions_shrub_cover$lower,
  upper_CI = predictions_shrub_cover$upper
)

ggplot(plot_data_shrub_cover, aes(x = shrub_cover, y = predicted_state, color = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill = treatment), alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

plot_data_shrub_cover_mine <- plot_data_shrub_cover %>% filter(treatment == "mine")
plot_data_shrub_cover_timber <- plot_data_shrub_cover %>% filter(treatment == "timber")
plot_data_shrub_cover_rx_fire_young <- plot_data_shrub_cover %>% filter(treatment == "rx_fire_young")
plot_data_shrub_cover_rx_fire_sec_growth <- plot_data_shrub_cover %>% filter(treatment == "rx_fire_sec_growth")

ggplot(plot_data_shrub_cover_mine, aes(x = shrub_cover, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/mine
dir.create("figures/predicted_abundance/mine/bacs")
ggsave("figures/predicted_abundance/mine/bacs/mine_shrub_cover_bacs.png", width = 4, height = 4)

ggplot(plot_data_shrub_cover_timber, aes(x = shrub_cover, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/timber
dir.create("figures/predicted_abundance/timber/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/timber/bacs/timber_shrub_cover_bacs.png", width = 4, height = 4)

ggplot(plot_data_shrub_cover_rx_fire_young, aes(x = shrub_cover, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/rx_fire_young
dir.create("figures/predicted_abundance/rx_fire_young/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/rx_fire_young/bacs/rx_fire_young_shrub_cover_bacs.png", width = 4, height = 4)

ggplot(plot_data_shrub_cover_rx_fire_sec_growth, aes(x = shrub_cover, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/rx_fire_sec_growth
dir.create("figures/predicted_abundance/rx_fire_sec_growth/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/rx_fire_sec_growth/bacs/rx_fire_sec_growth_shrub_cover_bacs.png", width = 4, height = 4)

#repeat for forb cover

forb_cover_range <- seq(min(umf_bacs@siteCovs$forb_cover), max(umf_bacs@siteCovs$forb_cover), length.out = 100)

newdata_forb_cover <- data.frame(
  treatment = factor(treatment_levels, levels = treatment_levels),
  forb_cover = forb_cover_range,
  grass_cover = mean(umf_bacs@siteCovs$grass_cover),
  shrub_cover = mean(umf_bacs@siteCovs$shrub_cover),
  grass_height = mean(umf_bacs@siteCovs$grass_height),
  shrub_height = mean(umf_bacs@siteCovs$shrub_height)
)

predictions_forb_cover <- predict(rn_model_t_bacs, newdata_forb_cover, type = "state")

plot_data_forb_cover <- data.frame(
  forb_cover = newdata_forb_cover$forb_cover,
  treatment = newdata_forb_cover$treatment,
  predicted_state = predictions_forb_cover$Predicted,
  lower_CI = predictions_forb_cover$lower,
  upper_CI = predictions_forb_cover$upper
)

ggplot(plot_data_forb_cover, aes(x = forb_cover, y = predicted_state, color = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill = treatment), alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

plot_data_forb_cover_mine <- plot_data_forb_cover %>% filter(treatment == "mine")
plot_data_forb_cover_timber <- plot_data_forb_cover %>% filter(treatment == "timber")
plot_data_forb_cover_rx_fire_young <- plot_data_forb_cover %>% filter(treatment == "rx_fire_young")
plot_data_forb_cover_rx_fire_sec_growth <- plot_data_forb_cover %>% filter(treatment == "rx_fire_sec_growth")

ggplot(plot_data_forb_cover_mine, aes(x = forb_cover, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/mine
dir.create("figures/predicted_abundance/mine/bacs")
ggsave("figures/predicted_abundance/mine/bacs/mine_forb_cover_bacs.png", width = 4, height = 4)

ggplot(plot_data_forb_cover_timber, aes(x = forb_cover, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/timber
dir.create("figures/predicted_abundance/timber/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/timber/bacs/timber_forb_cover_bacs.png", width = 4, height = 4)

ggplot(plot_data_forb_cover_rx_fire_young, aes(x = forb_cover, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/rx_fire_young
dir.create("figures/predicted_abundance/rx_fire_young/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/rx_fire_young/bacs/rx_fire_young_forb_cover_bacs.png", width = 4, height = 4)

ggplot(plot_data_forb_cover_rx_fire_sec_growth, aes(x = forb_cover, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/rx_fire_sec_growth
dir.create("figures/predicted_abundance/rx_fire_sec_growth/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/rx_fire_sec_growth/bacs/rx_fire_sec_growth_forb_cover_bacs.png", width = 4, height = 4)

#repeat for shrub height

shrub_height_range <- seq(min(umf_bacs@siteCovs$shrub_height), max(umf_bacs@siteCovs$shrub_height), length.out = 100)

newdata_shrub_height <- data.frame(
  treatment = factor(treatment_levels, levels = treatment_levels),
  shrub_height = shrub_height_range,
  grass_cover = mean(umf_bacs@siteCovs$grass_cover),
  shrub_cover = mean(umf_bacs@siteCovs$shrub_cover),
  forb_cover = mean(umf_bacs@siteCovs$forb_cover),
  grass_height = mean(umf_bacs@siteCovs$grass_height)
)

predictions_shrub_height <- predict(rn_model_t_bacs, newdata_shrub_height, type = "state")

plot_data_shrub_height <- data.frame(
  shrub_height = newdata_shrub_height$shrub_height,
  treatment = newdata_shrub_height$treatment,
  predicted_state = predictions_shrub_height$Predicted,
  lower_CI = predictions_shrub_height$lower,
  upper_CI = predictions_shrub_height$upper
)

ggplot(plot_data_shrub_height, aes(x = shrub_height, y = predicted_state, color = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill = treatment), alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

plot_data_shrub_height_mine <- plot_data_shrub_height %>% filter(treatment == "mine")
plot_data_shrub_height_timber <- plot_data_shrub_height %>% filter(treatment == "timber")
plot_data_shrub_height_rx_fire_young <- plot_data_shrub_height %>% filter(treatment == "rx_fire_young")
plot_data_shrub_height_rx_fire_sec_growth <- plot_data_shrub_height %>% filter(treatment == "rx_fire_sec_growth")

ggplot(plot_data_shrub_height_mine, aes(x = shrub_height, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/mine
dir.create("figures/predicted_abundance/mine/bacs")
ggsave("figures/predicted_abundance/mine/bacs/mine_shrub_height_bacs.png", width = 4, height = 4)

ggplot(plot_data_shrub_height_timber, aes(x = shrub_height, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/timber
dir.create("figures/predicted_abundance/timber/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/timber/bacs/timber_shrub_height_bacs.png", width = 4, height = 4)

ggplot(plot_data_shrub_height_rx_fire_young, aes(x = shrub_height, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/rx_fire_young
dir.create("figures/predicted_abundance/rx_fire_young/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/rx_fire_young/bacs/rx_fire_young_shrub_height_bacs.png", width = 4, height = 4)

ggplot(plot_data_shrub_height_rx_fire_sec_growth, aes(x = shrub_height, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/rx_fire_sec_growth
dir.create("figures/predicted_abundance/rx_fire_sec_growth/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/rx_fire_sec_growth/bacs/rx_fire_sec_growth_shrub_height_bacs.png", width = 4, height = 4)


# Repeat for grass height

grass_height_range <- seq(min(umf_bacs@siteCovs$grass_height), max(umf_bacs@siteCovs$grass_height), length.out = 100)

newdata_grass_height <- data.frame(
  treatment = factor(treatment_levels, levels = treatment_levels),
  grass_height = grass_height_range,
  grass_cover = mean(umf_bacs@siteCovs$grass_cover),
  shrub_cover = mean(umf_bacs@siteCovs$shrub_cover),
  forb_cover = mean(umf_bacs@siteCovs$forb_cover),
  shrub_height = mean(umf_bacs@siteCovs$shrub_height)
)

predictions_grass_height <- predict(rn_model_t_bacs, newdata_grass_height, type = "state")

plot_data_grass_height <- data.frame(
  grass_height = newdata_grass_height$grass_height,
  treatment = newdata_grass_height$treatment,
  predicted_state = predictions_grass_height$Predicted,
  lower_CI = predictions_grass_height$lower,
  upper_CI = predictions_grass_height$upper
)

ggplot(plot_data_grass_height, aes(x = grass_height, y = predicted_state, color = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill = treatment), alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

plot_data_grass_height_mine <- plot_data_grass_height %>% filter(treatment == "mine")
plot_data_grass_height_timber <- plot_data_grass_height %>% filter(treatment == "timber")
plot_data_grass_height_rx_fire_young <- plot_data_grass_height %>% filter(treatment == "rx_fire_young")
plot_data_grass_height_rx_fire_sec_growth <- plot_data_grass_height %>% filter(treatment == "rx_fire_sec_growth")

ggplot(plot_data_grass_height_mine, aes(x = grass_height, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/mine
dir.create("figures/predicted_abundance/mine/bacs")
ggsave("figures/predicted_abundance/mine/bacs/mine_grass_height_bacs.png", width = 4, height = 4)

ggplot(plot_data_grass_height_timber, aes(x = grass_height, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/timber
dir.create("figures/predicted_abundance/timber/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/timber/bacs/timber_grass_height_bacs.png", width = 4, height = 4)

ggplot(plot_data_grass_height_rx_fire_young, aes(x = grass_height, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/rx_fire_young
dir.create("figures/predicted_abundance/rx_fire_young/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/rx_fire_young/bacs/rx_fire_young_grass_height_bacs.png", width = 4, height = 4)

ggplot(plot_data_grass_height_rx_fire_sec_growth, aes(x = grass_height, y = predicted_state)) +
  geom_line(linewidth = 1, color = "black") +  # Customize the color if needed
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "darkgray", alpha = 0.2) +
  labs(x = "", y = "", title = "") +
  theme_classic()

#save to figures/predicted_abundance/rx_fire_sec_growth
dir.create("figures/predicted_abundance/rx_fire_sec_growth/bacs", showWarnings = FALSE)
ggsave("figures/predicted_abundance/rx_fire_sec_growth/bacs/rx_fire_sec_growth_grass_height_bacs.png", width = 4, height = 4)
                                     
                  