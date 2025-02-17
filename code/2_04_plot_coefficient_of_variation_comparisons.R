#Calculates and plots mean coefficient of variation (CV) for relative abundance estimates across treatments for Blue Grosbeak and Bachman's Sparrow 
#Last edited 2025-02-12

# Clear workspace
rm(list = ls())

# Load required libraries
library(ggplot2)
library(dplyr)
library(patchwork)

#read in in predicted relative abundance estimates across subsamples for blgr
newdata_blgr <- read.csv("data/means_abund_parameters/means_treatment_parameters_blgr_1.csv")
newdata_2_blgr <- read.csv("data/means_abund_parameters/means_treatment_parameters_blgr_2.csv")
newdata_3_blgr <- read.csv("data/means_abund_parameters/means_treatment_parameters_blgr_3.csv")
newdata_4_blgr <- read.csv("data/means_abund_parameters/means_treatment_parameters_blgr_4.csv")
newdata_7_blgr <- read.csv("data/means_abund_parameters/means_treatment_parameters_blgr_7.csv")
newdata_10_blgr <- read.csv("data/means_abund_parameters/means_treatment_parameters_blgr_10.csv")
#now bacs
newdata_bacs <- read.csv("data/means_abund_parameters/means_treatment_parameters_bacs_1.csv")
newdata_2_bacs <- read.csv("data/means_abund_parameters/means_treatment_parameters_bacs_2.csv")
newdata_3_bacs <- read.csv("data/means_abund_parameters/means_treatment_parameters_bacs_3.csv")
newdata_4_bacs <- read.csv("data/means_abund_parameters/means_treatment_parameters_bacs_4.csv")
newdata_7_bacs <- read.csv("data/means_abund_parameters/means_treatment_parameters_bacs_7.csv")
newdata_10_bacs <- read.csv("data/means_abund_parameters/means_treatment_parameters_bacs_10.csv")

# Function to calculate mean CV across treatments
calculate_mean_cv <- function(data_frame) {
  mean(data_frame$SE / data_frame$predicted_state, na.rm = TRUE)  # Compute mean CV
}

# Calculate mean CV for each survey interval for BLGR
mean_cv_every_day_blgr <- calculate_mean_cv(newdata_blgr)
mean_cv_2_days_blgr <- calculate_mean_cv(newdata_2_blgr)
mean_cv_3_days_blgr <- calculate_mean_cv(newdata_3_blgr)
mean_cv_4_days_blgr <- calculate_mean_cv(newdata_4_blgr)
mean_cv_7_days_blgr <- calculate_mean_cv(newdata_7_blgr)
mean_cv_10_days_blgr <- calculate_mean_cv(newdata_10_blgr)

# Combine results into a data frame for BLGR
cv_summary_blgr <- data.frame(
  Interval = factor(c("Everyday", "Every 2nd Day", "Every 3rd Day", 
                      "Every 4th Day", "Every 7th Day", "Every 10th Day"), 
                    levels = c("Everyday", "Every 2nd Day", "Every 3rd Day", 
                               "Every 4th Day", "Every 7th Day", "Every 10th Day")),
  Mean_CV = c(mean_cv_every_day_blgr, mean_cv_2_days_blgr, mean_cv_3_days_blgr, 
              mean_cv_4_days_blgr, mean_cv_7_days_blgr, mean_cv_10_days_blgr)
)

# Repeat for BACS
mean_cv_every_day_bacs <- calculate_mean_cv(newdata_bacs)
mean_cv_2_days_bacs <- calculate_mean_cv(newdata_2_bacs)
mean_cv_3_days_bacs <- calculate_mean_cv(newdata_3_bacs)
mean_cv_4_days_bacs <- calculate_mean_cv(newdata_4_bacs)
mean_cv_7_days_bacs <- calculate_mean_cv(newdata_7_bacs)
mean_cv_10_days_bacs <- calculate_mean_cv(newdata_10_bacs)

cv_summary_bacs <- data.frame(
  Interval = factor(c("Everyday", "Every 2nd Day", "Every 3rd Day", 
                      "Every 4th Day", "Every 7th Day", "Every 10th Day"), 
                    levels = c("Everyday", "Every 2nd Day", "Every 3rd Day", 
                               "Every 4th Day", "Every 7th Day", "Every 10th Day")),
  Mean_CV = c(mean_cv_every_day_bacs, mean_cv_2_days_bacs, mean_cv_3_days_bacs, 
              mean_cv_4_days_bacs, mean_cv_7_days_bacs, mean_cv_10_days_bacs)
)

# Combine BLGR and BACS results into one data frame
cv_summary_combined <- data.frame(
  Interval = factor(
    rep(c("Everyday", "Every 2nd Day", "Every 3rd Day", "Every 4th Day", 
          "Every 7th Day", "Every 10th Day"), 2),
    levels = c("Everyday", "Every 2nd Day", "Every 3rd Day", "Every 4th Day", 
               "Every 7th Day", "Every 10th Day")
  ),
  Mean_CV = c(mean_cv_every_day_blgr, mean_cv_2_days_blgr, mean_cv_3_days_blgr,
              mean_cv_4_days_blgr, mean_cv_7_days_blgr, mean_cv_10_days_blgr,
              mean_cv_every_day_bacs, mean_cv_2_days_bacs, mean_cv_3_days_bacs,
              mean_cv_4_days_bacs, mean_cv_7_days_bacs, mean_cv_10_days_bacs),
  Species = factor(rep(c("Blue Grosbeak", "Bachman's Sparrow"), each = 6),
                   levels = c("Blue Grosbeak", "Bachman's Sparrow")) # Reorder levels
)

# Define custom colors
species_colors <- c("Blue Grosbeak" = "#4F8FC4", "Bachman's Sparrow" = "#C4904E")

# Plot combined CV data
ggplot(cv_summary_combined, aes(x = Interval, y = Mean_CV, group = Species, color = Species)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = species_colors) +
  labs(x = "", 
       y = "Mean Coefficient of Variation (CV) for Relative Abundance Estimates") +
  theme_classic() +
  theme(
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = c(0.70, 0.1),
    legend.justification = c(0, 0),
    legend.box.background = element_rect(color = "black", size = 0.5),
    legend.box.margin = margin(5, 5, 5, 5),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 14),
    axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 2),
    axis.text.y = element_text(size = 14)
  )

# Save the combined CV summary
write.csv(cv_summary_combined, file = "data/cv_summary_combined.csv", row.names = FALSE)

