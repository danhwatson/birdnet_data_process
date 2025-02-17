#Plotting cumulative detection probabilities with 95% confidence intervals for Blue Grosbeak and Bachman's Sparrow across sub-samples
#Last edited: 2025-02-13

#Clear workspace
rm(list = ls())

# Load required libraries
library(ggplot2)
library(tidyverse)

#load models
load("data/rn_models/rn_subsample_models_blgr.RData")
load("data/rn_models/rn_subsample_models_bacs.RData")

# Define a function to calculate cumulative detection probability
calculate_cumulative_detection <- function(model) {
  # Extract detection probabilities from the model
  detection_probs <- fitted(model)
  
  # Check if the extraction worked correctly
  if (is.null(detection_probs)) {
    stop("Unable to extract detection probabilities from the model.")
  }
  
  # Calculate cumulative detection probability: P = 1 - prod(1 - p_i)
  cumulative_detection <- 1 - apply(1 - detection_probs, 1, prod, na.rm = TRUE)
  
  return(cumulative_detection)
}

# Define a function to calculate confidence intervals and standard error for cumulative detection probabilities
calculate_confidence_intervals <- function(cumulative_detection) {
  mean_value <- mean(cumulative_detection, na.rm = TRUE)
  se <- sd(cumulative_detection, na.rm = TRUE) / sqrt(length(na.omit(cumulative_detection)))
  lower_ci <- mean_value - 1.96 * se
  upper_ci <- mean_value + 1.96 * se
  return(c(mean = mean_value, se = se, lower_ci = lower_ci, upper_ci = upper_ci))
}

# List of models for each dataset and species
models <- list(
  every_day_blgr = every_day_blgr,
  every_2_days_blgr = every_2_days_blgr,
  every_3_days_blgr = every_3_days_blgr,
  every_4_days_blgr = every_4_days_blgr,
  every_7_days_blgr = every_7_days_blgr,
  every_10_days_blgr = every_10_days_blgr,
  every_day_bacs = every_day_bacs,
  every_2_days_bacs = every_2_days_bacs,
  every_3_days_bacs = every_3_days_bacs,
  every_4_days_bacs = every_4_days_bacs,
  every_7_days_bacs = every_7_days_bacs,
  every_10_days_bacs = every_10_days_bacs
)

# Add dataset column to each entry
datasets <- c(
  "Everyday",
  "Every 2nd Day",
  "Every 3rd Day",
  "Every 4th Day",
  "Every 7th Day",
  "Every 10th Day"
)

data_mapping <- c(
  rep("Blue Grosbeak", 6),
  rep("Bachman's Sparrow", 6)
)

# Create cumulative detection probabilities dataframe
cumulative_detection_df <- do.call(rbind, lapply(1:length(models), function(i) {
  model_name <- names(models)[i]
  cumulative_probs <- calculate_cumulative_detection(models[[model_name]])
  ci <- calculate_confidence_intervals(cumulative_probs)
  data.frame(
    Model = model_name,
    Dataset = factor(datasets[(i - 1) %% 6 + 1], levels = c(
      "Everyday",
      "Every 2nd Day",
      "Every 3rd Day",
      "Every 4th Day",
      "Every 7th Day",
      "Every 10th Day"
    )),
    Species = data_mapping[i],
    Mean_Cumulative_Detection = ci["mean"],
    SE = ci["se"],
    Lower_CI = ci["lower_ci"],
    Upper_CI = ci["upper_ci"]
  )
}))

# Reorder species so Blue Grosbeak appears first
cumulative_detection_df$Species <- factor(cumulative_detection_df$Species, 
                                          levels = c("Blue Grosbeak", "Bachman's Sparrow"))

# Save cumulative detection probabilities, confidence intervals, and SE to a CSV file
write.csv(cumulative_detection_df, "data/cumulative_detection_probabilities_with_CI_and_SE.csv", row.names = FALSE)

# Print the cumulative detection probabilities, confidence intervals, and SE
print(cumulative_detection_df)


# Plot cumulative detection probabilities with confidence intervals
combined_plot <- ggplot(cumulative_detection_df, aes(x = Dataset, y = Mean_Cumulative_Detection, color = Species, group = Species)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.3, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = seq(1.5, length(unique(cumulative_detection_df$Dataset)) - 0.5, by = 1), 
             linetype = "dashed", color = "gray", alpha = 0.7) +
  scale_color_manual(
    values = c("Blue Grosbeak" = "#4F8FC4", "Bachman's Sparrow" = "#C4904E")
  ) +
  labs(
    title = "",
    x = "",
    y = "Cumulative Detection Probability",
    color = "Species"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "black", angle = 45, hjust = .6, vjust = .70, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 13),
    legend.position = c(0.11, 0.1),  # Position legend inside plot (x, y)
    legend.justification = c(0, 0),  # Anchor legend to bottom-left
    legend.box.background = element_rect(color = "black", size = 0.5),  # Add box around legend
    legend.box.margin = margin(5, 5, 5, 5)  # Add padding inside legend box
  )

print(combined_plot)
