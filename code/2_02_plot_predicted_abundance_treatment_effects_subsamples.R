# Clear workspace
rm(list=ls())

# Load required libraries
library(ggplot2)
library(tidyverse)
library(viridis)
library(patchwork)

# Load all the CSV files from the 'means_parameter' folder
mean_blgr_1 <- read_csv("data/means_abund_parameters/means_treatment_parameters_blgr_1.csv")
mean_blgr_2 <- read_csv("data/means_abund_parameters/means_treatment_parameters_blgr_2.csv")
mean_blgr_3 <- read_csv("data/means_abund_parameters/means_treatment_parameters_blgr_3.csv")
mean_blgr_4 <- read_csv("data/means_abund_parameters/means_treatment_parameters_blgr_4.csv")
mean_blgr_7 <- read_csv("data/means_abund_parameters/means_treatment_parameters_blgr_7.csv")
mean_blgr_10 <- read_csv("data/means_abund_parameters/means_treatment_parameters_blgr_10.csv")

#Add dataset column to each df 
mean_blgr_1$dataset <- "Everyday"
mean_blgr_2$dataset <- "Every 2nd Day"
mean_blgr_3$dataset <- "Every 3rd Day"
mean_blgr_4$dataset <- "Every 4th Day"
mean_blgr_7$dataset <- "Every 7th Day"
mean_blgr_10$dataset <- "Every 10th Day"


# Combine all the dataframes into one
combined_data <- bind_rows(mean_blgr_1, mean_blgr_2, mean_blgr_3, mean_blgr_4, mean_blgr_7, mean_blgr_10)

#rename columns
#Rename the columns
combined_data <- combined_data %>%
  mutate(treatment = recode(treatment,
                            "mine" = "Mine",
                            "mature_savanna" = "Mature Savanna",
                            "young_savanna" = "Young Savanna",
                            "timber" = "Timber"
  )) 

# Ensure the treatments are in particular order
combined_data$treatment <- factor(combined_data$treatment, 
                                      levels = (c(
                            "Mine", "Timber", "Young Savanna", "Mature Savanna")))



# Ensure the datasets are in particular order
combined_data$dataset <- factor(combined_data$dataset, 
                                levels = c(
                                  "Everyday",
                                  "Every 2nd Day",
                                  "Every 3rd Day",
                                  "Every 4th Day",
                                  "Every 7th Day",
                                  "Every 10th Day"
                                ))


#save the combined data to a new CSV file
write_csv(combined_data, "data/means_abund_parameters/means_treatment_parameters_blgr_combined.csv")

# Create the plot
blgr <- ggplot(combined_data, aes(x = treatment, y = predicted_state, color = dataset)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), 
                width = 0.2, 
                position = position_dodge(width = 0.5)) + # Error bars aligned with points
scale_color_viridis(option = "viridis", begin = 0.05, end = 0.95, discrete = TRUE) +
geom_vline(xintercept = seq(1.5, length(unique(combined_data$treatment)) - 0.5, by = 1), 
             linetype = "dashed", color = "gray", alpha = 0.7) + # Add dashed vertical lines
  labs(
    title = "Blue Grosbeak",
    x = "Site Treatment",
    y = "Predicted Abundance",
    color = "Dataset"
  ) +
  theme_classic() + # Apply the classic theme
  theme(
    axis.text.x = element_text(angle = 0, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14)
  ) 
  #+
  #guides(color = guide_legend(reverse = TRUE))



# Load all the CSV files from the 'means_parameter' folder
mean_bacs_1 <- read_csv("data/means_abund_parameters/means_treatment_parameters_bacs_1.csv")
mean_bacs_2 <- read_csv("data/means_abund_parameters/means_treatment_parameters_bacs_2.csv")
mean_bacs_3 <- read_csv("data/means_abund_parameters/means_treatment_parameters_bacs_3.csv")
mean_bacs_4 <- read_csv("data/means_abund_parameters/means_treatment_parameters_bacs_4.csv")
mean_bacs_7 <- read_csv("data/means_abund_parameters/means_treatment_parameters_bacs_7.csv")
mean_bacs_10 <- read_csv("data/means_abund_parameters/means_treatment_parameters_bacs_10.csv")

#Add dataset column to each df 
mean_bacs_1$dataset <- "Everyday"
mean_bacs_2$dataset <- "Every 2nd Day"
mean_bacs_3$dataset <- "Every 3rd Day"
mean_bacs_4$dataset <- "Every 4th Day"
mean_bacs_7$dataset <- "Every 7th Day"
mean_bacs_10$dataset <- "Every 10th Day"

# Combine all the dataframes into one
combined_data <- bind_rows(mean_bacs_1, mean_bacs_2, mean_bacs_3, mean_bacs_4, mean_bacs_7, mean_bacs_10)

#rename columns
#Rename the columns
combined_data <- combined_data %>%
  mutate(treatment = recode(treatment,
                            "mine" = "Mine",
                            "mature_savanna" = "Mature Savanna",
                            "young_savanna" = "Young Savanna",
                            "timber" = "Timber"
  )) 

# Ensure the treatments are in particular order
combined_data$treatment <- factor(combined_data$treatment, 
                                  levels = (c(
                                    "Mine", "Timber", "Young Savanna", "Mature Savanna")))


# Ensure the datasets are in particular order
combined_data$dataset <- factor(combined_data$dataset, 
                                levels = c(
                                  "Everyday",
                                  "Every 2nd Day",
                                  "Every 3rd Day",
                                  "Every 4th Day",
                                  "Every 7th Day",
                                  "Every 10th Day"
                                ))



#save the combined data to a new CSV file
write_csv(combined_data, "data/means_abund_parameters/means_treatment_parameters_bacs_combined.csv")

# Create the plot
bacs <- ggplot(combined_data, aes(x = treatment, y = predicted_state, color = dataset)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), 
                width = 0.2, 
                position = position_dodge(width = 0.5)) + # Error bars aligned with points
  scale_color_viridis(option = "viridis", begin = 0.05, end = 0.95, discrete = TRUE) +
  geom_vline(xintercept = seq(1.5, length(unique(combined_data$treatment)) - 0.5, by = 1), 
             linetype = "dashed", color = "gray", alpha = 0.7) + # Add dashed vertical lines
  labs(
    title = "Bachman's Sparrow",
    x = "Site Treatment",
    y = "Predicted Abundance",
    color = "Dataset"
  ) +
  theme_classic() + # Apply the classic theme
  theme(
    axis.text.x = element_text(angle = 0, size = 14, face = "bold"),
    axis.text.y = element_text(size = 12,),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14)
  ) 
#  +
  #guides(color = guide_legend(reverse = TRUE))

print(bacs)

library(cowplot)
library(grid)

# Extract the legend from one of the plots
legend <- cowplot::get_legend(
  blgr + 
    theme(legend.position = "right", legend.box.background = element_rect(color = "black", size = 1))
)

# Adjust the combined plot to exclude individual legends
combined_plot <- (
  blgr + 
    labs(y = NULL, x = NULL) + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(),  # Remove x-axis text
          axis.ticks.x = element_blank())  # Remove x-axis ticks
) / (
  bacs + 
    labs(y = NULL, x = NULL) + 
    theme(legend.position = "none") 
)

# Add shared y-axis label with properly centered N̂
y_label <- grid::textGrob(
  label = "Predicted Relative Abundance", 
  gp = gpar(fontsize = 14, fontface = "bold"), 
  rot = 90
)



# Combine plots and place the legend in the top right corner
final_plot <- wrap_elements(full = y_label) + combined_plot + 
  plot_layout(widths = c(0.01, 1)) + # Adjust widths to make space for the y-axis label
  inset_element(
    legend, 
    left = 0.8, bottom = 0.8, right = 1, top = 1,  # Place the legend in the top right corner
    align_to = "panel"
  )

# Display the final plot
print(final_plot)

