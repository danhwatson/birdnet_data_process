# Clear Workspace
rm(list=ls())

#Load packages
library(ggplot2)
library(tidyverse)
library(viridis)


# Read in 
line_intercept_summary <- read.csv("data/line_intercept_summary.csv")


#Rename the columns
line_intercept_summary <- line_intercept_summary %>%
  rename(
    "Shrub Cover" = shrub_cover,
    "Forb Cover" = forb_cover,
    "Grass Cover" = grass_cover,
    "Conifer Cover" = conifer_cover,
    "Overall Vegetation Cover" = veg_cover_overall,
    "Overall Vegetation Height" = veg_height_overall,
    "Overall Vegetation Diversity" = veg_diversity_overall,
    "Shrub Height" = shrub_height,
    "Forb Height" = forb_height,
    "Grass Height" = grass_height,
    "Conifer Height" = conifer_height
  )



line_intercept_summary$treatment <- gsub("mine", "Mine Reclamation", line_intercept_summary$treatment)
line_intercept_summary$treatment <- gsub("timber", "Timber Production", line_intercept_summary$treatment)
line_intercept_summary$treatment <- gsub("rx_fire_young ", "Young Prescribed Fire", line_intercept_summary$treatment)
line_intercept_summary$treatment <- gsub("rx_fire_sec_growth", "Second-Growth Prescribed Fire", line_intercept_summary$treatment)

# Step 2: Check unique values to ensure all names are correct
print(unique(line_intercept_summary$treatment))

# Step 3: Convert 'treatment' to a factor with the desired order
line_intercept_summary$treatment <- factor(line_intercept_summary$treatment, 
                                           levels = c("Mine Reclamation", "Timber Production", "Young Prescribed Fire", "Second-Growth Prescribed Fire"))

# Step 4: Plotting
p1 <- line_intercept_summary %>%
  pivot_longer(cols = c(`Shrub Cover`, `Forb Cover`, `Grass Cover`, `Conifer Cover`), 
               names_to = "cover_type", values_to = "percent_cover") %>%
  ggplot(aes(x = treatment, y = percent_cover, fill = cover_type)) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "viridis", begin = 0.40, end = 1) +  # Adjusting the color range
  labs(
    title = "",
    x = "Site Treatment",
    y = "Percent Cover",
    fill = "Cover"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12)
  )

# Display the plot
p1



# Repeat for height 
p2 <- line_intercept_summary %>%
  pivot_longer(cols = c(`Shrub Height`, `Grass Height`, `Conifer Height`), 
               names_to = "cover_type", values_to = "mean_height") %>%
  ggplot(aes(x = treatment, y = mean_height, fill = cover_type)) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "viridis", begin = 0.40, end = 1) +  # Adjusting the color range
  labs(
    title = "",
    x = "Site Treatment",
    y = "Height (dm)",
    fill = "Height",
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12), 
  ) 

#panel together 
p1
p2

