# Clear workspace
rm(list=ls())

# Load required libraries
library(ggplot2)
library(tidyverse)
library(viridis)

# Read in acoustic_dat_24_07.csv
aru_data_master <- read.csv('~/Downloads/acoustic_dat_24_07.csv')

# Filter out dates before 2024-03-10 and after 2024-06-28
aru_data <- aru_data_master %>% 
  filter(date >= '2024-05-01' & date <= '2024-06-28')

#Remove sites labelled as 'NA' - remnants from having to switch ARUs in the field
aru_data <- aru_data %>% filter(site != 'NA')

# Format date as date
aru_data$date <- as.Date(aru_data$date)

# Remove observations with a confidence score or logit score 
aru_data <- aru_data %>% filter(confidence >= 0.98)

# Create a new column for the site treatment grouped by the first letter of the site name
aru_data <- aru_data %>%
  mutate(treatment = case_when(
    str_detect(site, "^M") ~ "mine",
    str_detect(site, "^T") ~ "timber",
    str_detect(site, "^R") ~ "rx_fire_young ",
    str_detect(site, "^O") ~ "rx_fire_sec_growth",
    TRUE ~ "Unknown"
  ))

# Remove observations with "nocall" in common_name
aru_data <- aru_data %>% filter(!str_detect(common_name, "nocall"))

#Remove "Florida Scrub-Jay" from the dataset
aru_data <- aru_data %>% 
  filter(unique_species != "Florida Scrub-Jay")

# Calculate the unique number of species for each treatment in aru_data
unique_species_by_treatment <- aru_data %>%
  group_by(treatment) %>%
  summarise(unique_species = n_distinct(common_name))

# Print the result
print(unique_species_by_treatment)

# Step 1: Calculate the total number of detections for each species overall and assign ranks
species_ranked <- aru_data %>%
  group_by(common_name) %>%
  summarise(total_detections = n(), .groups = 'drop') %>%
  arrange(desc(total_detections)) %>%
  mutate(rank = dense_rank(desc(total_detections)))  # Rank species overall by most detections

# Step 2: Join the species ranking with the treatment data and set desired column order
species_by_treatment <- aru_data %>%
  distinct(treatment, common_name) %>%
  left_join(species_ranked, by = "common_name") %>%
  mutate(treatment = case_when(
    treatment == "mine" ~ "Mine Reclamation",
    treatment == "timber" ~ "Timber Production",
    treatment == "rx_fire_young" ~ "Rx Fire (Young)",
    treatment == "rx_fire_sec_growth" ~ "Rx Fire (Second Growth)",
    TRUE ~ treatment)) %>%
  mutate(treatment = factor(treatment, levels = c("Mine Reclamation", "Timber Production", "Rx Fire (Young)", "Rx Fire (Second Growth)"))) %>%
  arrange(rank, treatment)

# Step 3: Reverse the y-axis order for species and arrange based on rank
species_by_treatment <- species_by_treatment %>%
  group_by(treatment) %>%
  mutate(y_pos = row_number()) %>%
  ungroup() %>%
  mutate(y_pos = max(y_pos) - y_pos + 1)  # Reverse the y_pos for displaying species from top to bottom

# Step 4: Plot the species list with treatments
ggplot(species_by_treatment, aes(x = treatment, y = y_pos)) +  # Use reversed y_pos for correct vertical order
  geom_text(aes(label = common_name), size = 4, hjust = 0.5, nudge_y = -0.1) +  # Move species list closer to x-axis labels
  geom_vline(xintercept = seq(1.5, length(unique(species_by_treatment$treatment)) - 0.5, by = 1), 
             linetype = "solid", color = "black") +  # Solid black lines between treatments
  theme_minimal() +
  labs(x = "", y = "") +  # Simplified labels, removed y-axis title
  theme(axis.text.y = element_blank(),  # Hide y-axis text (no need for it)
        axis.ticks.y = element_blank(), # Hide y-axis ticks
        panel.grid.major.y = element_blank(),  # Hide y-axis grid lines
        panel.grid.major.x = element_blank(),  # Remove vertical grid lines
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold", size = 15, vjust = -.5, color = "black"),  # Bold and lower treatment labels closer to species lists
        axis.title.x = element_text(face = "bold", size = 14)) +
  scale_x_discrete(position = "top") +  # Place treatment headers at the top
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))  # Evenly space species vertically, adjust spacing to x-axis

