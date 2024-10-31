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

# Step 1: Calculate the total number of unique detections for each species (remove duplicates)
species_ranked <- aru_data %>%
  distinct(treatment, common_name) %>%  # Ensure uniqueness by species and treatment
  group_by(common_name) %>%
  summarise(total_detections = n()) %>%
  arrange(desc(total_detections)) %>%
  mutate(rank = row_number())

# Step 2: Rename the treatments manually
species_by_treatment <- aru_data %>%
  distinct(treatment, common_name) %>%
  left_join(species_ranked, by = "common_name") %>%
  mutate(treatment = recode(treatment,
                            "mine" = "Mine Reclamation",
                            "timber" = "Timber Production",
                            "rx_fire_young" = "Rx Fire (Young)",
                            "rx_fire_sec_growth" = "Rx Fire (Second Growth)")) %>%
  arrange(treatment, rank)

# Step 3: Create a rank for the species within each treatment to position them on the y-axis
species_by_treatment <- species_by_treatment %>%
  group_by(treatment) %>%
  mutate(y_pos = row_number())

# Step 4: Adjust plotting to ensure even spacing and centering
ggplot(species_by_treatment, aes(x = treatment, y = -rank)) +  # Use negative rank to show most detected at the top
  geom_text(aes(label = common_name), size = 4, hjust = 0.5) +  # Set hjust to 0.5 for centered text
  geom_vline(xintercept = seq(1.5, length(unique(species_by_treatment$treatment)) - 0.5, by = 1), linetype = "dashed", color = "gray") +  # Lines between treatments
  theme_minimal() +
  labs(x = "Treatment", y = "Species Rank",
       title = "Ranked Unique Species by Treatment") +  # Changed "TREATMENT" to "Treatment"
  theme(axis.text.y = element_blank(),  # Hide y-axis text
        axis.ticks.y = element_blank(), # Hide y-axis ticks
        panel.grid.major.y = element_blank(), # Hide y-axis grid lines
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12),  # Bold and size for treatment labels
        axis.title.x = element_text(face = "bold", size = 14)) +
  scale_x_discrete(position = "top") +  # Place treatment headers at the top
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))  # Evenly space species ranks and center them

