## Code to clean and format line intercept vegetation data for use as covariates 
## Collected across 36 sites/ARUs in 4 site treatments (M - Mine Reclamation, T - Timber Production, R - Young Savanna = , O = Mature Savanna) across GA
## Each site had four replicates of a 20 meter (m) line intercept established around each ARU, 
## Summing to a 80 m total line intercept sample (unless otherwise noted in site_adjustments) 
## Start and end lengths for vegetation was measured in m and height in decimeters (dm). 
## Last edited: 2025-02-13

# Clear Workspace  
rm(list=ls())

# Load packages 
library(tidyverse)

# Read in data
veg <- read.csv("data/line_intercept_veg_data_24.csv")

# Accounting for sites with reduced survey lengths 
site_adjustments <- data.frame(
  site = c("R-2", "M-12", "R-11"),
  survey_length = c(76, 40, 70) # The actual survey length for these particular sites
)

# Calculate total coverage by grouping and summing the line lengths
cover_type_summary <- veg %>%
  group_by(site, cover_type) %>%
  summarise(total_coverage = sum(line_length, na.rm = TRUE), .groups = 'drop') %>%
  left_join(site_adjustments, by = "site") %>%
  mutate(survey_length = ifelse(is.na(survey_length), 80, survey_length)) %>%
  mutate(percent_cover = (total_coverage / survey_length) * 100)

# Convert in wide format
cover_type_summary <- cover_type_summary %>%
  select(-total_coverage) %>%
  group_by(site) %>%
  pivot_wider(names_from = cover_type, values_from = percent_cover, values_fill = 0, names_glue = "{cover_type}_cover")

# Add columns for average heights of each cover type
height_cover_type_summary <- veg %>%
  group_by(site, cover_type) %>%
  summarise(mean_height = mean(height, na.rm = TRUE))

# Convert in wide format
height_cover_type_summary <- height_cover_type_summary %>%
  group_by(site) %>%
  pivot_wider(names_from = cover_type, values_from = mean_height, values_fill = 0, names_glue = "{cover_type}_height")

# Repeat for cover_type_sub
cover_type_sub_summary <- veg %>%
  group_by(site, cover_type_sub) %>%
  summarise(total_coverage = sum(line_length, na.rm = TRUE), .groups = 'drop') %>%
  left_join(site_adjustments, by = "site") %>%
  mutate(survey_length = ifelse(is.na(survey_length), 80, survey_length)) %>%
  mutate(percent_cover = (total_coverage / survey_length) * 100)

# Make new columns for cover in wide format
cover_type_sub_summary <- cover_type_sub_summary %>%
  select(-total_coverage) %>%
  group_by(site) %>%
  pivot_wider(names_from = cover_type_sub, values_from = percent_cover, values_fill = 0, names_glue = "{cover_type_sub}_cover")

# Add columns for average heights of each cover sub type
height_cover_type_sub_summary <- veg %>%
  group_by(site, cover_type_sub) %>%
  summarise(mean_height = mean(height, na.rm = TRUE))

# Make new columns for height in wide format
height_cover_type_sub_summary <- height_cover_type_sub_summary %>%
  group_by(site) %>%
  pivot_wider(names_from = cover_type_sub, values_from = mean_height, values_fill = 0, names_glue = "{cover_type_sub}_height")

# Calculate overall vegetation coverage by summing 'line_length' across all cover_types for each site
veg_cover_overall_summary <- veg %>%
  group_by(site) %>%
  summarise(sum_veg_cover = sum(line_length, na.rm = TRUE), .groups = 'drop') %>%
  left_join(site_adjustments, by = "site") %>%
  mutate(survey_length = ifelse(is.na(survey_length), 80, survey_length)) %>%
  mutate(veg_cover_overall = (sum_veg_cover / survey_length) * 100) %>%
  select(site, veg_cover_overall)

# Calculate overall vegetation height by averaging 'height' across all cover_types for each site- exclude conifer
veg_height_overall_summary <- veg %>%
  filter(cover_type != "conifer") %>%
  group_by(site) %>%
  summarise(veg_height_overall = mean(height, na.rm = TRUE)) %>%
  select(site, veg_height_overall)

# Calculate overall vegetation diversity by counting the number of unique cover_type_sub for each site
veg_diversity_overall <- veg %>%
  group_by(site) %>%
  summarise(veg_diversity_overall= n_distinct(cover_type_sub)) %>%
  select(site, veg_diversity_overall)

# Combine all summaries into one dataframe
line_intercept_summary <- left_join(cover_type_summary, height_cover_type_summary, by = "site") %>%
  left_join(cover_type_sub_summary, by = "site") %>%
  left_join(height_cover_type_sub_summary, by = "site") %>%
  left_join(veg_cover_overall_summary, by = "site") %>%
  left_join(veg_height_overall_summary, by = "site") %>%
  left_join(veg_diversity_overall, by = "site")
  
# Remove survey_length.y and survey_length.x column
line_intercept_summary <- line_intercept_summary %>%
  select(-survey_length.x, -survey_length.y) 


# Create a new column for the site treatment grouped by the first letter of the site name
line_intercept_summary <- line_intercept_summary %>%
  mutate(treatment = case_when(
    str_detect(site, "^M") ~ "mine",
    str_detect(site, "^T") ~ "timber",
    str_detect(site, "^R") ~ "young_savanna",
    str_detect(site, "^O") ~ "mature_savanna",
    TRUE ~ "Unknown"
  ))

#save
write.csv(line_intercept_summary, "data/line_intercept_summary.csv", row.names = FALSE)
