## Dan Watson
## Last edit: 08/02/2024
## Code to clean and format line intercept vegetation data for use as covariates in occupancy models 
## Collected across 36 sites/ARUs in 4 site treatments (T = Timber, M = Mine, R = Sansavilla, O = Okefenokee) across GA
## Each site had four replicates of a 20 meter (m) line intercept established around each ARU, 
## Summing to a 80 m total line intercept sample (unless otherwise noted in site_adjustments) 
## Start and end lengths for vegetation was measured in m and height in decimeters (dm). 

# Clear workspace  
rm(list=ls())

# Load packages 
library(tidyverse)

# Read in data, "data/line_intercept_veg_data_24_raw.csv"
raw_veg <- read.csv("data/line_intercept_veg_data_24_raw.csv")

#check data structure
str(raw_veg)

#format column structure correctly 
raw_veg$start <- as.numeric(raw_veg$start)
raw_veg$end <- as.numeric(raw_veg$end)
raw_veg$height <- as.numeric(raw_veg$height)
raw_veg$date <- as.Date(raw_veg$date)

## Checking data for errors

# Identify every lowest 'start' value for every unique 'site_line'
raw_veg_min <- raw_veg %>% 
  group_by(site_line) %>% 
  filter(start == min(start))

# Identify every highest 'end' value for every unique 'site_line'
raw_veg_max <- raw_veg %>% 
  group_by(site_line) %>% 
  filter(end == max(end))

# Identify unique values for site in a df
unique_site <- unique(raw_veg_min$site)

# Show unique values for site_line
unique_site_line <- unique(raw_veg_min$site_line)

# Clean up the 'cover_type' column 
raw_veg <- raw_veg %>%
  mutate(cover_type = case_when(
    cover_type == "shrubs" ~ "shrub",
    cover_type == "grasses" ~ "grass",
    cover_type == "forbs" ~ "forb",
    cover_type == " conifer" ~ "conifer",
    cover_type == "confier" ~ "conifer",
    cover_type == "forb " ~ "forb",
    cover_type == "forb/grass" ~ "forb",
    cover_type == "grass/forb" ~ "grass",
    TRUE ~ cover_type
  ))
    
#Display values for unique 'cover_type'
unique_cover_type <- unique(raw_veg$cover_type)
unique_cover_type

# Clean up the 'cover_type_sub' column
raw_veg <- raw_veg %>%
  mutate(cover_type_sub = case_when(
    cover_type_sub == "saw" ~ "saw_palmetto",
    cover_type_sub == "broadleaf_shrubs" ~ "broadleaf_shrub",
    cover_type_sub == "broadleaf_shrub " ~ "broadleaf_shrub",
    cover_type_sub == "shrub" ~ "broadleaf_shrub",
    cover_type_sub == "ferns" ~ "fern",
    cover_type_sub == "forb_fern" ~ "fern",
    cover_type_sub == "forbs_other" ~ "forb_other",
    cover_type_sub == "forn_other" ~ "forb_other",
    cover_type_sub == "forb_other " ~ "forb_other",
    cover_type_sub == "bluestem " ~ "bluestem",
    cover_type_sub == "forb_other/grass_other" ~ "forb_other",
    cover_type_sub == "non.native" ~ "non_native",
    cover_type_sub == "non-native" ~ "non_native",
    # Add more conditions if there are other spelling variations
    TRUE ~ cover_type_sub
  ))

# Display values for unique 'cover_type_sub' 
unique_cover_type_sub <- unique(raw_veg$cover_type_sub)
unique_cover_type_sub

# Create a new column for line length and check values 
raw_veg <- raw_veg %>%
  group_by(site_line) %>%
  mutate(line_length = end - start)

# Once values are checked for mistakes, remove any accidental observations with a line_length < .19
raw_veg <- raw_veg %>%
  filter(line_length > .19)

# Good point to save work and clean environment
# Write the data to a new csv file
write.csv(raw_veg, "data/line_intercept_veg_data_24_clean.csv", row.names = FALSE)
# Clear workspace  
rm(list=ls())

# Read in the cleaned data
veg <- read.csv("data/line_intercept_veg_data_24_clean.csv")

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

 #Convert in wide format
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

# Write the data to a new csv file
write.csv(line_intercept_summary, "data/line_intercept_summary.csv", row.names = FALSE)

# Read in 
line_intercept_summary <- read.csv("data/line_intercept_summary.csv")


# Create a new column for the site treatment grouped by the first letter of the site name
line_intercept_summary <- line_intercept_summary %>%
  mutate(treatment = case_when(
    str_detect(site, "^M") ~ "mine",
    str_detect(site, "^T") ~ "timber",
    str_detect(site, "^R") ~ "rx_fire_young ",
    str_detect(site, "^O") ~ "rx_fire_sec_growth",
    TRUE ~ "Unknown"
  ))

# Calculate average cover types by treatment just to look at
line_intercept_summary %>%
  group_by(treatment) %>%
  summarise(
    mean_shrub_cover = mean(shrub_cover),
    mean_forb_cover = mean(forb_cover),
    mean_grass_cover = mean(grass_cover),
    mean_conifer_cover = mean(conifer_cover),
    mean_veg_cover_overall = mean(veg_cover_overall),
    mean_veg_height_overall = mean(veg_height_overall),
    mean_veg_diversity_overall = mean(veg_diversity_overall)
  )


#save
write.csv(line_intercept_summary, "data/line_intercept_summary.csv", row.names = FALSE)