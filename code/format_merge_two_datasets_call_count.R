#Combining ARU datasets 
#Clear Workspace
rm(list=ls())

#Load packages
library(tidyverse)

#load datasets 
# Read in acoustic_dat_24_07.csv
aru_data_24_07 <- read.csv('~/Downloads/acoustic_dat_24_07.csv')
# Read in acoustic_dat_24_09.csv
aru_data_24_09 <- read.csv('~/Downloads/acoustic_dat_24_09.csv')

#remove "X", "Y", "utm_zone_3" columns from 07
aru_data_24_07 <- aru_data_24_07 %>% select(-X, -Y, -utm_zone_3)

#merge the two datasets
aru_data <- rbind(aru_data_24_07, aru_data_24_09)

#Remove sites labelled as 'NA' - remnants from having to switch ARUs in the field
aru_data <- aru_data %>% filter(site != 'NA')

#write 
write.csv(aru_data, '~/Downloads/aoustic_dat.csv')



#read in the data
aru_data <- read.csv('~/Downloads/aoustic_dat.csv')
# Read in timeline 
aru_timeline <- read_csv("data/aru_timeline.csv")
# load thresholds 
aru_thresholds <- read_csv("data/thresholds_table_master.csv")

# Format date as date
aru_data$date <- as.Date(aru_data$date)

# Create a sequence of all dates and extract all unique sites
all_dates <- seq.Date(from = min(aru_data$date), to = max(aru_data$date), by = "day")
all_sites <- unique(aru_data$site)

# Filter down to a focal species
aru_data <- aru_data %>% filter(sp_code == 'NOBO')

# Add a column for the logit score
aru_data$logit <- log(aru_data$confidence/(1-aru_data$confidence))

# Create column for count == 1
aru_data <- aru_data %>% mutate(count = 1)

# Remove observations with a confidence score or logit score 
aru_data <- aru_data %>% filter(confidence >= 0.7731874)
#aru_data <- aru_data %>% filter(logit >= -2.09271153)

# Create call_count_long with complete date and site combinations filled with zeros
call_count_long <- aru_data %>%
  group_by(site, date) %>%
  summarise(count = n(), .groups = 'drop') %>%
  complete(site = all_sites, date = all_dates, fill = list(count = 0)) 

# Now pivot call_count_long to wide format
call_count_wide <- call_count_long %>%
  pivot_wider(names_from = date, values_from = count, values_fill = list(count = 0)) 

# visualize heatmap
ggplot(call_count_long, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "NOBO Call Count Across Sites and Dates (99% confidence threshold)", x = "Date", y = "Site", fill = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Save call_count_long to csv
write.csv(call_count_long, 'data/count_data_24_nobo_phenology.csv', row.names = FALSE)














