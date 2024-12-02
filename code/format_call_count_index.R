# Clear workspace 
rm(list=ls()) 

# Load packages 
library(tidyverse)
library(ggplot2)

# Read in acoustic_dat_24_07.csv
aru_data_master <- read.csv('~/Downloads/acoustic_dat_24_07.csv')

# Read in probabilistic_thresholds_master.csv
p_thresholds <- read.csv('data/thresholds_table_master.csv')

# Filter out dates before 2024-03-10 and after 2024-06-28
aru_data <- aru_data_master %>% 
  filter(date >= '2024-03-10' & date <= '2024-06-28')

#Remove sites labelled as 'NA' - remnants from having to switch ARUs in the field
aru_data <- aru_data %>% filter(site != 'NA')

# Format date as date
aru_data$date <- as.Date(aru_data$date)

# Create a sequence of all dates and extract all unique sites
all_dates <- seq.Date(from = min(aru_data$date), to = max(aru_data$date), by = "day")
all_sites <- unique(aru_data$site)

# Filter down to a focal species
aru_data <- aru_data %>% filter(sp_code == 'WEVI')

# Filter out any rows at times >13:00:00
# aru_data <- aru_data %>% filter(time <= '13:00:00')

# Add a column for the logit score
aru_data$logit <- log(aru_data$confidence/(1-aru_data$confidence))

# Create column for count == 1
aru_data <- aru_data %>% mutate(count = 1)

# Remove observations with a confidence score or logit score 
#aru_data <- aru_data %>% filter(confidence >= 0.7731874)
aru_data <- aru_data %>% filter(logit >= 0.88569517)

# Create call_count_long with complete date and site combinations filled with zeros
call_count_long <- aru_data %>%
  group_by(site, date) %>%
  summarise(count = n(), .groups = 'drop') %>%
  complete(site = all_sites, date = all_dates, fill = list(count = 0)) 

# Save call_count_long to csv
write.csv(call_count_long, 'data/count_data_24_wevi.csv', row.names = FALSE)

# Now pivot call_count_long to wide format
call_count_wide <- call_count_long %>%
  pivot_wider(names_from = date, values_from = count, values_fill = list(count = 0)) 

# Change date from yyyy-mm-dd to mm-dd
call_count_long$date <- format(as.Date(call_count_long$date), "%m-%d")

# visualize heatmap
ggplot(call_count_long, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "WEVI Call Count Across Sites and Dates (99% confidence threshold)", x = "Date", y = "Site", fill = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# visualize line plot
#ggplot(call_count_long, aes(x = date, y = count, group = site, color = site)) +
#  geom_line() +
#  labs(title = "Counts Across Sites and Dates", x = "Date", y = "Count", color = "Site") +
#  theme_minimal() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#load count data 

count_data <- read.csv('data/count_data/count_data_24_cwwi.csv')



# Filter for only site "O-2"
count_data <- count_data %>% filter(site == "R-12")

#Plot as points for the call count for site "R-12" over date 
ggplot(count_data, aes(x = date, y = count)) +
  geom_point() +
  labs(title = "Call Count for Site R-12", x = "Date", y = "Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

















