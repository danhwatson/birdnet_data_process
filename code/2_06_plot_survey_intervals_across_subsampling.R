#Plots survey intervals across subsampling 
#Last updated 2025-02-12


#clear workspace
rm(list=ls())

#load libraries
library(tidyverse)
library(ggplot2)
library(viridis)

# Read the CSV file
df <- read.csv("data/survey_intervals_across_subsampling.csv")

# Transform the dataframe to long format
df_long <- pivot_longer(df, cols = -dataset, names_to = "date", values_to = "survey")

# Fix the date format from 'XYYYY.MM.DD' to 'YYYY-MM-DD'
df_long <- df_long %>%
  mutate(date = gsub("X", "", date)) %>%
  mutate(date = gsub("\\.", "-", date))

# Replace NA with "NA" for 'survey' column
df_long$survey[is.na(df_long$survey)] <- "NA"


# Replace NA with 0 for 'survey' column
df_long$survey[is.na(df_long$survey)] <- 0

# Custom order for datasets (flipped order)
df_long$dataset <- factor(df_long$dataset, levels = c("M-1_every_10_days", "M-1_every_7_days", "M-1_every_4_days", "M-1_every_3_days", "M-1_every_2_days", "M-1_every_day"))

# Define your desired date range
start_date <- as.Date("2024-05-03")
end_date <- as.Date("2024-06-21")

# Generate daily sequence of dates
date_seq <- seq.Date(from = start_date, to = end_date, by = "1 day")

# Filter in 5-day intervals
filtered_dates <- date_seq[seq(1, length(date_seq), by = 5)]

#ensure date is in date format 
df_long$date <- as.Date(df_long$date)


# Create the heatmap
ggplot(df_long, aes(x = date, y = dataset, fill = factor(survey))) +
  geom_tile(color = "black", size = .25) +
  scale_fill_manual(values = c("darkgrey", "white"),
                    labels = c("NA" = "No Surveys", "1" = "Surveys")) +
       labs(title = "",
       x = "Date",
       y = "Dataset",
       fill = "Schedule") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12, face = "bold"),
        aspect.ratio = 1) +
  scale_y_discrete(labels = c("M-1_every_day" = "Everyday",
                              "M-1_every_2_days" = "Every 2nd Day",
                              "M-1_every_3_days" = "Every 3rd Day",
                              "M-1_every_4_days" = "Every 4th Day",
                              "M-1_every_7_days" = "Every 7th Day",
                              "M-1_every_10_days" = "Every 10th Day")) +
  scale_x_date(breaks = filtered_dates, date_labels = "%m-%d")


# Display the first few rows of the long dataframe
head(df_long)
