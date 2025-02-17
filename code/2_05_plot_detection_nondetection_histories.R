#Plot detection / non-detection histories from count data
#Last updated 2025-02-13

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(ggplot2)
library(unmarked)
library(cowplot)

# Load count data for Blue Grosbeak
count_blgr <- read.csv("data/count_data/count_data_24_blgr.csv")

# Ensure date column is in Date format
count_blgr$date <- as.Date(count_blgr$date)

# Restrict to 2024-05-03 - 2024-06-21
count_blgr <- count_blgr %>%
  filter(date >= as.Date("2024-05-03") & date <= as.Date("2024-06-21"))


# Define your desired date range
start_date <- as.Date("2024-05-03")
end_date <- as.Date("2024-06-21")

# Generate daily sequence of dates
date_seq <- seq.Date(from = start_date, to = end_date, by = "1 day")

# Filter in 5-day intervals
filtered_dates <- date_seq[seq(1, length(date_seq), by = 5)]

# Heatmap for Blue Grosbeak call count indices 
heatmap_blgr_call_count <- ggplot(count_blgr, aes(x = date, y = site, fill = count)) +
  geom_tile(color = "black", size = 0.25) +
  labs(title = "Blue Grosbeak", x = "Date", y = "Monitoring Point", fill = "Call Counts") +
  scale_fill_gradientn(
    colors = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"), 
    values = scales::rescale(c(0, 1, 2, 5, max(count_blgr$count))),
    name = "Detections"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, face = "bold", size = 10),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, face = "bold", size = 10),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  scale_x_date(breaks = filtered_dates, date_labels = "%m-%d") +
  guides(fill = guide_legend(reverse = TRUE))

# Print the heatmap
print(heatmap_blgr_call_count)

# Convert values greater than 0 or NA to 1 in "count_blgr"
count_blgr$count <- ifelse(count_blgr$count > 0 | is.na(count_blgr$count), 1, 0)

# Heatmap of just detection/non-detection
heatmap_blgr <- ggplot(count_blgr, aes(x = date, y = site, fill = as.factor(count))) +
  geom_tile(color = "black", size = 0.25) +
  scale_fill_manual(
    values = c("0" = "white", "1" = "skyblue"),
    labels = c("0" = "Non-Detection", "1" = "Detection"),
    name = "Detection History"
  ) +
  labs(title = "Blue Grosbeak", x = "Date", y = "Monitoring Point", fill = "History") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, face = "bold", size = 10),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, face = "bold", size = 10),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  scale_x_date(breaks = filtered_dates, date_labels = "%m-%d") +
  guides(fill = guide_legend(reverse = TRUE))

# Print the heatmap
print(heatmap_blgr)

# Repeat for Bachman's Sparrow
count_bacs <- read.csv("data/count_data/count_data_24_bacs.csv")

# Ensure date column is in Date format
count_bacs$date <- as.Date(count_bacs$date)

# Convert values greater than 0 or NA to 1 in "count_bacs"
count_bacs$count <- ifelse(count_bacs$count > 0 | is.na(count_bacs$count), 1, 0)

# Restrict to 2024-05-03 - 2024-06-21
count_bacs <- count_bacs %>%
  filter(date >= as.Date("2024-05-03") & date <= as.Date("2024-06-21"))

# Define your desired date range
start_date <- as.Date("2024-05-03")
end_date <- as.Date("2024-06-21")

# Generate daily sequence of dates
date_seq <- seq.Date(from = start_date, to = end_date, by = "1 day")

# Filter in 5-day intervals
filtered_dates <- date_seq[seq(1, length(date_seq), by = 5)]

# Heatmap for Bachman's Sparrow
heatmap_bacs <- ggplot(count_bacs, aes(x = date, y = site, fill = as.factor(count))) +
  geom_tile(color = "black", size = 0.25) +
  scale_fill_manual(
    values = c("0" = "white", "1" = "skyblue"),
    labels = c("0" = "No Detection", "1" = "Detection")
  ) +
  labs(title = "Bachman's Sparrow", x = "Date", y = "Monitoring Point", fill = "History") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, face = "bold", size = 10),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, face = "bold", size = 10),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  scale_x_date(breaks = filtered_dates, date_labels = "%m-%d") +
  guides(fill = guide_legend(reverse = TRUE))

# Print the heatmap
print(heatmap_bacs)

# Shared y-axis label
y_label <- ggdraw() + 
  draw_label("Monitoring Point", fontface = "bold", angle = 90, size = 14)

# Shared legend and scale
common_scale <- scale_fill_manual(
  values = c("0" = "white", "1" = "darkgrey"),
  labels = c("0" = "Non-Detection", "1" = "Detection"),
  name = "History" # Legend title
)

# Blue Grosbeak heatmap
heatmap_blgr <- ggplot(count_blgr, aes(x = date, y = site, fill = as.factor(count))) +
  geom_tile(color = "black", size = 0.25) +
  common_scale +
  labs(title = "Blue Grosbeak", x = NULL, y = NULL) +
  theme_minimal_hgrid() +
  theme(
    axis.text.x = element_blank(), # Remove x-axis text
    axis.ticks.x = element_line(color = "grey", size = 0.5), 
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, face = "bold", size = 12),
    legend.position = "none", # Remove legend for this plot
    plot.title = element_text(face = "bold", size = 14)
  ) +
  scale_x_date(breaks = filtered_dates, date_labels = "%m-%d")

# Bachman's Sparrow heatmap
heatmap_bacs <- ggplot(count_bacs, aes(x = date, y = site, fill = as.factor(count))) +
  geom_tile(color = "black", size = 0.25) +
  common_scale +
  labs(title = "Bachman's Sparrow", x = NULL, y = NULL, fill = NULL) +
  theme_minimal_hgrid() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, face = "bold", size = 12),
    axis.ticks.x = element_line(color = "grey", size = 0.5), 
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right", # Legend only for this plot
    plot.title = element_text(face = "bold", size = 14)
  ) +
  scale_x_date(breaks = filtered_dates, date_labels = "%m-%d") +
  guides(fill = guide_legend(reverse = TRUE))

# Extract legend from Bachman's Sparrow plot
legend <- get_legend(heatmap_bacs)

# Remove legend from Bachman's Sparrow plot
heatmap_bacs <- heatmap_bacs +
  theme(legend.position = "none")

# Combine plots using cowplot
combined_plots <- plot_grid(
  heatmap_blgr, heatmap_bacs, 
  ncol = 1, # Stack vertically
  align = "v",
  rel_heights = c(1, 1) # Equal height
)

# Add shared y-axis label and legend
final_plot <- plot_grid(
  y_label, combined_plots, legend,
  ncol = 3, # Arrange in 3 columns
  rel_widths = c(0.03, 1, 0.15) # Adjust widths for y-axis, plots, and legend
)

# Display the final plot
print(final_plot)
