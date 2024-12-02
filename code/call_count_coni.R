#Clear Workspace
rm(list = ls())

#Load packages
library(tidyverse)
library(ggplot2)
library(viridis)
library(segmented)
library(mgcv)


#load count data 
###Data set before adding July/August/September 2024 data
call_count_coni <- read.csv('data/count_data/count_data_24_coni.csv')


# List all call count dataframes
call_count_dfs <- list(call_count_coni)

# Define the treatment assignment function
assign_treatment <- function(df) {
  df %>%
    mutate(treatment = case_when(
      str_detect(site, "^M") ~ "Mine",
      str_detect(site, "^T") ~ "Timber",
      str_detect(site, "^R") ~ "Rx Fire (Young)",
      str_detect(site, "^O") ~ "Rx Fire (Sec Growth)",
      TRUE ~ "Unknown"
    ))
}

# Apply the function to each dataframe in the list
call_count_dfs <- lapply(call_count_dfs, assign_treatment)

#filter out certain sites by site name 
call_count_dfs <- lapply(call_count_dfs, function(df) {
  df %>%
    filter(!site %in% c("M-6", "M-9", "M-11", "T-16", "T-21", "T-22", "R-1", "R-2", "R-6", "R-9", "R-11"))
})


#Remove observations past 2024-06-26
call_count_dfs <- lapply(call_count_dfs, function(df) {
  df %>%
    filter(date <= as.Date("2024-06-26"))
})

# Re-assign each modified dataframe back to the original name
call_count_coni <- call_count_dfs[[1]]

#turn date into date format
call_count_coni$date <- as.Date(call_count_coni$date)



# Create a sequence of dates for only the 1st and 15th of each month
date_breaks_custom <- seq(from = as.Date("2024-01-01"), to = as.Date("2024-12-31"), by = "1 day")
date_breaks_custom <- date_breaks_custom[format(date_breaks_custom, "%d") %in% c("01", "15")]


###Fitting with different model types###

# Define the segmented regression fitting function
fit_segmented <- function(data) {
  data <- data %>% mutate(date_numeric = as.numeric(date))
  lm_fit <- lm(count ~ date_numeric, data = data) # Fit initial linear model
  tryCatch({
    segmented_fit <- segmented(lm_fit, seg.Z = ~date_numeric, npsi = 2) # Add breakpoints
    data.frame(date = data$date, count = data$count, fit = fitted(segmented_fit))
  }, error = function(e) {
    message("Error in segmented fitting: ", e)
    NULL
  })
}

# Define a GAM fitting function with a smooth term for date
fit_gam <- function(data) {
  gam_fit <- gam(count ~ s(as.numeric(date)), data = data)
  data.frame(date = data$date, count = data$count, fit = fitted(gam_fit))
}


#CONI

# Apply segmented regression to each treatment group
segmented_fits <- call_count_coni %>%
  group_by(treatment) %>%
  group_modify(~ fit_segmented(.))

# Filter out any NULL fits
segmented_fits <- segmented_fits %>% filter(!is.null(fit))

# Plot the results
ggplot(call_count_coni, aes(x = date, y = count, color = treatment)) +
  geom_point(alpha = 0.3) +
  geom_line(data = segmented_fits, aes(x = date, y = fit, color = treatment), size = 1.2) +
  labs(title = "CONI - Breakpoint Regression Call Count for Treatments", x = "Date", y = "Call Count", color = "Site Treatment") +
  scale_x_date(date_labels = "%b-%d") +
  scale_colour_viridis_d(option = "viridis") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold"))


# Apply the GAM model to each treatment group
gam_fits <- call_count_coni %>%
  group_by(treatment) %>%
  group_modify(~ fit_gam(.))

# Plot the raw data with smooth GAM fits
ggplot(call_count_coni, aes(x = date, y = count, color = treatment)) +
  geom_point(alpha = 0.3) +
  geom_line(data = gam_fits, aes(x = date, y = fit, color = treatment), size = 1.2) +
  labs(title = "CONI - GAM Call Count for Treatments", x = "Date", y = "Call Count", color = "Site Treatment") +
  scale_x_date(date_labels = "%b-%d") +
  scale_colour_viridis_d(option = "viridis") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold"))



count_data <- read.csv("data/count_data_expanded_24_coni.csv")


