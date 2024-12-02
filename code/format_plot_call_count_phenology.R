#clear workspace 
rm(list=ls())

#load packages 
library(tidyverse)
library(ggplot2)
library(viridis)
library(segmented)
library(mgcv)

#load data
count_bacs <- read.csv("data/phenology_24_bacs.csv")
count_cwwi <- read.csv("data/phenology_24_cwwi.csv")
count_nobo <- read.csv("data/phenology_24_nobo.csv")
count_coni <- read.csv("data/phenology_24_coni.csv")

#ensure date is dat format 
count_bacs$date <- as.Date(count_bacs$date)
count_cwwi$date <- as.Date(count_cwwi$date)
count_nobo$date <- as.Date(count_nobo$date)
count_coni$date <- as.Date(count_coni$date)

#Add treatments 

count_bacs <- count_bacs %>% 
  mutate(treatment = case_when(
    str_detect(site, "^M") ~ "Mine",
    str_detect(site, "^T") ~ "Timber",
    str_detect(site, "^R") ~ "Rx Fire (Young)",
    str_detect(site, "^O") ~ "Rx Fire (Sec Growth)",
    TRUE ~ "Unknown"
  ))

count_cwwi <- count_cwwi %>% 
  mutate(treatment = case_when(
    str_detect(site, "^M") ~ "Mine",
    str_detect(site, "^T") ~ "Timber",
    str_detect(site, "^R") ~ "Rx Fire (Young)",
    str_detect(site, "^O") ~ "Rx Fire (Sec Growth)",
    TRUE ~ "Unknown"
  ))

count_nobo <- count_nobo %>% 
  mutate(treatment = case_when(
    str_detect(site, "^M") ~ "Mine",
    str_detect(site, "^T") ~ "Timber",
    str_detect(site, "^R") ~ "Rx Fire (Young)",
    str_detect(site, "^O") ~ "Rx Fire (Sec Growth)",
    TRUE ~ "Unknown"
  ))

count_coni <- count_coni %>% 
  mutate(treatment = case_when(
    str_detect(site, "^M") ~ "Mine",
    str_detect(site, "^T") ~ "Timber",
    str_detect(site, "^R") ~ "Rx Fire (Young)",
    str_detect(site, "^O") ~ "Rx Fire (Sec Growth)",
    TRUE ~ "Unknown"
  ))


#remove observations before March 15th 
count_bacs <- count_bacs %>% filter(date >= "2024-03-15")
count_cwwi <- count_cwwi %>% filter(date >= "2024-03-15")
count_nobo <- count_nobo %>% filter(date >= "2024-03-15")
count_coni <- count_coni %>% filter(date >= "2024-03-15")
#remove observations after August 1st
count_bacs <- count_bacs %>% filter(date <= "2024-08-01")
count_cwwi <- count_cwwi %>% filter(date <= "2024-08-01")
count_nobo <- count_nobo %>% filter(date <= "2024-08-01")
count_coni <- count_coni %>% filter(date <= "2024-08-01")


# Create new dfs seperated out by treatment type 
# BACS
count_bacs_mine <- count_bacs %>% filter(treatment == "Mine")
count_bacs_timber <- count_bacs %>% filter(treatment == "Timber")
count_bacs_rx_fire_young <- count_bacs %>% filter(treatment == "Rx Fire (Young)")
count_bacs_rx_fire_sec_growth <- count_bacs %>% filter(treatment == "Rx Fire (Sec Growth)")

# CWWI
count_cwwi_mine <- count_cwwi %>% filter(treatment == "Mine")
count_cwwi_timber <- count_cwwi %>% filter(treatment == "Timber")
count_cwwi_rx_fire_young <- count_cwwi %>% filter(treatment == "Rx Fire (Young)")
count_cwwi_rx_fire_sec_growth <- count_cwwi %>% filter(treatment == "Rx Fire (Sec Growth)")


# NOBO
count_nobo_mine <- count_nobo %>% filter(treatment == "Mine")
count_nobo_timber <- count_nobo %>% filter(treatment == "Timber")
count_nobo_rx_fire_young <- count_nobo %>% filter(treatment == "Rx Fire (Young)")
count_nobo_rx_fire_sec_growth <- count_nobo %>% filter(treatment == "Rx Fire (Sec Growth)")

# CONI
count_coni_mine <- count_coni %>% filter(treatment == "Mine")
count_coni_timber <- count_coni %>% filter(treatment == "Timber")
count_coni_rx_fire_young <- count_coni %>% filter(treatment == "Rx Fire (Young)")
count_coni_rx_fire_sec_growth <- count_coni %>% filter(treatment == "Rx Fire (Sec Growth)")



# Create a sequence of dates for only the 1st and 15th of each month
date_breaks_custom <- seq(from = as.Date("2024-01-01"), to = as.Date("2024-12-31"), by = "1 day")
date_breaks_custom <- date_breaks_custom[format(date_breaks_custom, "%d") %in% c("01", "15")]



#BACS heatmaps
heatmap_bacs_mine <- ggplot(count_bacs_mine, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "BACS - Mine", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_bacs_mine$site)) - 0.5, color = "grey80", linewidth = 0.3)


heatmap_bacs_mine

#Remove sites M-9, M-6, M-5, M-2, M-12, M-11, M-10x
count_bacs_mine <- count_bacs_mine %>% filter(site != "M-9" & site != "M-6" & site != "M-5" & 
                                                site != "M-2" & site != "M-12" & site != "M-11" & site != "M-10")



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

segmented_fits_bacs_mine <- count_bacs_mine %>% 
  group_by(site) %>%
  group_modify(~ fit_segmented(.))

#identify peak in for each site 
peak_date_bacs_mine <- segmented_fits_bacs_mine %>% 
  group_by(site) %>%
  filter(fit == max(fit))

#plot segmented regression
ggplot(count_bacs_mine, aes(x = date, y = count, color = site)) +
  geom_point(alpha = 0.3) +
  geom_line(data = segmented_fits_bacs_mine, aes(x = date, y = fit, color = site), size = 1.2) +
  labs(title = "BACS - Mine - Breakpoint Regression Call Count for Sites", x = "Date", y = "Call Count", color = "Site") +
  scale_x_date(date_labels = "%b-%d") +
  scale_colour_viridis_d(option = "viridis") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold"))


# Apply the GAM model 
gam_fits_bacs_mine <- count_bacs_mine %>% 
  group_by(site) %>%
  group_modify(~ fit_gam(.))

#plot GAM
ggplot(count_bacs_mine, aes(x = date, y = count, color = site)) +
  geom_point(alpha = 0.3) +
  geom_line(data = gam_fits_bacs_mine, aes(x = date, y = fit, color = site), size = 1.2) +
  labs(title = "BACS - Mine - GAM Call Count for Sites", x = "Date", y = "Call Count", color = "Site") +
  scale_x_date(date_labels = "%b-%d") +
  scale_colour_viridis_d(option = "viridis") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold"))


heatmap_bacs_timber <- ggplot(count_bacs_timber, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "BACS - Timber", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_bacs_timber$site)) - 0.5, color = "grey80", linewidth = 0.3)

heatmap_bacs_timber

#remove observations after June 27th
count_bacs_timber <- count_bacs_timber %>% filter(date <= "2024-06-27")



heatmap_bacs_rx_fire_young <- ggplot(count_bacs_rx_fire_young, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "BACS - Rx Fire Young", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_bacs_rx_fire_young$site)) - 0.5, color = "grey80", linewidth = 0.3)

heatmap_bacs_rx_fire_young

#Remove sites 
count_bacs_rx_fire_young <- count_bacs_rx_fire_young %>% filter(site != "R-9" & site != "R-6" & site != "R-2" & 
                                                site != "R-13" & site != "R-11" & site != "R-10", site != "R-1")

#apply the segmented regression model 
segmented_fits_bacs_rx_fire_young <- count_bacs_rx_fire_young %>% 
  group_by(site) %>%
  group_modify(~ fit_segmented(.))

#identify peak in for each site
peak_date_bacs_rx_fire_young <- segmented_fits_bacs_rx_fire_young %>% 
  group_by(site) %>%
  filter(fit == max(fit))

#plot segmented regression
ggplot(count_bacs_rx_fire_young, aes(x = date, y = count, color = site)) +
  geom_point(alpha = 0.3) +
  geom_line(data = segmented_fits_bacs_rx_fire_young, aes(x = date, y = fit, color = site), size = 1.2) +
  labs(title = "BACS - Rx Fire Young - Breakpoint Regression Call Count for Sites", x = "Date", y = "Call Count", color = "Site") +
  scale_x_date(date_labels = "%b-%d") +
  scale_colour_viridis_d(option = "viridis") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold"))

# Apply the GAM model
gam_fits_bacs_rx_fire_young <- count_bacs_rx_fire_young %>% 
  group_by(site) %>%
  group_modify(~ fit_gam(.))


#plot GAM
ggplot(count_bacs_rx_fire_young, aes(x = date, y = count, color = site)) +
  geom_point(alpha = 0.3) +
  geom_line(data = gam_fits_bacs_rx_fire_young, aes(x = date, y = fit, color = site), size = 1.2) +
  labs(title = "BACS - Rx Fire Young - GAM Call Count for Sites", x = "Date", y = "Call Count", color = "Site") +
  scale_x_date(date_labels = "%b-%d") +
  scale_colour_viridis_d(option = "viridis") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold"))


heatmap_bacs_rx_fire_sec_growth <- ggplot(count_bacs_rx_fire_sec_growth, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "BACS - Rx Fire Sec Growth", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_bacs_rx_fire_sec_growth$site)) - 0.5, color = "grey80", linewidth = 0.3)

heatmap_bacs_rx_fire_sec_growth

#remove sites
count_bacs_rx_fire_sec_growth <- count_bacs_rx_fire_sec_growth %>% filter(site != "O-4" & site != "O-3")

#apply the segmented regression model
segmented_fits_bacs_rx_fire_sec_growth <- count_bacs_rx_fire_sec_growth %>% 
  group_by(site) %>%
  group_modify(~ fit_segmented(.))

#identify peak in for each site
peak_date_bacs_rx_fire_sec_growth <- segmented_fits_bacs_rx_fire_sec_growth %>% 
  group_by(site) %>%
  filter(fit == max(fit))

#plot segmented regression
ggplot(count_bacs_rx_fire_sec_growth, aes(x = date, y = count, color = site)) +
  geom_point(alpha = 0.3) +
  geom_line(data = segmented_fits_bacs_rx_fire_sec_growth, aes(x = date, y = fit, color = site), size = 1.2) +
  labs(title = "BACS - Rx Fire Sec Growth - Breakpoint Regression Call Count for Sites", x = "Date", y = "Call Count", color = "Site") +
  scale_x_date(date_labels = "%b-%d") +
  scale_colour_viridis_d(option = "viridis") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold"))

# Apply the GAM model
gam_fits_bacs_rx_fire_sec_growth <- count_bacs_rx_fire_sec_growth %>% 
  group_by(site) %>%
  group_modify(~ fit_gam(.))

#plot GAM
ggplot(count_bacs_rx_fire_sec_growth, aes(x = date, y = count, color = site)) +
  geom_point(alpha = 0.3) +
  geom_line(data = gam_fits_bacs_rx_fire_sec_growth, aes(x = date, y = fit, color = site), size = 1.2) +
  labs(title = "BACS - Rx Fire Sec Growth - GAM Call Count for Sites", x = "Date", y = "Call Count", color = "Site") +
  scale_x_date(date_labels = "%b-%d") +
  scale_colour_viridis_d(option = "viridis") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold"))

#CWWI heatmaps
heatmap_cwwi_mine <- ggplot(count_cwwi_mine, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "CWWI - Mine", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_cwwi_mine$site)) - 0.5, color = "grey80", linewidth = 0.3)


heatmap_cwwi_mine

#Remove sites 
count_cwwi_mine <- count_cwwi_mine %>% filter(site != "M-9" & site != "M-6" & site != "M-5" & 
                                              site != "M-11", site != "M-7")

#apply the segmented regression model 
segmented_fits_cwwi_mine <- count_cwwi_mine %>% 
  group_by(site) %>%
  group_modify(~ fit_segmented(.))

#identify peak in for each site
peak_date_cwwi_mine <- segmented_fits_cwwi_mine %>% 
  group_by(site) %>%
  filter(fit == max(fit))

#plot segmented regression
ggplot(count_cwwi_mine, aes(x = date, y = count, color = site)) +
  geom_point(alpha = 0.3) +
  geom_line(data = segmented_fits_cwwi_mine, aes(x = date, y = fit, color = site), size = 1.2) +
  labs(title = "CWWI - Mine - Breakpoint Regression Call Count for Sites", x = "Date", y = "Call Count", color = "Site") +
  scale_x_date(date_labels = "%b-%d") +
  scale_colour_viridis_d(option = "viridis") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold"))

# Apply the GAM model
gam_fits_cwwi_mine <- count_cwwi_mine %>% 
  group_by(site) %>%
  group_modify(~ fit_gam(.))

#plot GAM
ggplot(count_cwwi_mine, aes(x = date, y = count, color = site)) +
  geom_point(alpha = 0.3) +
  geom_line(data = gam_fits_cwwi_mine, aes(x = date, y = fit, color = site), size = 1.2) +
  labs(title = "CWWI - Mine - GAM Call Count for Sites", x = "Date", y = "Call Count", color = "Site") +
  scale_x_date(date_labels = "%b-%d") +
  scale_colour_viridis_d(option = "viridis") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold"))


heatmap_cwwi_timber <- ggplot(count_cwwi_timber, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "CWWI - Timber", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_cwwi_timber$site)) - 0.5, color = "grey80", linewidth = 0.3)

heatmap_cwwi_timber


heatmap_cwwi_rx_fire_young <- ggplot(count_cwwi_rx_fire_young, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "CWWI - Rx Fire Young", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_cwwi_rx_fire_young$site)) - 0.5, color = "grey80", linewidth = 0.3)

heatmap_cwwi_rx_fire_young

heatmap_cwwi_rx_fire_sec_growth <- ggplot(count_cwwi_rx_fire_sec_growth, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "CWWI - Rx Fire Sec Growth", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_cwwi_rx_fire_sec_growth$site)) - 0.5, color = "grey80", linewidth = 0.3)

heatmap_cwwi_rx_fire_sec_growth


#NOBO heatmaps
heatmap_nobo_mine <- ggplot(count_nobo_mine, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "NOBO - Mine", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_nobo_mine$site)) - 0.5, color = "grey80", linewidth = 0.3)


heatmap_nobo_mine


heatmap_nobo_timber <- ggplot(count_nobo_timber, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "NOBO - Timber", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_nobo_timber$site)) - 0.5, color = "grey80", linewidth = 0.3)

heatmap_nobo_timber


heatmap_nobo_rx_fire_young <- ggplot(count_nobo_rx_fire_young, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "NOBO - Rx Fire Young", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_nobo_rx_fire_young$site)) - 0.5, color = "grey80", linewidth = 0.3)

heatmap_nobo_rx_fire_young

heatmap_nobo_rx_fire_sec_growth <- ggplot(count_nobo_rx_fire_sec_growth, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "NOBO - RX Fire Sec Growth", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_nobo_rx_fire_sec_growth$site)) - 0.5, color = "grey80", linewidth = 0.3)

heatmap_nobo_rx_fire_sec_growth


#CONI heatmaps
heatmap_coni_mine <- ggplot(count_coni_mine, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "CONI - Mine", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_coni_mine$site)) - 0.5, color = "grey80", linewidth = 0.3)


heatmap_coni_mine


heatmap_coni_timber <- ggplot(count_coni_timber, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "CONI - Timber", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_coni_timber$site)) - 0.5, color = "grey80", linewidth = 0.3)

heatmap_coni_timber


heatmap_coni_rx_fire_young <- ggplot(count_coni_rx_fire_young, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "CONI - Rx Fire Young", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_coni_rx_fire_young$site)) - 0.5, color = "grey80", linewidth = 0.3)

heatmap_coni_rx_fire_young

heatmap_coni_rx_fire_sec_growth <- ggplot(count_coni_rx_fire_sec_growth, aes(x = date, y = site, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "CONI - Rx Fire Sec Growth", x = "Date", y = "Site", fill = "Count") +
  scale_x_date(breaks = date_breaks_custom, date_labels = "%b-%d") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, face = "bold", size = 10),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, face = "bold", size = 10),
        legend.title = element_text(face = "bold")) +
  geom_hline(yintercept = seq_along(unique(count_coni_rx_fire_sec_growth$site)) - 0.5, color = "grey80", linewidth = 0.3)

heatmap_coni_rx_fire_sec_growth




