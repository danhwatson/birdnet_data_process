# Clear workspace
rm(list=ls())

# Load required libraries
library(ggplot2)
library(tidyverse)
library(viridis)

# Load all the CSV files from the 'means_parameter' folder
mean_praw <- read_csv("data/means_abund_parameters/means_treatment_parameters_praw.csv")
mean_nobo <- read_csv("data/means_abund_parameters/means_treatment_parameters_nobo.csv")
mean_inbu <- read_csv("data/means_abund_parameters/means_treatment_parameters_inbu.csv")
mean_cwwi <- read_csv("data/means_abund_parameters/means_treatment_parameters_cwwi.csv")
mean_coye <- read_csv("data/means_abund_parameters/means_treatment_parameters_coye.csv")
mean_coni <- read_csv("data/means_abund_parameters/means_treetment_parameters_coni.csv")
#mean_cgdo <- read_csv("data/means_abund_parameters/means_treatment_parameters_cgdo.csv") #decided to remove due to lack of differences in abundance between treatments
mean_blgr <- read_csv("data/means_abund_parameters/means_treatment_parameters_blgr.csv")
mean_bacs <- read_csv("data/means_abund_parameters/means_treatment_parameters_bacs.csv")
mean_wevi <- read_csv("data/means_abund_parameters/means_treatment_parameters_wevi.csv")


# Add species column to each dataframe
mean_praw$species <- "PRAW"
mean_nobo$species <- "NOBO"
mean_inbu$species <- "INBU"
mean_cwwi$species <- "CWWI"
mean_coye$species <- "COYE"
mean_coni$species <- "CONI"
#mean_cgdo$species <- "CGDO"
mean_blgr$species <- "BLGR"
mean_bacs$species <- "BACS"
mean_wevi$species <- "WEVI"

# Combine all dataframes into one
combined_data <- bind_rows(
  mean_praw, 
  mean_nobo, 
  mean_inbu, 
  mean_cwwi, 
  mean_coye, 
  mean_coni, 
  #mean_cgdo, 
  mean_blgr, 
  mean_bacs, 
  mean_wevi
)

#Rename the columns
combined_data <- combined_data %>%
  mutate(treatment = recode(treatment,
                            "mine" = "Mine",
                            "rx_fire_sec_growth" = "Rx Fire (Sec Growth)",
                            "rx_fire_young" = "Rx Fire (Young)",
                            "timber" = "Timber"
  )) #%>%
#filter(treatment != "Rx Fire (Sec Growth)")  # to exclude this treatment

# Ensure the treatments are in particular order
combined_data$treatment <- factor(combined_data$treatment, 
                                  levels = c("Mine", "Timber", "Rx Fire (Young)", "Rx Fire (Sec Growth)"))


# Ensure species are in particular order 
combined_data$species <- factor(combined_data$species, 
                                levels = rev(c("BACS", "NOBO", "BLGR", "PRAW", 
                                           "COYE", "WEVI", "INBU", "CONI", "CWWI")))


# Plot 
ggplot(combined_data, aes(x = predicted_state, y = species, color = treatment)) +
  geom_point(position = position_dodge(width = 0.7), size = 4) +  
  geom_errorbarh(aes(xmin = lower_CI, 
                     xmax = upper_CI), 
                     height = 0.3, 
                     position = position_dodge(width = 0.7), 
                     size = .4) + 
  scale_color_viridis(option = "cividis", begin = 0.1, end = 0.8, discrete = TRUE, 
                      name = "Treatment") +  
  theme_classic() +  
  labs(
    title = "Predicted Relative Abundance in Response to Site Treatment",
    x = "Relative Abundance",
    y = "Species Code"
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12, face = "bold"),  
    axis.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  # Add dashed horizontal lines between species
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5), 
             linetype = "dashed", color = "lightgrey", size = 0.5)


