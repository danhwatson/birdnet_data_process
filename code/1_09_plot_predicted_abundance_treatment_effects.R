# Clear workspace
rm(list=ls())

# Load required libraries
library(ggplot2)
library(tidyverse)
library(viridis)

# Load all the CSV files from the 'means_parameter' folder
#Only Top Models 

mean_bacs_top <- read_csv("data/means_abund_parameters/means_treatment_parameters_bacs.csv")
mean_coni_top <- read_csv("data/means_abund_parameters/means_treatment_parameters_coni.csv")
mean_praw_top <- read_csv("data/means_abund_parameters/means_treatment_parameters_praw.csv")
mean_coye_top <- read_csv("data/means_abund_parameters/means_treatment_parameters_coye.csv")
mean_blgr_top <- read_csv("data/means_abund_parameters/means_treatment_parameters_blgr.csv")
mean_inbu_top <- read_csv("data/means_abund_parameters/means_treatment_parameters_inbu_t.csv")
mean_cwwi_top <- read_csv("data/means_abund_parameters/means_treatment_parameters_cwwi_t.csv")
mean_nobo_top <- read_csv("data/means_abund_parameters/means_treatment_parameters_nobo_t.csv")

#remove trailing edges from treatment column values
mean_bacs_top$treatment <- gsub("_$", "", mean_bacs_top$treatment)
mean_coni_top$treatment <- gsub("_$", "", mean_coni_top$treatment)
mean_praw_top$treatment <- gsub("_$", "", mean_praw_top$treatment)
mean_coye_top$treatment <- gsub("_$", "", mean_coye_top$treatment)
mean_blgr_top$treatment <- gsub("_$", "", mean_blgr_top$treatment)
mean_inbu_top$treatment <- gsub("_$", "", mean_inbu_top$treatment)
mean_cwwi_top$treatment <- gsub("_$", "", mean_cwwi_top$treatment)
mean_nobo_top$treatment <- gsub("_$", "", mean_nobo_top$treatment)


# Add species column to each dataframe
mean_bacs_top$species <- "Bachman's Sparrow"
mean_coni_top$species <- "Common Nighthawk"
mean_praw_top$species <- "Prairie Warbler"
mean_coye_top$species <- "Common Yellowthroat"
mean_blgr_top$species <- "Blue Grosbeak"
mean_inbu_top$species <- "Indigo Bunting"
mean_cwwi_top$species <- "Chuck-will's-widow"
mean_nobo_top$species <- "Northern Bobwhite"


# Combine all dataframes into one
combined_data_top <- bind_rows(
  mean_bacs_top, 
  mean_coni_top, 
  mean_praw_top, 
  mean_coye_top, 
  mean_blgr_top, 
  mean_inbu_top, 
  mean_cwwi_top, 
  mean_nobo_top
)

#Rename the columns
combined_data_top <- combined_data_top %>%
  mutate(treatment = case_match(
    treatment,
    "mine" ~ "Mine Reclamation",
    "mature_savanna" ~ "Mature Savanna",
    "young_savanna" ~ "Young Savanna",
    "timber" ~ "Timber Production",
    .default = treatment 
  ))


# Ensure the treatments are in particular order
combined_data_top$treatment <- factor(combined_data_top$treatment, 
                                    levels = rev(c("Mine Reclamation", "Timber Production", "Young Savanna", "Mature Savanna")))


# Ensure species are in particular order
combined_data_top$species <- factor(combined_data_top$species, 
                                  levels = rev(c("Bachman's Sparrow", "Common Nighthawk", "Prairie Warbler", "Common Yellowthroat", 
                                                 "Blue Grosbeak", "Indigo Bunting", "Chuck-will's-widow", "Northern Bobwhite")))

# Plot
ggplot(combined_data_top, aes(x = predicted_state, y = species, color = treatment)) +
  geom_point(position = position_dodge(width = 0.7), size = 4) +  
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), 
                 height = 0.3, 
                 position = position_dodge(width = 0.7), 
                 size = .4) + 
  scale_color_viridis(option = "viridis", begin = 0.05, end = 0.95, discrete = TRUE, 
                      name = "Site Treatment") +  
  theme_classic() +  
  labs(
    title = "",
    x = "Predicted Relative Abundance",
    y = ""
  ) +
  theme(
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),  
    axis.title.x = element_text(size = 16, face = "bold", vjust=-1.15),
    axis.title.y = element_text(size = 16, face = "bold", vjust=1.15),
    axis.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.key.height = unit(.75, "cm"),
    legend.key.width = unit(1, "cm"),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.position = c(0.85, 0.6),  
    legend.background = element_rect(fill = "white", color = "black", linetype = "solid"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  # Add dashed horizontal lines between species
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5), 
             linetype = "dashed", color = "lightgrey", size = 0.5) +
  guides(color = guide_legend(reverse = TRUE))


