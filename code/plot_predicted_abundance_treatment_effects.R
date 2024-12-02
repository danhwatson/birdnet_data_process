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
mean_coni <- read_csv("data/means_abund_parameters/means_treatment_parameters_coni.csv")
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
                                  levels = rev(c("Mine", "Timber", "Rx Fire (Young)", "Rx Fire (Sec Growth)")))

# Ensure species are in particular order 
combined_data$species <- factor(combined_data$species, 
                                levels = rev(c("BACS", "PRAW", "BLGR", "CONI", "COYE", "INBU", "NOBO", "WEVI", "CWWI")))


# Plot  
ggplot(combined_data, aes(x = predicted_state, y = species, color = treatment)) +
  geom_point(position = position_dodge(width = 0.7), size = 4) +  
  geom_errorbarh(aes(xmin = lower_CI, 
                     xmax = upper_CI), 
                     height = 0.3, 
                     position = position_dodge(width = 0.7), 
                     size = .4) + 
  scale_color_viridis(option = "viridis", begin = 0.05, end = 0.95, discrete = TRUE, 
                      name = "Site Treatment") +  
  theme_classic() +  
  labs(
    title = "",
    x = "Predicted Relative Abundance (Site Treatment + Vegetation)",
    y = ""
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
             linetype = "dashed", color = "lightgrey", size = 0.5) +
  guides(color = guide_legend(reverse = TRUE))

ggsave("figures/treatment_veg_predicted_abundance.png", width = 10, height = 6)


# Treatment Only 
mean_praw_t <- read_csv("data/means_abund_parameters/means_treatment_parameters_praw_t.csv")
mean_nobo_t <- read_csv("data/means_abund_parameters/means_treatment_parameters_nobo_t.csv")
mean_inbu_t <- read_csv("data/means_abund_parameters/means_treatment_parameters_inbu_t.csv")
mean_cwwi_t <- read_csv("data/means_abund_parameters/means_treatment_parameters_cwwi_t.csv")
mean_coye_t <- read_csv("data/means_abund_parameters/means_treatment_parameters_coye_t.csv")
mean_coni_t <- read_csv("data/means_abund_parameters/means_treatment_parameters_coni_t.csv")
mean_blgr_t <- read_csv("data/means_abund_parameters/means_treatment_parameters_blgr_t.csv")
mean_bacs_t <- read_csv("data/means_abund_parameters/means_treatment_parameters_bacs_t.csv")
mean_wevi_t <- read_csv("data/means_abund_parameters/means_treatment_parameters_wevi_t.csv")

# Add species column to each dataframe
mean_praw_t$species <- "PRAW"
mean_nobo_t$species <- "NOBO"
mean_inbu_t$species <- "INBU"
mean_cwwi_t$species <- "CWWI"
mean_coye_t$species <- "COYE"
mean_coni_t$species <- "CONI"
mean_blgr_t$species <- "BLGR"
mean_bacs_t$species <- "BACS"
mean_wevi_t$species <- "WEVI"

# Combine all dataframes into one
combined_data_t <- bind_rows(
  mean_praw_t, 
  mean_nobo_t, 
  mean_inbu_t, 
  mean_cwwi_t, 
  mean_coye_t, 
  mean_coni_t, 
  mean_blgr_t, 
  mean_bacs_t, 
  mean_wevi_t
)

#Rename the columns
combined_data_t <- combined_data_t %>%
  mutate(treatment = recode(treatment,
                            "mine" = "Mine",
                            "rx_fire_sec_growth" = "Rx Fire (Sec Growth)",
                            "rx_fire_young" = "Rx Fire (Young)",
                            "timber" = "Timber"
  )) #%>%

# Ensure the treatments are in particular order
combined_data_t$treatment <- factor(combined_data_t$treatment, 
                                  levels = rev(c("Mine", "Timber", "Rx Fire (Young)", "Rx Fire (Sec Growth)")))


# Ensure species are in particular order
combined_data_t$species <- factor(combined_data_t$species, 
                                  levels = rev(c("BACS", "PRAW", "BLGR", "CONI", "COYE", "INBU", "NOBO", "WEVI", "CWWI")))


# Plot  
ggplot(combined_data_t, aes(x = predicted_state, y = species, color = treatment)) +
  geom_point(position = position_dodge(width = 0.7), size = 4) +  
  geom_errorbarh(aes(xmin = lower_CI, 
                     xmax = upper_CI), 
                 height = 0.3, 
                 position = position_dodge(width = 0.7), 
                 size = .4) + 
  scale_color_viridis(option = "viridis", begin = 0.05, end = 0.95, discrete = TRUE, 
                      name = "Site Treatment") +  
  theme_classic() +  
  labs(
    title = "",
    x = "Predicted Relative Abundance (Site Treatment Only)",
    y = ""
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
             linetype = "dashed", color = "lightgrey", size = 0.5) +
  guides(color = guide_legend(reverse = TRUE))

ggsave("figures/treatment_only_predicted_abundance.png", width = 10, height = 6)


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
  mutate(treatment = recode(treatment,
                            "mine" = "Mine Reclamation",
                            "rx_fire_sec_growth" = "Second-Growth Prescribed Fire",
                            "rx_fire_young" = "Young Prescribed Fire",
                            "timber" = "Timber Production"
  )) #%>%

# Ensure the treatments are in particular order
combined_data_top$treatment <- factor(combined_data_top$treatment, 
                                    levels = rev(c("Mine Reclamation", "Timber Production", "Young Prescribed Fire", "Second-Growth Prescribed Fire")))


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
    y = "Species"
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


