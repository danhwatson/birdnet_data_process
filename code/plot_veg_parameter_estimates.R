# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(viridis)
library(ggplot2)

# Load data
bacs_mod <- readRDS("data/rn_model_files/rn_model_t_bacs.rds")
coni_mod <- readRDS("data/rn_model_files/rn_model_t_coni.rds")
praw_mod <- readRDS("data/rn_model_files/rn_model_t_praw.rds")
coye_mod <- readRDS("data/rn_model_files/rn_model_t_coye.rds")
blgr_mod <- readRDS("data/rn_model_files/rn_model_t_blgr.rds")

# Function to extract estimates and SE from model
extract_params <- function(model, species_code) {
  est <- model@estimates@estimates$state@estimates
  se <- sqrt(diag(model@estimates@estimates$state@covMat))
  
  data.frame(
    species = species_code,
    parameter = names(est),
    estimate = est,
    std_error = se,
    stringsAsFactors = FALSE
  )
}

# Apply the function to each model with alpha codes
bacs_params <- extract_params(bacs_mod, "BACS")
coni_params <- extract_params(coni_mod, "CONI")
praw_params <- extract_params(praw_mod, "PRAW")
coye_params <- extract_params(coye_mod, "COYE")
blgr_params <- extract_params(blgr_mod, "BLGR")

# Combine all into a single dataframe
params_df <- bind_rows(bacs_params, coni_params, praw_params, coye_params, blgr_params)

# Rename species to common names after combining
params_df <- params_df %>%
  mutate(
    species = case_when(
      species == "BACS" ~ "Bachman's Sparrow",
      species == "CONI" ~ "Common Nighthawk",
      species == "PRAW" ~ "Prairie Warbler",
      species == "COYE" ~ "Common Yellowthroat",
      species == "BLGR" ~ "Blue Grosbeak",
      TRUE ~ species
    )
  )

# Filter and rename parameters of interest
renamed_params_df <- params_df %>%
  filter(parameter %in% c("scale(shrub_cover)", "scale(grass_cover)", "scale(forb_cover)",
                          "scale(shrub_height)", "scale(grass_height)")) %>%
  mutate(
    parameter = case_when(
      parameter == "scale(shrub_cover)" ~ "Shrub Cover",
      parameter == "scale(grass_cover)" ~ "Grass Cover",
      parameter == "scale(forb_cover)" ~ "Forb Cover",
      parameter == "scale(shrub_height)" ~ "Shrub Height",
      parameter == "scale(grass_height)" ~ "Grass Height"
    ),
    ci_lower = estimate - 1.96 * std_error,
    ci_upper = estimate + 1.96 * std_error
  )

# Specify the order of species
species_order <- c("Bachman's Sparrow", "Common Nighthawk", "Prairie Warbler", "Common Yellowthroat", "Blue Grosbeak")
renamed_params_df <- renamed_params_df %>%
  mutate(species = factor(species, levels = rev(species_order)))

# Define position dodge
position_dodge_val <- position_dodge(width = 0.5)

# Set the desired order of parameter levels
renamed_params_df <- renamed_params_df %>%
  mutate(parameter = factor(parameter, levels = c("Shrub Height", "Grass Height", "Shrub Cover", "Grass Cover", "Forb Cover")))

# Plot with colors in the specified order and matching legend
ggplot(renamed_params_df, aes(x = estimate, y = species, fill = parameter)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3, position = position_dodge_val) +
  geom_point(position = position_dodge_val, shape = 21, size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", size = 0.8) +
  labs(x = "Parameter Estimate", y = "Species", fill = "Parameter") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),  
    axis.title.x = element_text(size = 16, face = "bold", vjust = -1.15),
    axis.title.y = element_text(size = 16, face = "bold", vjust = 1.15),
    axis.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.key.height = unit(0.75, "cm"),
    legend.key.width = unit(1, "cm"),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.position = c(0.85, 0.75),  
    legend.background = element_rect(fill = "white", color = "black", linetype = "solid"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  scale_fill_viridis_d(option = "viridis", begin = 0, end = 1, 
                       labels = c("Shrub Height", "Grass Height", "Shrub Cover", "Grass Cover", "Forb Cover")) +
  
  # Add horizontal lines between species without the top line
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "solid", color = "lightgrey", size = 0.3) +
  
  guides(fill = guide_legend(reverse = TRUE))



# Save the combined plot
ggsave("figures/combined_parameter_estimates.png", width = 10, height = 6)

