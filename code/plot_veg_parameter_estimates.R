# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)
library(viridis)
library(ggplot2)

# Load data
bacs_mod <- readRDS("data/rn_model_files/rn_model_t_bacs.rds")
blgr_mod <- readRDS("data/rn_model_files/rn_model_t_blgr.rds")
coni_mod <- readRDS("data/rn_model_files/rn_model_t_coni.rds")
coye_mod <- readRDS("data/rn_model_files/rn_model_t_coye.rds")
inbu_mod <- readRDS("data/rn_model_files/rn_model_t_inbu.rds")
praw_mod <- readRDS("data/rn_model_files/rn_model_t_praw.rds")

# Function to extract estimates and SE from model
extract_params <- function(model, species_code) {
  est <- coef(model)
  se <- sqrt(diag(vcov(model)))
  
  data.frame(
    species = species_code,
    parameter = names(est),
    estimate = est,
    std_error = se,
    stringsAsFactors = FALSE
  )
}

# Apply the function to each model
bacs_params <- extract_params(bacs_mod, "BACS")
blgr_params <- extract_params(blgr_mod, "BLGR")
coni_params <- extract_params(coni_mod, "CONI")
coye_params <- extract_params(coye_mod, "COYE")
inbu_params <- extract_params(inbu_mod, "INBU")
praw_params <- extract_params(praw_mod, "PRAW")

# Combine all into a single dataframe
params_df <- bind_rows(
  bacs_params, blgr_params, coni_params, coye_params, inbu_params, praw_params
)

# Filter for the specific parameters you're interested in
filtered_params_df <- params_df %>%
  filter(
    parameter %in% c(
      "lam(scale(shrub_cover))", "lam(scale(grass_cover))", "lam(scale(forb_cover))",
      "lam(scale(shrub_height))", "lam(scale(grass_height))"
    )
  )

# Rename parameters
renamed_params_df <- filtered_params_df %>%
  mutate(
    parameter = case_when(
      parameter == "lam(scale(shrub_cover))" ~ "Shrub Cover",
      parameter == "lam(scale(grass_cover))" ~ "Grass Cover",
      parameter == "lam(scale(forb_cover))" ~ "Forb Cover",
      parameter == "lam(scale(shrub_height))" ~ "Shrub Height",
      parameter == "lam(scale(grass_height))" ~ "Grass Height",
      TRUE ~ parameter
    )
  )

# Reset row names to default numeric sequence
rownames(renamed_params_df) <- NULL

# Calculate 95% confidence intervals
renamed_params_df <- renamed_params_df %>%
  mutate(
    ci_lower = estimate - 1.96 * std_error,
    ci_upper = estimate + 1.96 * std_error
  )

# Specify the order of species (without NOBO, WEVI, CWWI)
species_order <- c(
  "BACS",  # Bachman's Sparrow
  "PRAW",  # Prairie Warbler
  "CONI",   # Common Nighthawk
  "BLGR",  # Blue Grosbeak
  "COYE",  # Common Yellowthroat
  "INBU"  # Indigo Bunting
)

# Reverse the order of species in the factor
renamed_params_df <- renamed_params_df %>%
  mutate(species = factor(species, levels = rev(species_order)))

# Filter data for the cover estimates (shrub cover, grass cover, forb cover)
cover_estimates_data <- renamed_params_df %>%
  filter(parameter %in% c("Shrub Cover", "Grass Cover", "Forb Cover"))

# Filter data for the height estimates (shrub height, grass height)
height_estimates_data <- renamed_params_df %>%
  filter(parameter %in% c("Shrub Height", "Grass Height"))

# Define position dodge to stagger the points and error bars
position_dodge_val <- position_dodge(width = 0.5)

# Create the combined plot with staggered heights for cover estimates
ggplot(cover_estimates_data, aes(x = estimate, y = species, fill = parameter)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, position = position_dodge_val) +
  geom_point(position = position_dodge_val, shape = 21, size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  labs(title = "", x = "Parameter Estimate", y = "Species Alpha Code") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(color = "black", face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 11)  # Make legend values bold
  ) +
  scale_fill_viridis_d(option = "viridis", begin = 0.40, end = 1)

# Optionally save the plot
ggsave("figures/staggered_cover_estimates.png", width = 10, height = 6)

# Create the combined plot with staggered heights for height estimates
ggplot(height_estimates_data, aes(x = estimate, y = species, fill = parameter)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, position = position_dodge_val) +
  geom_point(position = position_dodge_val, shape = 21, size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  labs(title = "", x = "Parameter Estimate", y = "Species Alpha Code") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(color = "black", face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 11)  # Make legend values bold
  ) +
  scale_fill_viridis_d(option = "viridis", begin = 0.40, end = 1)

# Optionally save the plot
ggsave("figures/staggered_height_estimates.png", width = 10, height = 6)
