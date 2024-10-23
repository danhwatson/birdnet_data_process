# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse)

# Load data
bacs_mod <- readRDS("data/rn_model_files/rn_model_t_bacs.rds")
blgr_mod <- readRDS("data/rn_model_files/rn_model_t_blgr.rds")
coni_mod <- readRDS("data/rn_model_files/rn_model_t_coni.rds")
coye_mod <- readRDS("data/rn_model_files/rn_model_t_coye.rds")
cwwi_mod <- readRDS("data/rn_model_files/rn_model_t_cwwi.rds")
inbu_mod <- readRDS("data/rn_model_files/rn_model_t_inbu.rds")
nobo_mod <- readRDS("data/rn_model_files/rn_model_t_nobo.rds")
praw_mod <- readRDS("data/rn_model_files/rn_model_t_praw.rds")
wevi_mod <- readRDS("data/rn_model_files/rn_model_t_wevi.rds")

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
cwwi_params <- extract_params(cwwi_mod, "CWWI")
inbu_params <- extract_params(inbu_mod, "INBU")
nobo_params <- extract_params(nobo_mod, "NOBO")
praw_params <- extract_params(praw_mod, "PRAW")
wevi_params <- extract_params(wevi_mod, "WEVI")

# Combine all into a single dataframe
params_df <- bind_rows(
  bacs_params, blgr_params, coni_params, coye_params, 
  cwwi_params, inbu_params, nobo_params, praw_params, wevi_params
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

# Assuming renamed_params_df is the dataframe with species (alpha codes), parameter, estimate, and std_error
# Calculate 95% confidence intervals
renamed_params_df <- renamed_params_df %>%
  mutate(
    ci_lower = estimate - 1.96 * std_error,
    ci_upper = estimate + 1.96 * std_error
  )

# Manually specify the order of species using their alpha codes
species_order <- c(
  "BACS",  # Bachman's Sparrow
  "NOBO",  # Northern Bobwhite
  "BLGR",  # Blue Grosbeak
  "PRAW",  # Prairie Warbler
  "COYE",  # Common Yellowthroat
  "WEVI",  # White-eyed Vireo
  "INBU",  # Indigo Bunting
  "CONI",  # Common Nighthawk
  "CWWI"   # Chuck-will's-widow
)

# Reverse the order of species in the factor
renamed_params_df <- renamed_params_df %>%
  mutate(species = factor(species, levels = rev(species_order)))

# Filter data for the cover estimates (shrub cover, grass cover, forb cover)
cover_estimates_data <- renamed_params_df %>%
  filter(parameter %in% c("Shrub Cover", "Grass Cover", "Forb Cover"))

# Define position dodge to stagger the points and error bars
position_dodge_val <- position_dodge(width = 0.5)

# Create the combined plot with staggered heights
ggplot(cover_estimates_data, aes(x = estimate, y = species, color = parameter)) +
  geom_point(position = position_dodge_val) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, position = position_dodge_val) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +  # Vertical line at x = 0
  labs(title = "", x = "Parameter Estimate", y = "Species Alpha Code") +
  theme_minimal() +  # Use default minimal theme
  theme(
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(color = "black", face = "bold")
  ) +
  scale_color_manual(values = c("Forb Cover" = "lightblue", "Grass Cover" = "blue", "Shrub Cover" = "darkblue" ),
                     guide = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_blank())  # Remove legend title

# Optionally save the plot
ggsave("staggered_cover_estimates.png", width = 10, height = 6)


# Filter data for the cover estimates (shrub cover, grass cover, forb cover)
height_estimates_data <- renamed_params_df %>%
  filter(parameter %in% c("Shrub Height", "Grass Height"))

# Define position dodge to stagger the points and error bars
position_dodge_val <- position_dodge(width = 0.5)

# Create the combined plot with staggered heights
ggplot(height_estimates_data, aes(x = estimate, y = species, color = parameter)) +
  geom_point(position = position_dodge_val) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, position = position_dodge_val) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +  # Vertical line at x = 0
  labs(title = "", x = "Parameter Estimate", y = "Species Alpha Code") +
  theme_minimal() +  # Use default minimal theme
  theme(
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(color = "black", face = "bold")
  ) +
  scale_color_manual(values = c("Grass Height" = "blue", "Shrub Height" = "darkblue" ),
                     guide = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_blank())  # Remove legend title

# Optionally save the plot
ggsave("staggered_height_estimates.png", width = 10, height = 6)
