# Clear Workspace
rm(list = ls())

# Load Packages 
library(tidyverse)
library(spAbundance)
library(ggplot2)

# Load data object
load("data/abundance_data/abundance_24_cwwi.RData")

# Define formulas without treatment
abund.formula <- ~ scale(day) + I(scale(day)^2) + factor(treatment)
det.formula <- ~ scale(day) + I(scale(day)^2) + scale(wind_speed) + scale(precipitation)


out <- NMix(
  abund.formula = abund.formula, 
  det.formula = det.formula, 
  data = abund_data_cwwi, 
  n.batch = 1600,
  batch.length = 25, 
  n.omp.threads = 1,
  n.report = 400,
  family = 'NB',
  verbose = TRUE,
  n.burn = 20000,
  n.thin = 20, 
  n.chains = 3
)

summary(out)
plot(out, 'beta', density = FALSE) # Good way to visualize convergence 
plot(out, 'alpha', density = FALSE)


# Posterior predictive checks for the best model (sp_out)
ppc_sp_out_group_1 <- ppcOcc(object = t_out, fit.stat = "freeman-tukey", group = 1)

ppc_sp_out_group_2 <- ppcOcc(object = sp_out, fit.stat = "freeman-tukey", group = 2)
summary(ppc_sp_out_group_1)
summary(ppc_sp_out_group_2)



diff.fit <- ppc_sp_out_group_1$fit.y.rep.group.quants[3, ] - ppc_sp_out_group_1$fit.y.group.quants[3, ] 
plot(diff.fit, pch = 19, xlab = 'Sites (group 1 PPC)', ylab = 'Replicate - True Discrepancy')

diff.fit <- ppc_sp_out_group_2$fit.y.rep.group.quants[3, ] - ppc_sp_out_group_2$fit.y.group.quants[3, ]
plot(diff.fit, pch = 19, xlab = 'Replicate Surveys (group 2 PPC)', ylab = 'Replicate - True Discrepancy')







