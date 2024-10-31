# clear workspace 
rm(list=ls())

# load packages
library(spOccupancy)

# Load coni abundance data
load("data/abundance_data/abundance_24_coni.RData")

# Convert values greater than 0 or NA to 1 in "y"
abund_data_coni$y <- ifelse(abund_data_coni$y > 0 | is.na(abund_data_coni$y), 1, 0)

# Remove the first 55 visits (columns) from "y" and detection covariates as RN model can't handle NAs
abund_data_coni$y <- abund_data_coni$y[, -c(1:55)]
abund_data_coni$det.covs$day <- abund_data_coni$det.covs$day[, -c(1:55)]
abund_data_coni$det.covs$temp <- abund_data_coni$det.covs$temp[, -c(1:55)]
abund_data_coni$det.covs$wind_speed <- abund_data_coni$det.covs$wind_speed[, -c(1:55)]
abund_data_coni$det.covs$precipitation <- abund_data_coni$det.covs$precipitation[, -c(1:55)]

#rename abund.covs to occ.covs
abund_data_coni$occ.covs <- abund_data_coni$abund.covs

# Treatment occurence formula 
occ.formula.v1 <- ~ factor(treatment) + scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height)

# Without treatment occurence formula 
occ.formula.v2 <- ~ scale(shrub_cover) + scale(grass_cover) + scale(forb_cover) + scale(shrub_height) + scale(grass_height)

# Detection formula 
det.formula <- ~  I(scale(day)^2) + scale(day)+ scale(wind_speed) + scale(precipitation)

# Null model
null_out_coni <- PGOcc(
  occ.formula = ~1,
  det.formula = det.formula,
  data = abund_data_coni,
  inits =,
  n.samples = 10000,
  priors =,
  n.omp.threads = 1,
  verbose = TRUE,
  n.report = 1000,
  n.burn = 1000,
  n.thin = 1,
  n.chains = 5
)
summary(null_out_coni)
save(null_out_coni, file = "data/spOcc_null_out_coni.RData")

# Null model with NNGP
#sp_null_out <- spPGOcc(
#occ.formula = ~1,
#  det.formula = det.formula,
#  data = abund_data_coni,
#  inits =,
# priors =,
#  n.batch = 400,
#  batch.length = 25,
#  NNGP = TRUE,
#  n.neighbors = 5,
#  n.burn = 5000,
#  n.thin = 10,
#  n.chains = 5
#)
#summary(sp_null_out)

# Run the non-spatial model with treatment
t_out_coni <- PGOcc(
  occ.formula = occ.formula.v1,
  det.formula = det.formula,
  data = abund_data_coni,
  inits =,
  n.samples = 10000,
  priors =,
  n.omp.threads = 1,
  verbose = TRUE,
  n.report = 1000,
  n.burn = 1000,
  n.thin = 1,
  n.chains = 5
)
summary(t_out_coni)
#save model 
save(t_out_coni, file = "data/spOcc_t_out_coni.RData")

# Run the non-spatial model without treatment
out_coni <- PGOcc(
  occ.formula = occ.formula.v2,
  det.formula = det.formula,
  data = abund_data_coni,
  inits =,
  n.samples = 10000,
  priors =,
  n.omp.threads = 1,
  verbose = TRUE,
  n.report = 1000,
  n.burn = 1000,
  n.thin = 1,
  n.chains = 5
)
summary(out_coni)
save(out_coni, file = "data/spOcc_out_coni.RData")

plot(t_out_coni, 'beta', density = FALSE) # Good way to visualize convergence 
plot(t_out_coni, 'alpha', density = FALSE)

# Calculate WAIC
waic_null_out_coni <- waicOcc(null_out_coni)
waic_t_out_coni <- waicOcc(t_out_coni)
waic_out_coni <- waicOcc(out_coni) 

print(waic_null_out_coni)
print(waic_t_out_coni)
print(waic_out_coni)

# Posterior predictive checks for the best model (sp_out)
ppc_t_out_coni_group_1 <- ppcOcc(object = t_out_coni, fit.stat = "freeman-tukey", group = 1)

ppc_t_out_coni_group_2 <- ppcOcc(object = t_out_coni, fit.stat = "freeman-tukey", group = 2)

summary(ppc_t_out_coni_group_1)

summary(ppc_t_out_coni_group_2)

diff.fit <- ppc_t_out_coni_group_1$fit.y.rep.group.quants[3, ] - ppc_t_out_coni_group_1$fit.y.group.quants[3, ] 
plot(diff.fit, pch = 19, xlab = 'Sites (group 1 PPC)', ylab = 'Replicate - True Discrepancy')

diff.fit <- ppc_t_out_coni_group_2$fit.y.rep.group.quants[3, ] - ppc_t_out_coni_group_2$fit.y.group.quants[3, ]
plot(diff.fit, pch = 19, xlab = 'Replicate Surveys (group 2 PPC)', ylab = 'Replicate - True Discrepancy')

