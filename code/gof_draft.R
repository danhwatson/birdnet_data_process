# Function returning fit-statistics
#fitstats_nobo <- function(rn_model_t_nobo) {
#  observed <- getY(rn_model_t_nobo@data)
#  expected <- fitted(rn_model_t_nobo)
#  resids <- residuals(rn_model_t_nobo)
#  sse <- sum(resids^2, na.rm=TRUE)
#  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
#  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
#  fit_out <- c(cSSE = sse, Chisq = chisq, freemanTukey = freeTuke)
#  return(fit_out)
#}


#parametric bootstrapping 
# Assuming 'mod' is your fitted model
#pb_nobo <- parboot(rn_model_t_nobo, fitstats_nobo, nsim=100)
#print(pb_nobo)
