# Clear workspace
rm(list=ls())


data = read.csv('data/bacs_val.csv')
data$logit=log(data$confidence/(1-data$confidence))

null.model=glm(correct~1, data, family = 'binomial')
conf.model=glm(correct~confidence, data, family = 'binomial')
logit.model=glm(correct~logit, data, family = 'binomial')
aic.table = AIC(null.model, conf.model, logit.model)
aic.table[order(aic.table$AIC, decreasing=F), ]

# create a df for aic.table, naming columns null_aic, c_aic, and l_aic
aic_table = data.frame(null_aic = aic.table[1,2], c_aic = aic.table[2,2], l_aic = aic.table[3,2])

prediction.range.conf=seq(0,1,.001) 
prediction.range.logit=seq(-3,7,.1) # this is the approximate range of the logit scores

predictions.conf=predict(conf.model, list(confidence=prediction.range.conf), type='r')
predictions.logit=predict(logit.model, list(logit=prediction.range.logit), type='r')

# thresholds for pr(tp)= .90, 0.95, and 0.99
cutoff90_c=(log(.90/(1-.90))-conf.model$coefficients[1])/conf.model$coefficients[2]
cutoff95_c=(log(.95/(1-.95))-conf.model$coefficients[1])/conf.model$coefficients[2]
cutoff99_c=(log(.99/(1-.99))-conf.model$coefficients[1])/conf.model$coefficients[2]

cutoff90_l=(log(.90/(1-.90))-logit.model$coefficients[1])/logit.model$coefficients[2]
cutoff95_l=(log(.95/(1-.95))-logit.model$coefficients[1])/logit.model$coefficients[2]
cutoff99_l=(log(.99/(1-.99))-logit.model$coefficients[1])/logit.model$coefficients[2]

# create a df with columns for each cutoffs
cutoffs = data.frame(cutoff90_c, cutoff95_c, cutoff99_c, cutoff90_l, cutoff95_l, cutoff99_l)

# create a df for sp_code
sp_code = data.frame(sp_code = 'BACS')

# merge the dataframes
thresholds = cbind(sp_code, aic_table, cutoffs)

# read in master .csv
master = read.csv('data/probabilistic_thresholds_master.csv')

# merge the dataframes
thresholds = rbind(master, thresholds)

# write the data to a csv
write.csv(thresholds, 'data/probabilistic_thresholds_master.csv', row.names = FALSE)

#plotting
par(mfrow=c(1,2))

plot(correct~confidence, data, main = 'Confidence scores',
     ylab = 'pr(BirdNET prediction is correct)', xlab = 'confidence score',
     xlim=range(prediction.range.conf), pch=16, cex=1.5, col=rgb(0,0,0,.2))

lines(predictions.conf~prediction.range.conf, lwd=4, col=rgb(0,.75,1,.5))
abline(v=cutoff90_c, col='orange', lwd=4)
abline(v=cutoff95_c, col='red', lwd=4)
abline(v=cutoff99_c, col='magenta', lwd=4)

plot(correct~logit, data, main = 'Logit scores',
     ylab = 'pr(BirdNET prediction is correct)', xlab = 'logit score',
     xlim=range(prediction.range.logit), pch=16, cex=1.5, col=rgb(0,0,0,.2))

lines(predictions.logit~prediction.range.logit, lwd=4, col=rgb(0,.75,1,.5))
abline(v=cutoff90_l, col='orange', lwd=4)
abline(v=cutoff95_l, col='red', lwd=4)
abline(v=cutoff99_l, col='magenta', lwd=4)



