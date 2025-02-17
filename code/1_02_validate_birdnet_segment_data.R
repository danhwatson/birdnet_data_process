#Code is based on supplementary material from Wood and Kahl 2024 - Guidelines for appropriate use of BirdNET scores and other detector outputs
#last updated: 2025-02-12

# Clear workspace
rm(list=ls())

library(ggplot2)

# read in data
data = read.csv('data/validation_csv/val_nobo.csv')  

# Remove the '.wav' extension from the 'full_id' column
data$full_id <- gsub(".wav$", "", data$full_id)

# find Duplicates based on the cleaned 'full_id' column
duplicates <- data[duplicated(data$full_id), ]
print(duplicates)

# Remove duplicates based only on the 'full_id' column
# This keeps only the first occurrence of each 'full_id'
data <- data[!duplicated(data$full_id), ]

# create columns for confidence and logit scores
data$confidence <- substr(data$full_id, 1, 5)

# format "confidence" columns as numeric
data$confidence = as.numeric(data$confidence)
data$logit=log(data$confidence/(1-data$confidence))

# format "logit" and "correct" columns as numeric
data$logit = as.numeric(data$logit)
data$correct = as.numeric(data$correct)

#Remove columns titled X
data = data[, !grepl("X", names(data))]


citation("stats")

# run models 
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
sp_code = data.frame(sp_code = 'NOBO')

# merge the dataframes
thresholds = cbind(sp_code, aic_table, cutoffs)

# Count all the correct predictions (1)
thresholds$correct_predicitons = sum(data$correct)

# Count number of rows in data
thresholds$total_predictions = nrow(data)

# read in master .csv
master = read.csv('data/thresholds_master.csv')

# Merge the dataframes
thresholds = rbind(master, thresholds)

# write the data to a csv
write.csv(thresholds, 'data/thresholds_master.csv', row.names = FALSE)

# Set plot output parameters
output_width <- 7   # Width in inches
output_height <- 6   # Height in inches
output_dpi <- 300    # Resolution in DPI

# Set up plotting parameters to save as PNG
png("figures/aos_pres/thresholds_conf_nobo.png", width = output_width, height = output_height, units = "in", res = output_dpi)

# Plot confidence scores
par(mfrow=c(1, 1))
plot(correct ~ confidence, data, main = 'Confidence scores',
     ylab = 'pr(BirdNET prediction is correct)', xlab = 'confidence score',
     xlim = range(prediction.range.conf), pch = 16, cex = 1.5, col = rgb(0, 0, 0, .2))
lines(predictions.conf ~ prediction.range.conf, lwd = 4, col = rgb(0, .75, 1, .5))

# Add the 99% threshold line
abline(v = cutoff99_c, col = 'red', lwd = 4)

# For 90% and 95% thresholds
#abline(v = cutoff90_c, col = 'magenta', lwd = 4)  # 90% threshold
#abline(v = cutoff95_c, col = 'orange', lwd = 4)   # 95% threshold

# Close the plotting device
dev.off()

# Now for the logit scores plot
png("figures/aos_pres/thresholds_logit_nobo.png", width = output_width, height = output_height, units = "in", res = output_dpi)

# Plot logit scores
plot(correct ~ logit, data, main = 'Logit scores',
     ylab = 'pr(BirdNET prediction is correct)', xlab = 'logit score',
     xlim = range(prediction.range.logit), pch = 16, cex = 1.5, col = rgb(0, 0, 0, .2))
lines(predictions.logit ~ prediction.range.logit, lwd = 4, col = rgb(0, .75, 1, .5))
abline(v = cutoff99_l, col = 'red', lwd = 4)

#abline(v = cutoff90_l, col = 'magenta', lwd = 4)  # 90% threshold
#abline(v = cutoff95_l, col = 'orange', lwd = 4)   # 95% threshold

dev.off()

#plotting
#par(mfrow=c(1,2))
#plot(correct~confidence, data, main = 'Confidence scores',
# ylab = 'pr(BirdNET prediction is correct)', xlab = 'confidence score',
# xlim=range(prediction.range.conf), pch=16, cex=1.5, col=rgb(0,0,0,.2))
#lines(predictions.conf~prediction.range.conf, lwd=4, col=rgb(0,.75,1,.5))
#abline(v=cutoff90_c, col='magenta', lwd=4)
#abline(v=cutoff95_c, col='orange', lwd=4)
#abline(v=cutoff99_c, col='red', lwd=4)

#plot(correct~logit, data, main = 'Logit scores',
# ylab = 'pr(BirdNET prediction is correct)', xlab = 'logit score',
# xlim=range(prediction.range.logit), pch=16, cex=1.5, col=rgb(0,0,0,.2))
#lines(predictions.logit~prediction.range.logit, lwd=4, col=rgb(0,.75,1,.5))
#abline(v=cutoff90_l, col='magenta', lwd=4)
#abline(v=cutoff95_l, col='orange', lwd=4)
#abline(v=cutoff99_l, col='red', lwd=4)

