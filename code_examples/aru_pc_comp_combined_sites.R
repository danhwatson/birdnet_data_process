library(unmarked)
library(jagsUI)
library(AICcmodavg)
library(mrds)
library(ggplot2)


########################################################################################################################
###################################################ARU MODELS###########################################################
########################################################################################################################

####################################################WOOD THRUSH#######################################################

detections <- read.csv('examples/WOTH_det_hist_2020-2021_dryad.csv')
woth_aru_models <- data.frame(species=character(),data=character(),model=character(),formulation=character(),
                              site=character(),lambda=numeric(),
                              lambda_lower=numeric(),lambda_upper=numeric(),
                              p=numeric(),p_lower=numeric(),
                              p_upper=numeric(),c.hat=numeric(),convergence=logical())

### ROYLE-NICHOLS ###

#extract just the 0/1 detection history from the detection table
dets <- cbind(detections[,4:13])
dets <- ifelse(dets>0,1,0)

detcovs<-list(wind=detections[,25:34], date=detections[,15:24])

#create "unmarked frame" for use with unmarked occuRN function
umf_rn <- unmarkedFrameOccu(y=dets, obsCovs = detcovs, siteCovs=data.frame("site"=factor(detections$site)))

#fit Royle Nichols model using unmarked
m0 <-occuRN(~ 1 ~ 1, umf_rn, K=150)
m1 <-occuRN(~ 1 ~ -1+site, umf_rn, K=150)
m2 <-occuRN(~ wind ~ -1+site, umf_rn, K=150)
m3 <-occuRN(~ date ~ -1+site, umf_rn, K=150)
m4 <-occuRN(~ wind+date ~ -1+site, umf_rn, K=150)

wothARU_RNtable<-aictab(c(m0, m1, m2, m3, m4), c("Null", "Abu(Site)", "Abu(Site) Det(wind)", "Abu(Site) Det(Date)",
                                                 "Abu(Site) Det(Wind+Date)"))
#view
print(wothARU_RNtable)

#summary of best model
summary(m4)

write.csv(wothARU_RNtable, "woth_aru_RN_aic.csv")

#convergence check
c1<-checkConv(m4)
conv<-c1$converged

#GOF
woth_ARU_RN_GF<-mb.gof.test(m4, nsim=n.sim)
c.hat<-woth_ARU_RN_GF$c.hat.est

#correct for possible overdispersion
m1_summary<-summaryOD(m4, c.hat = c.hat, conf.level = 0.95, 
                      out.type = "confint")

#save fit model results in dataframe

for (i in 1:nsites){
  # lambda, with + and - two standard error
  lambda <- exp(m1_summary$outMat[i,1])
  lambda_p2SE <- exp(m1_summary$outMat[i,1]+1.96*m1_summary$outMat[i,2]) #plus two std err
  lambda_m2SE <- exp(m1_summary$outMat[i,1]-1.96*m1_summary$outMat[i,2]) #minus two std err
  
  # p, with + and - two standard error
  p <- plogis(m1_summary$outMat[nsites+1,1])
  p_p2SE <- plogis(m1_summary$outMat[nsites+1,1]+1.96*m1_summary$outMat[nsites+1,2]) #plus two std err
  p_m2SE <- plogis(m1_summary$outMat[nsites+1,1]-1.96*m1_summary$outMat[nsites+1,2]) #minus two std err
  
  woth_aru_models[nrow(woth_aru_models)+1,]<-list('WOTH','ARU','RN','~Wind+Date ~Site', site[i],
                                                  lambda,lambda_m2SE,lambda_p2SE,p,p_m2SE,p_p2SE,c.hat,conv)
}