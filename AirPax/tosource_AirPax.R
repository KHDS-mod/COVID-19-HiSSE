source("getPaxRateProps.R")

## get air passenger intra- and intercontinental proportions
load("AirPax.RData")
lres_pax<-f_getPaxProps(pFeb2020) ## to obtain proportions for other time use the other stored proportion matrices: pJan2020, pDec2019, pNov2019, pYear2019
mInterContinentalPaxProps<-lres_pax$Minter
vIntraContinentalPaxProps<-lres_pax$vintra



## MCMC samples of best found model by BIC
## to check for other models, one needs to change the name of the *.log file
## all *.log files containing the posterior MCMC sample from all the considred models are in this directory
lres_rates<-f_getRateProps("model_QPAIR.log",mInterContinentalPaxProps,burnin=40,b_diag=FALSE,b_cfmean=FALSE)
## Transition rates proportions
EmqProps<-lres_rates$EmProps ## Posterior mean
mLowerQuantile95q<-lres_rates$mLowerQuantile95 ## lower 95% CI
mUpperQuantile95q<-lres_rates$mUpperQuantile95 ## upper 95% CI

## Diversification rates proportions
vlambdas_prop<-lres_rates$vlambdas_prop ## Posterior mean
mLowerQuantile95lambda<-lres_rates$mLowerQuantile95lambda ## lower 95% CI
mUpperQuantile95lambda<-lres_rates$mUpperQuantile95lambda ## upper 95% CI

vlambdas<-lres_rates$vlambdas ## posterior mean diversifictation rates
vdouble_time<-lres_rates$vdouble_time ## posterior mean doubling times

## make scatterplot and calculate correlations
v_q<-c(EmqProps)
v_q<-v_q[-which(is.na(v_q))]
vQL<-c(v_q,vlambdas_prop)
names(vQL)<-NULL

v_inter<-c(mInterContinentalPaxProps)
v_inter<-v_inter[-which(is.na(v_inter))]
vPAX<-c(v_inter,vIntraContinentalPaxProps)
names(vPAX)<-NULL

mQLPAX<-cbind(vQL,vPAX)
sink("Correlation_output.txt")
print("Pearson, all data")
print(cor.test(mQLPAX[,1],mQLPAX[,2],method="pearson"))
print("=======================")
print("Pearson, transition rates")
print(cor.test(mQLPAX[1:15,1],mQLPAX[1:15,2],method="pearson"))
print("=======================")
print("Pearson, diversification rates")
print(cor.test(mQLPAX[16:21,1],mQLPAX[16:21,2],method="pearson"))
print("=======================")
print("Kendall, all rates")
print(cor.test(mQLPAX[,1],mQLPAX[,2],method="kendall"))
print("=======================")
print("Kendall, transition rates")
print(cor.test(mQLPAX[1:15,1],mQLPAX[1:15,2],method="kendall"))
print("=======================")
print("Kendall, diversification rates")
print(cor.test(mQLPAX[16:21,1],mQLPAX[16:21,2],method="kendall"))
sink()

pdf("QLpaxProps.pdf")
plot(mQLPAX[,c(2,1)],pch=19,xlab="passenger proportions",ylab="rate proportions",cex.lab=1.5,cex.axis=1.5,col=c(rep("blue",15),rep("red",6)),cex=1.5)
legend("bottomright",legend=c("intercontinental/transitions","intracontinental/diversifications"),col=c("blue","red"),pch=19,bty="n",cex=1.3)
dev.off()