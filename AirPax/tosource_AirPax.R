source("getPaxRateProps.R")
f_cf_matrices<-function(M1){
    abs(M1-t(M1))
}


## get air passenger intra- and intercontinental proportions
load("AirPax.RData")
lres_pax<-f_getPaxProps(pFeb2020) ## to obtain proportions for other time use the other stored proportion matrices: pJan2020, pDec2019, pNov2019, pYear2019
mInterContinentalPaxProps<-lres_pax$Minter
vIntraContinentalPaxProps<-lres_pax$vintra

## Here we look at proportions for all time periods and the symmetricity (Tab. S.2.)
lAllTimePeriods<-list(y2019=pYear2019,Jan2019=pJan2019,Feb2019=pFeb2019,Nov2019=pNov2019,Dec2019=pDec2019,Jan2020=pJan2020,Feb2020=pFeb2020)
num_TimePeriods<-length(lAllTimePeriods)
lAllTimePeriodsInterIntraContinentalPaxProps<-vector("list",length=num_TimePeriods)
names(lAllTimePeriodsInterIntraContinentalPaxProps)<-names(lAllTimePeriods)

for (i in 1:num_TimePeriods){
    lAllTimePeriodsInterIntraContinentalPaxProps[[i]]<-vector("list",length=7)
    names(lAllTimePeriodsInterIntraContinentalPaxProps[[i]])<-c("mInterContinentalPaxProps","vIntraContinentalPaxProps","mInterContinentalPaxProps_Raw","AvgMdifftM","SdMdifftM","MinMdifftM","MaxMdifftM")
    lres_pax_tmp<-f_getPaxProps(lAllTimePeriods[[i]]) ## to obtain proportions for other time use the other stored proportion matrices: pJan2020, pDec2019, pNov2019, pYear2019
    lAllTimePeriodsInterIntraContinentalPaxProps[[i]]$vIntraContinentalPaxProps<-lres_pax_tmp$Minter
    lAllTimePeriodsInterIntraContinentalPaxProps[[i]]$mInterContinentalPaxProps_raw<-lres_pax_tmp$Minter_raw
    lAllTimePeriodsInterIntraContinentalPaxProps[[i]]$mInterContinentalPaxProps<-lres_pax_tmp$vintra
    mCFmats_cf<-f_cf_matrices(lAllTimePeriodsInterIntraContinentalPaxProps[[i]]$mInterContinentalPaxProps_raw)
    lAllTimePeriodsInterIntraContinentalPaxProps[[i]]$AvgMdifftM<-mean(c(mCFmats_cf))
    lAllTimePeriodsInterIntraContinentalPaxProps[[i]]$SdMdifftM<-sd(c(mCFmats_cf))
    lAllTimePeriodsInterIntraContinentalPaxProps[[i]]$MinMdifftM<-min(c(mCFmats_cf))
    lAllTimePeriodsInterIntraContinentalPaxProps[[i]]$MaxMdifftM<-max(c(mCFmats_cf))
}
## ==================================================================================


## MCMC samples of best found model by BIC
## to check for other models, one needs to change path to the folder which contains all
## the model01.log ... model32.log. This folders are named as RES_HISSE_ST6NOULTRA_*_WITHMU_2
## in the ZIP file.
##
## Link to ZIP file (same as in the top level README):
##    https://liuonline-my.sharepoint.com/:u:/g/personal/haoki85_liu_se/EcPaj_4NIKRKk9LWPBFOrSQBooTUe5bGY1OU7Jdq6_YYkw?e=YNdD7V
##

cfile_name<-Sys.glob(file.path('..','chain_data','RES_HISSE_ST6NOULTRA_QPAIR_WITHMU_2','model*.log'))

if (length(cfile_name) == 0  || ! all(file.exists(cfile_name))){
    stop("Please download and extract the ZIP file containing the MCMC sample for the parameters from https://liuonline-my.sharepoint.com/:u:/g/personal/haoki85_liu_se/EcPaj_4NIKRKk9LWPBFOrSQBooTUe5bGY1OU7Jdq6_YYkw?e=YNdD7V . The folders with the posterior samples are named as RES_HISSE_ST6NOULTRA_*_WITHMU_2. The folder in the script here is the one with the best found by BIC model. If you are interested in another model, then please change RES_HISSE_ST6NOULTRA_QPAIR_WITHMU_2 in the cfile_name variable to what is desired.",call.=FALSE)
}

lres_rates <- f_getRateProps(cfile_name,mInterContinentalPaxProps,burnin=368,b_diag=FALSE,b_cfmean=FALSE)

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

pdf("QLpaxProps.pdf") ## Fig. S.19.
plot(mQLPAX[,c(2,1)],pch=19,xlab="passenger proportions",ylab="rate proportions",cex.lab=1.5,cex.axis=1.5,col=c(rep("blue",15),rep("red",6)),cex=1.5)
legend("bottomright",legend=c("intercontinental/transitions","intracontinental/diversifications"),col=c("blue","red"),pch=19,bty="n",cex=1.3)
dev.off()
