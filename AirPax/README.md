# AIRPAX
AirPax: directory containing data and scripts for comparing estimated rates and intra- and intercontinental airline passenger travel proportions. 

Contents of AirPax
AirPax.RData: data concerning the airline passenger proportions is in the file AirPax.RData and has been obtained from SABRE

.log files: posterior MCMC samples under different model constraints, the best found model by BIC is the file model_QPAIR.log

tosource_AirPax.R: script to run to extract the relevant air travel and rates proportions from the posterior and air travel data

getPaxRateProps.R: functions to calculate posterior mean and CIs of rates, and proportions of intra- and intercontinental air passenger volumes.


# Before running `tosource_AirPax.R`

The scripts in this folder assumes that the numerical result (which contains the MCMC chain data) has already been downloaded and extracted in the correct place. Please see the top-level README for more information about downloading and extracting it.

