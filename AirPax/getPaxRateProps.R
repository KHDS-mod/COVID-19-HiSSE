f_getPaxProps<-function(M){
    Minter<-M
## [1] "africa"          "asia"            "central america" "europe"         
## [5] "middle east"     "north america"   "oceania"         "south america" 
    Minter["north america",]<-Minter["north america",]+Minter["central america",]
    Minter[,"north america"]<-Minter[,"north america"]+Minter[,"central america"]
    Minter["north america","north america"]<-Minter["north america","north america"]-Minter["central america","central america"]
    
    Minter["asia",]<-Minter["asia",]+Minter["middle east",]
    Minter[,"asia"]<-Minter[,"asia"]+Minter[,"middle east"]
    Minter["asia","asia"]<-Minter["asia","asia"]-Minter["middle east","middle east"]

    Minter<-Minter[-which(rownames(Minter)=="central america"),]
    Minter<-Minter[,-which(colnames(Minter)=="central america")]
    Minter<-Minter[-which(rownames(Minter)=="middle east"),]
    Minter<-Minter[,-which(colnames(Minter)=="middle east")]
        
    vintra<-diag(Minter)
    vintra<-vintra/sum(vintra)
    diag(Minter)<-0
    
    Minter_raw<-Minter
    Minter<-Minter+t(Minter)
    Minter[lower.tri(Minter,diag=TRUE)]<-NA
    Minter<-Minter/sum(Minter,na.rm=TRUE)
    list(Minter=Minter,vintra=vintra,Minter_raw=Minter_raw)
}

f_getRateProps<-function(filenameposterior,mData,burnin=368,b_diag=FALSE,b_cfmean=FALSE){
    filenames <- as.list(filenameposterior)
    if (length(filenames) == 0) stop('MCMC chain data not found. Please download the chain data according the instruction in the top-level README file before running this script.')
    burnin    <- burnin * rep(1, length(filenames))
    dfmodelpost <- do.call(rbind, mapply(function (name, B) {
	read.table(name, header=T)[-(1:B),,drop=F]
    }, name=filenames, B=burnin, SIMPLIFY=F))
    
    ##chain = read.csv("model.log", sep="")
    ##nstates = 6
    dfmodelpost <- dfmodelpost[,c(
		"lambda_hid.1.",  "lambda_hid.2.",
                "lambda_obs.1.",  "lambda_obs.2.",
                "lambda_obs.3.",  "lambda_obs.4.",
                "lambda_obs.5.",  "lambda_obs.6.",
                "q_hid.1.",       "q_hid.2.",
                "q_obs.1.",       "q_obs.2.",       "q_obs.3.",       "q_obs.4.",
                "q_obs.5.",       "q_obs.6.",       "q_obs.7.",       "q_obs.8.",
                "q_obs.9.",       "q_obs.10.",      "q_obs.11.",      "q_obs.12.",
                "q_obs.13.",       "q_obs.14.",      "q_obs.15.",      "q_obs.16.",
                "q_obs.17.",       "q_obs.18.",      "q_obs.19.",      "q_obs.20.",
                "q_obs.21.",       "q_obs.22.",      "q_obs.23.",      "q_obs.24.",
                "q_obs.25.",       "q_obs.26.",      "q_obs.27.",      "q_obs.28.",
                "q_obs.29.",       "q_obs.30.")]
    colnames(dfmodelpost) = c(
		    "lambda_hid.A", "lambda_hid.B.",
                    "lambda_obs.1",  "lambda_obs.2",
                    "lambda_obs.3",  "lambda_obs.4",
                    "lambda_obs.5",  "lambda_obs.6",
                    "q_hid.A",       "q_hid.B",
                   "q_obs.12",      "q_obs.13",       "q_obs.14",      "q_obs.15",      "q_obs.16",
                   "q_obs.21",                       "q_obs.23",       "q_obs.24",      "q_obs.25",      "q_obs.26",
                   "q_obs.31",      "q_obs.32",                        "q_obs.34",      "q_obs.35",      "q_obs.36",
                   "q_obs.41",      "q_obs.42",      "q_obs.43",                        "q_obs.45",      "q_obs.46",
                   "q_obs.51",      "q_obs.52",      "q_obs.53",       "q_obs.54",                       "q_obs.56",
                   "q_obs.61",      "q_obs.62",      "q_obs.63",       "q_obs.64",      "q_obs.65")
    
    lmvProps<-sapply(1:nrow(dfmodelpost),function(k,dfmodelpost,nregions,b_diag,b_cfmean){
	mres<-matrix(NA,nregions,nregions)
	for (i in 1:nregions){
	    for (j in 1:nregions){
		if (i==j){
		    mres[i,i]<-dfmodelpost[k,which(colnames(dfmodelpost)==paste0("lambda_obs.",i))]  
		}else{
		    mres[i,j]<-dfmodelpost[k,which(colnames(dfmodelpost)==paste0("q_obs.",i,j))]  
		}
	    }
	}
	vlambdas<-diag(mres)
	vlambdas_prop<-diag(mres)
	vdouble_time<-1/diag(mres)
	vlambdas_prop<-vlambdas_prop/sum(vlambdas_prop)
	if (!b_diag){
	    diag(mres)<-NA	    
	}	
	if (!b_cfmean){
	    #mres<-mres+t(mres)
	    mres[lower.tri(mres,diag=TRUE)]<-NA
	    mres<-mres/sum(mres,na.rm=TRUE)
	}
	list(mres=mres,vlambdas=vlambdas,vlambdas_prop=vlambdas_prop,vdouble_time=vdouble_time)
    },dfmodelpost=dfmodelpost,nregions=ncol(mData),b_diag=b_diag,b_cfmean=b_cfmean,simplify=FALSE)
    lmProps<-sapply(lmvProps,function(x){x$mres},simplify=FALSE)
    lvlambdaProps<-sapply(lmvProps,function(x){x$vlambdas_prop},simplify=FALSE)
    lvlambda<-sapply(lmvProps,function(x){x$vlambdas},simplify=FALSE)
    lvdouble_time<-sapply(lmvProps,function(x){x$vdouble_time},simplify=FALSE)
    EmProps<-Reduce("+",lmProps)/length(lmProps)
    mLowerQuantile95<-matrix(NA,nrow=nrow(EmProps),ncol=ncol(EmProps))
    mUpperQuantile95<-matrix(NA,nrow=nrow(EmProps),ncol=ncol(EmProps))
    mLowerQuantile95lambda<-rep(NA,length=nrow(EmProps))
    mUpperQuantile95lambda<-rep(NA,length=nrow(EmProps))
    
    for (i in 1:nrow(EmProps)){
	for (j in 1:ncol(EmProps)){
	    v_obs<-sapply(lmProps,function(x){x[i,j]},simplify=TRUE)
	    mLowerQuantile95[i,j]<-quantile(v_obs, probs = 0.025, na.rm = TRUE)
	    mUpperQuantile95[i,j]<-quantile(v_obs, probs = 0.975, na.rm = TRUE)
	}
	v_obs<-sapply(lvlambdaProps,function(x){x[[i]]},simplify=TRUE)
	mLowerQuantile95lambda[i]<-quantile(v_obs, probs = 0.025, na.rm = TRUE)
	mUpperQuantile95lambda[i]<-quantile(v_obs, probs = 0.975, na.rm = TRUE)
    }
    vlambdas_prop<-Reduce("+",lvlambdaProps)/length(lvlambdaProps)
    vlambdas<-Reduce("+",lvlambda)/length(lvlambda)
    vdouble_time<-30*Reduce("+",lvdouble_time)/length(lvdouble_time)
    if(b_cfmean){
	EmProps<-EmProps+t(EmProps)
	EmProps[lower.tri(EmProps)]<-NA
	EmProps<-EmProps/sum(EmProps,na.rm=TRUE)
    }
    colnames(EmProps)<-colnames(mData)
    rownames(EmProps)<-rownames(mData)
    colnames(mLowerQuantile95)<-colnames(mData)
    rownames(mLowerQuantile95)<-rownames(mData)
    colnames(mUpperQuantile95)<-colnames(mData)
    rownames(mUpperQuantile95)<-rownames(mData)
    names(mLowerQuantile95lambda)<-rownames(mData)
    names(mUpperQuantile95lambda)<-colnames(mData)


    names(vlambdas_prop)<-rownames(EmProps)
    names(vlambdas)<-rownames(EmProps)
    names(vdouble_time)<-rownames(EmProps)
    list(EmProps=EmProps,vlambdas=vlambdas,vlambdas_prop=vlambdas_prop,vdouble_time=vdouble_time,mLowerQuantile95=mLowerQuantile95,mUpperQuantile95=mUpperQuantile95,mLowerQuantile95lambda=mLowerQuantile95lambda,mUpperQuantile95lambda=mUpperQuantile95lambda)
}

