ci_report = function (dfmerged, colselect, qpair=F) {
    orig_colselect = colselect
    colselect = colselect[!(colselect %in% c("Iteration", "Likelihood", "Prior"))]
    
    ## Compute CI for each statistics
    myquant = function (x) c(mean(x), quantile(x, prob=c(0.05,0.95)))
    rawcolci = sapply(dfmerged[,colselect], myquant)
    rownames(rawcolci) = c('posmean', '5%', '95%')
    
    ## Compute lambda proportion CIs
    which_lambdas = colselect[startsWith(colselect, "lambda_obs[")]
    howmany_lambdas = length(which_lambdas)
    if (! (howmany_lambdas == 1)) {
        lambsums = rowSums(dfmerged[,which_lambdas])
        lamb_prop = dfmerged[,which_lambdas]/lambsums
        lamb_prop_ci = sapply(lamb_prop, function (x) {
            r = c(mean(x), quantile(x, prob=c(0.05,0.95)))
            names(r) = NULL
            r
        })
        rownames(lamb_prop_ci) = c('posmean', '5%', '95%')
        colnames(lamb_prop_ci) = sapply(colnames(lamb_prop_ci), function (x)
	    humanise_colnames(colselect, x, qpair = qpair)
	)
    } else {
	lamb_prop_ci = rep(NA,6)
    }
    
    ## Compute q proportion CIs
    which_qs = colselect[startsWith(colselect, "q_obs[")]
    howmany_qs = length(which_qs)
    if (! (howmany_qs == 1)) {
	qsums = rowSums(dfmerged[,which_qs])
	q_prop = dfmerged[,which_qs]/qsums
	q_prop_ci = sapply(q_prop, function (x) {
	    r = c(mean(x), quantile(x, prob=c(0.05,0.95)))
	    names(r) = NULL
	    r
	})
	rownames(q_prop_ci) = c('posmean', '5%', '95%')
	colnames(q_prop_ci) = sapply(colnames(q_prop_ci), function (x)
	    humanise_colnames(colselect, x, qpair = qpair)
	)
    } else {
	q_prop_ci = NA
    }
    cat('------- Parameter C.I. --------\n')
    print(rawcolci)
    
    cat('------- Diversification Proportion C.I. ---------\n')
    print(lamb_prop_ci)
    
    cat('------- Transition Proportion C.I. ---------\n')
    print(q_prop_ci)
}