##
## Compatible with pretty-quick-R.
##
## DO NOT READ NEWICK FILES WITH PRETTY-QUICK-R!
## Because both read.tree() in older versions of `ape` returns a wrong tree silently.
## 
options(digits = 5)
suppressPackageStartupMessages({
  library("diversitree")
  library("tictoc")
  library("phytools")
  library("snow")
  library("snowfall")
  library("R.utils")
#  library("dfoptim")
#  library('subplex')
  library('Rcgmin')
})
source("utils.R")
source("constraints.R")

respath = function (file) sprintf("RES_MUSSE_MLE_WITHMU_1/%s", file)
datpath = function (file) sprintf("ALL_CLEANED_DAT_3/%s", file)
##datpath = function (file) sprintf("ALL_CLEANED_DAT_ST6NOULTRA/%s", file)
stord = readRDS(datpath("stord.rds"))
tr = readRDS(datpath("tru.rds"))

## Settings
nstates = length(stord)
root_state = which(stord == "Asia")
constraints = function (likobj) nop_ct(likobj, nstates)


mle = function (tr, add_constraints, init = NULL) {
  ## We need to set this epsilon, otherwise the DE solver won't terminate, it seems. This is "instability"
  ## also documented in the package help page.
  ##
  likobj = make.musse(tr, tr$tip.state, nstates, sampling.f=NULL, strict=T, control = list(eps=1e-7))
  rootp = numeric(nstates)
  rootp[root_state] = 1
  likobj = set.defaults(likobj, defaults = list(root = quote(ROOT.GIVEN), root.p = rootp))
  ## likobj = add_constraints(likobj)
  ## 1. Even if subplex does seem to always reach `maxit`, it usually returns better likelihood
  ##    than optim, which reports convergence, but maybe stuck at some local maximum.
  ## 2. Even if we set reltol to a big value like 1e-5 or 1e-4, it does not make subplex stop in
  ##    a reasonable time. And looking at the Fortran code of the subplex package it doesn't seem
  ##    to take any absolute tolerance threshold.
  ## 3. If you print out the values at each iteration, you'll find that after 5000 iterations
  ##    the likelihood won't change much afterward, so we figured 7500 iterations are
  ##    more than enough. Some q's will move much more (in percentage) than others, after 5000
  ##    iterations, but without influencing the likelihood a lot; this means they have big
  ##    uncertainty. As we do bootstrap anyway there's not much point of getting 0.1% accuracy or
  ##    so in the MLE.
  ##
  cnt = 1
  init = starting.point.musse(tr, nstates)##[-(((nstates+1):(2*nstates)))]
  #init[7:12]=0.05
  init = runif(length(init), min=0, max=1.5)
# init = c(3.2791e+00,  3.0120e-16,  3.3993e-16,  2.7061e-02,
#          3.5042e+00,  4.5646e+00,  3.5777e-20,  2.8748e-20,  1.1612e-22,  7.0491e-17,
#          9.2661e-16,  8.8549e-20,  3.9950e+00,  1.2779e-16,  2.8209e-02,  1.4087e-17,
#          2.6776e-01,  4.7652e-02,  5.1615e-19,  1.4963e-19,  4.6569e-02,  2.7758e-03,
#          6.4063e-02,  1.3203e-19,  2.9178e-20,  2.0372e-18,  1.2965e-01,  2.2318e-02,
#          9.7286e-19,  2.1243e-15,  2.4344e-01,  5.1554e-02,  1.7668e-18,  4.0515e-19,
#          4.1177e+00,  3.8283e-14,  1.2528e-18,  4.8523e-02,  2.7005e-15,  2.2394e-14,
#          5.8985e+00,  1.3186e-01)
  npar = length(init)
#  r = hjk(log(init), function (x) {
#            l = likobj(exp(x));
#            if (cnt %% 100 == 0) print(c(lik = l, cnt = cnt, exp(x)))
#            cnt <<- cnt+1
#            -l
#  })
#  r = subplex(log(init),
#          function (x) {
#            l = likobj(exp(x));
#            if (cnt %% 100 == 0) print(c(lik = l, cnt = cnt, exp(x)))
#            cnt <<- cnt+1
#            -l
#          },
#          control=list(maxit=150000, parscale=rep(1,npar)), hessian=F)
  r = Rcgmin(log(init),
          function (x) {
            y = exp(x)
            if (all(is.finite(y))) {
                tryCatch({
                    l = likobj(exp(x));
                    if (cnt %% 800 == 0) print(c(lik = l, cnt = cnt, exp(x)))
                    cnt <<- cnt+1
                    -l
                }, error= function (e,...) {return(NA);});
            } else {
                return(NA);
            }
          }, control=list(maxit=1500000,trace=T))
 
 
  r$par = exp(r$par)
  r$lnLik = likobj(r$par)
  r
}

tic()
mleest = mle(tr, constraints)
toc()
saveRDS(mleest, respath(sprintf("mleest-%02d.rds",initnum)))
