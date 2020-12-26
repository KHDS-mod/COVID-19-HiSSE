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
})
source("utils.R")
source("constraints.R")

respath = function (file) sprintf("RES_MLESUBPLEX_3/%s", file)
##datpath = function (file) sprintf("ALL_CLEANED_DAT_3/%s", file)
datpath = function (file) sprintf("ALL_CLEANED_DAT_ST6NOULTRA/%s", file)
stord = readRDS(datpath("stord.rds"))
tr = readRDS(datpath("tru.rds"))

## Settings
nstates = length(stord)
root_state = which(stord == "Asia")
constraints = function (likobj) nullmu_ct(likobj, nstates)

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
  dumbinit = starting.point.musse(tr, nstates)##[-(((nstates+1):(2*nstates)))]
  npar = length(dumbinit)
  ## dumbinit[] = 1
  r = subplex(log(dumbinit), #1/starting.point.musse(tr, nstates)[-(((nstates+1):(2*nstates)))],
          function (x) {
            l = likobj(exp(x));
            if (cnt %% 100 == 0) print(c(lik = l, cnt = cnt, x))
            cnt <<- cnt+1
            -l
          },
          control=list(maxit=50000, parscale=rep(1, npar)), hessian=F)
  r$par = exp(r$par)
  r$lnLik = likobj(r$par)
  r
}

tic()
mleest = mle(tr, constraints)
toc()
saveRDS(mleest, respath("mleest.rds"))
