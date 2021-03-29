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
  library('Rcgmin')
})
source("utils.R")

respath = function (file) sprintf("RES_MLESUBPLEX_WITHOUTMU_1/%s", file)
datpath = function (file) sprintf("ALL_CLEANED_DAT_3/%s", file)
stord = readRDS(datpath("stord.rds"))
tr = readRDS(datpath("tru.rds"))

## Settings
nstates = length(stord)
root_state = which(stord == "Asia")

mle = function (tr, init = NULL) {
  ## We need to set this epsilon, otherwise the DE solver won't terminate, it seems. This is "instability"
  ## also documented in the package help page.
  ##
  likobj = make.musse(tr, tr$tip.state, nstates, sampling.f=NULL, strict=T, control = list(eps=1e-7))
  rootp = numeric(nstates)
  rootp[root_state] = 1
  likobj = set.defaults(likobj, defaults = list(root = quote(ROOT.GIVEN), root.p = rootp))

  cnt = 1
  init = starting.point.musse(tr, nstates)[-(((nstates+1):(2*nstates)))]
  npar = length(init)
  init = runif(npar, min=0, max=1.5)
  r = Rcgmin(log(init),
          function (x) {
            y = exp(x)
            y = c(y[1:nstates], rep(0,nstates),y[(nstates+1):(length(y))])
            if (all(is.finite(y))) {
                tryCatch({
                    l = likobj(y);
                    if (cnt %% 800 == 0) print(c(lik = l, cnt = cnt, y))
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
mleest = mle(tr)
toc()
saveRDS(mleest, respath(sprintf("mleest-%02d.rds",initnum)))
