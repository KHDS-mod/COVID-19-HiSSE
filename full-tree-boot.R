##
## Compatible with pretty-quick-R.
##
## DO NOT READ NEWICK FILES WITH PRETTY-QUICK-R!
## Because both read.tree() in older versions of `ape` returns a wrong tree silently.
## 

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

respath = function (file) sprintf("RES_MLESUBPLEX_2/%s", file)
datpath = function (file) sprintf("ALL_CLEANED_DAT_2/%s", file)
stord = readRDS(datpath("stord.rds"))
tr = readRDS(datpath("tru.rds"))
mleest = readRDS(respath("mleest.rds"))
mleest$par = abs(mleest$par)

## Settings
nstates = length(stord)
root_state = which(stord == "Asia")
constraints = function (likobj) nullmu_ct(likobj, nstates)
to_euclid = function (par) nullmu_to_euclid(par, nstates)
from_euclid = function (par) nullmu_from_euclid(par, nstates)
cpus = 40
B = 10000

mle = function (tr, add_constraints, to_euclid, init = NULL) {
  ## We need to set this epsilon, otherwise the DE solver won't terminate, it seems. This is "instability"
  ## also documented in the package help page.
  ##
  likobj = make.musse(tr, tr$tip.state, nstates, sampling.f=NULL, strict=T, control = list(eps=1e-7))
  rootp = numeric(nstates)
  rootp[root_state] = 1
  likobj = set.defaults(likobj, defaults = list(root = quote(ROOT.GIVEN), root.p = rootp))
  likobj = add_constraints(likobj)
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
  npar = length(starting.point.musse(tr, nstates)) - nstates
  subplex(to_euclid(starting.point.musse(tr, nstates)),
          function (x) { -likobj(abs(x)) },
          control=list(maxit=30000, parscale=rep(1/4, npar)), hessian=F)
}

parbootest = function (tr, mlepars_full, add_constraints, to_euclid, B=100) {
  ntips = length(tr$tip.label)
  
  est = function (dummy) {
    phy = NULL
    while (is.null(phy)) phy = tree.musse(mlepars_full, ntips, x0=root_state)
    saveRDS((r = mle(phy, add_constraints, to_euclid)),
            respath(paste0("B",strsplit(as.character(runif(1, 0, 1)),"\\.")[[1]][2], ".rds")))
    r
  }

  sfLibrary(ape)
  sfLibrary(phytools)
  sfLibrary(diversitree)
  sfLibrary(deSolve)
  sfLibrary(R.utils)
  sfLibrary(subplex)
  sfExport("mlepars_full", local=T)
  sfExport("ntips", local=T)
  sfExport("add_constraints", local=T)
  
  sfExport("nullmu_from_euclid", local=F)
  sfExport("nullmu_to_euclid", local=F)
  sfExport("nullmu_ct", local=F)
  
  sfExport("nstates", local=F)
  sfExport("to_euclid", local=F)
  sfExport("root_state", local=F)
  sfExport("respath", local=F)
  sfExport("mle", local=F)
  clusterMap(sfGetCluster(), est, rep(T, B))
}
print(root_state)

sfInit(parallel=T, cpus=cpus, type="SOCK")
parbootest(tr, from_euclid(mleest$par), constraints, to_euclid, B)
sfStop()
