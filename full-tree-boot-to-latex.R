suppressPackageStartupMessages({
  library("diversitree")
  library("tictoc")
  library("phytools")
  library("R.utils")
  library("reporttools")
})
source("utils.R")
source("constraints.R")

options(width = 180)

respath = function (file) sprintf("RES_MLESUBPLEX_2/%s", file)
datpath = function (file) sprintf("ALL_CLEANED_DAT_2/%s", file)
stord = readRDS(datpath("stord.rds"))
tr = readRDS(datpath("tru.rds"))
mleest = readRDS(respath("mleest.rds"))
mleest$par = abs(mleest$par)

nstates = length(stord)
root_state = which(stord == "Asia")
constraints = function (likobj) nullmu_ct(likobj, nstates)
to_euclid = function (par) nullmu_to_euclid(par, nstates)
from_euclid = function (par) nullmu_from_euclid(par, nstates)


mconfint = function (mlepars, parM, alpha = 0.05) {
  q = t(apply(parM, 1, function (x) {
    quantile(x, c(alpha/2, 1-alpha/2))
  }))
  cbind(2 * mlepars - q[,2], 2 * mlepars - q[,1])
}

mconfint_ln = function (mlepars, parM, alpha = 0.05) {
  q = t(apply(log(parM), 1, function (x) {
    quantile(x, c(alpha/2, 1-alpha/2))
  }))
  exp(cbind(2 * log(mlepars) - q[,2], 2 * log(mlepars) - q[,1]))
}

bcaconfint = function (parM, alpha = 0.05) {
  bca = function (theta, conf.level = 0.95) {
    low <- (1 - conf.level)/2
    high <- 1 - low
    sims <- length(theta)
    z.inv <- length(theta[theta < mean(theta)])/sims
    z <- qnorm(z.inv)
    U <- (sims - 1) * (mean(theta, na.rm = TRUE) - theta)
    top <- sum(U^3)
    under <- 6 * (sum(U^2))^{
      3/2
    }
    a <- top/under
    lower.inv <- pnorm(z + (z + qnorm(low))/(1 - a * (z + qnorm(low))))
    lower <- quantile(theta, lower.inv, names = FALSE)
    upper.inv <- pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
    upper <- quantile(theta, upper.inv, names = FALSE)
    return(c(lower, upper))
  }

  t(apply(parM, 1, function (x) bca(x, 1-alpha) ))
}

qconfint = function (parM, alpha = 0.05) {
  t(apply(parM, 1, function (x) quantile(x, c(alpha/2, 1-alpha/2))))
}

name_ci = function (ci) {
  colnames(ci) = c("ci_low", "ci_up")
  ci
}

res = Filter(Negate(is.null), Map(function (f) {
  tryCatch(readRDS(f),
           error = function (e) {
             print(e)
             NULL
           })
}, Sys.glob(respath("B*.rds"))))

parM = abs(sapply(res, function (x) x$par))
ci_basic = name_ci(mconfint(mleest$par, parM))
ci_ln    = name_ci(mconfint_ln(mleest$par, parM))
ci_bca   = name_ci(bcaconfint(parM))
ci_quant = name_ci(qconfint(parM))

sink(respath("musse_latextab_basic_ci.txt"))
ci_basic_2 = cbind(mleest$par, ci_basic)
colnames(ci_basic_2) = c("MLE", "lower", "upper")
xtable(ci_basic_2, digits=3, cap="MuSSE model: MLE and basic bootstrap C.I. (alpha=0.95)")
sink()
sink(respath("musse_latextab_log_ci.txt"))
ci_ln_2 = cbind(mleest$par, ci_ln)
colnames(ci_ln_2) = c("MLE", "lower", "upper")
xtable(ci_ln_2, digits=3, cap="MuSSE model: MLE and basic bootstrap C.I. on log-scale parametrisation (alpha=0.95)")
sink()
sink(respath("musse_latextab_bca.txt"))
ci_bca_2 = cbind(mleest$par, ci_bca)
colnames(ci_bca_2) = c("MLE", "lower", "upper")
xtable(ci_bca_2, digits=3, cap="MuSSE model: MLE and bootstrap BCa C.I. (alpha=0.95)")
sink()
sink(respath("musse_latextab_quant.txt"))
ci_quant_2 = cbind(mleest$par, ci_quant)
colnames(ci_quant_2) = c("MLE", "lower", "upper")
xtable(ci_quant_2, digits=3, cap="MuSSE model: MLE and bootstrap quantile C.I. (alpha=0.95)")
sink()

sink(respath("ci.txt"))
cat(sprintf("Note: Number of bootstrap samples = %d.\n\n", dim(parM)[2]))
cat("Region numbers\n")
cat("--------------\n\n")
print(data.frame(State_no = 1:length(stord), Regions = stord), row.names=F)
cat(sprintf("\n\nMLE, suplex optimisation, %s iterations\n", mleest$counts))
cat("--------------------------------------------------\n\n")
p = matrix(mleest$par)
rownames(p) = rownames(parM)
print(p)
cat("\n\nBASIC BOOTSTRAP CI, alpha=0.05\n")
cat("------------------------------\n\n")
print(ci_basic)
cat("\n\nBASIC BOOTSTRAP CI REPARAMETRISED ON LOG-SCALE AND CONVERTED BACK, alpha = 0.05\n")
cat("------------------------------\n\n")
print(ci_ln)
cat("\n\nBCa INTERVALS, alpha = 0.05\n")
cat("---------------------------\n")
print(ci_bca)
cat("\n\nQUANTILE BOOSTRAP, alpha = 0.05\n")
cat("-------------------------------\n\n")
print(ci_quant)
sink()
