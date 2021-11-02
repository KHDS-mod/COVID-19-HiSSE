#!/usr/bin/env Rscript

colours = c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
            "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
            "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
            "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
            "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
            "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
            "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
            "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
            "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
            "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
            "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
            "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
            "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
            "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
            "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
            "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58")
suppressPackageStartupMessages({
  library(functional)
  library(lattice)
  library(latticeExtra)
  library(reshape2)
  library(vioplot)
  library(data.table)
  library(mcmcr)
  library(psych)
})
cmdargs = commandArgs(trailingOnly=TRUE)
chnfiles= cmdargs[-(length(cmdargs)-(0:1))]
plotdir = cmdargs[length(cmdargs)-1]
namestub = cmdargs[length(cmdargs)]
qpair = grepl("qpair", namestub, fixed=TRUE)
chns = lapply(chnfiles, Curry(read.csv, sep="", check.names=F))
maxlen = max(sapply(chns, nrow))
minlen = min(sapply(chns, nrow))
burnin = round(minlen * 0.25)

## Only plot some columns
source("select_columns.R")
allcols = colnames(chns[[1]])
colselect = coldict[[namestub]]

## We trim the chains to the same length because the chains were stopped with a time-budget.
## Except in a non-convergence situation like the unrestricted model, the number of iterations
## contained in each chains are very similar, often varying less than only a dozen of iterations
## or two.
chns_lentrimmed = lapply(chns, function (M) M[1:minlen,])
chns_burnintrimmed = lapply(chns_lentrimmed, function(M) M[-seq_len(burnin),])

## Print the length of each chain...
cat(paste("----", namestub,"----\n"))
cat("The length of the original chain:\n")
print(sapply(chns, nrow))
cat(sprintf("Trimmed it to %d cycles\n", nrow(chns_burnintrimmed[[1]])))

source("varratio.R")
source("plot_rhat.R")
RhatPerColumn = GelmanRubinRhat(chns_burnintrimmed, colselect)
names(RhatPerColumn) = sapply(colselect, function(col)humanise_colnames(colselect,col,qpair=qpair))
cairo_pdf(paste0(plotdir,"/",namestub,"-Rhat.pdf"), family="DejaVu Sans", width=12,height=12)
plot_rhat(RhatPerColumn)
dev.off()

## ------------------------- Trace plot
plts = list()
for (parname in colselect) {
  pl = list()
  if (! parname %in% c("Iteration", "Likelihood", "Prior")) {
    ylim = Reduce(
      function (minmax1, minmax2) {
        c(min(minmax1[1],minmax2[1]), max(minmax1[2],minmax2[2]))
      },
      x = lapply(seq_along(chns),
                 function(j)
                   sapply(list(min,max),
                          function (f) f(chns[[1]][[parname]][-seq_len(burnin)]))),
      init=c(Inf,-Inf))
    panelfn = list( function (...) {
      panel.polygon(c(burnin,burnin,minlen,minlen),
                    c(ylim[1],ylim[2],ylim[2],ylim[1]),
                    col='#cccccc', border="#cccccc")
    } )
    for (j in seq_along(chns)) {
      pl[[j]] = xyplot(chns[[j]][[parname]]~seq_len(nrow(chns[[j]])),
                       xlim=c(1,maxlen),
                       ylim=ylim,
                       type='l',
                       xlab = "Iteration",
                       ylab = "",
                       col=colours[j],
                       panel = function(...) {
                         panel.xyplot(...)
                       })
    }
    basepolygon = xyplot(0~0,
                         xlim=c(1,maxlen),
                         ylim=ylim,
                         xlab = "Iteration",
                         ylab = "",
                         panel=function(...) {
                           panel.polygon(c(burnin,burnin,minlen,minlen),c(-1e6,1e6,1e6,-1e6),
                                         col='#cccccc',
                                         border="#cccccc")
                         })
    plts[[humanise_colnames(colselect,parname,qpair=qpair)]] = Reduce( `+`, pl, basepolygon)
  }
}

plttrace = do.call(c, plts)
cairo_pdf(paste0(plotdir,"/",namestub,"-traceplot.pdf"), width=32, height=32, family="DejaVu Sans")
##trellis.device(device="png", filename=paste0(plotdir,"/",namestub,"-traceplot.png"), width=2560, height=1440)
print(plttrace)
dev.off()


##--------------- Box plot of marginal distribution per chain
i = 1
merged = rbindlist(lapply(chns_burnintrimmed, function (chn) {
  chn[['chainno']] = rep(i, nrow(chn))
  i <<- i + 1
  chn
}))
plts = list()
for (parname in colselect) {
  plts[[humanise_colnames(colselect,parname,qpair=qpair)]] = bwplot(as.formula(sprintf("`%s`~chainno", parname)),
                                                                  horizontal=F,
                                                                  data = merged,
                                                                  pch = '.')
}
plttrace = do.call(c, plts)
cairo_pdf(paste0(plotdir,"/",namestub,"-perchain-box.pdf"), width=32, height=18, family="DejaVu Sans")
print(plttrace)
dev.off()



##--------------- Box plot of all-chain marginal distribution
source('boxmarginal.R')
boxmerged = as.list(merged)[colselect[!(colselect %in% c('Posterior','sig0'))]]
names(boxmerged) = sapply(names(boxmerged), function(col) humanise_colnames(colselect, col, qpair=qpair))
cairo_pdf(paste0(plotdir,"/",namestub,"-marginal-box.pdf"), width=16, height=12, family="DejaVu Sans")
boxmarginal(boxmerged)
dev.off()


##--------------- Histogram of the diversification rates
#xlim = quantile(do.call(c, as.list(merged)[sprintf('lambda_obs[%d]',1:6)]), prob=c(.001,.999))
#cairo_pdf(paste0(plotdir,"/",namestub,"-diversification-rates-hist.pdf"), width=10, height=8)
#par(mar=c(0,0,0,0))
#multi.hist(do.call(cbind, as.list(merged)[sprintf('lambda_obs[%d]',1:6)]),
#           xlim = xlim, main = c("Lambda[Africa]", "Lambda[Asia]", "Lambda[Europe]", "Lambda[N. America]", "Lambda[Oceania]", "Lambda[S. Amer]"))
#dev.off()


##--------------- Compute the BIC
bic = kdict[[namestub]]*log(3585) -2*max(merged[['Posterior']])
sink(paste0(plotdir,"/",namestub,"-bic.txt"))
cat(sprintf("%.03f\n", bic))
sink()
cat(sprintf("BIC=%0.03f\n", bic))

##--------------- Box plot of all-chain marginal distribution but with inset zoom in.
source("zoominbox.R")
cairo_pdf(paste0(plotdir,"/",namestub,"-marginal-box-zoomed.pdf"), width=9, height=8, family="IBM Plex Sans")
zoominbox(namestub, boxmerged)
dev.off()


##--------------- 95% C.I. of various stuff
source('cireport.R')
dfmerged = as.data.frame(merged)
sink(paste0(plotdir,'/',namestub,'-ci.txt'))
ci_report(dfmerged, colselect, qpair=qpair)
sink()
