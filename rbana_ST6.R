library("matrixStats")
library("lattice")
library("reshape")
library("psych")

## Number of parameters:
##          lambda_hid    lambda_obs      q_hid      q_obs        total
## CORRECT:          1             6          1         30           38
## LAMBSAME:         1             1          1         30           33
## QPAIR:            1             6          1         15           23
## LAMBSAMEQPAIR:    1             1          1         15           18
## QEQ:              1             6          1          1            9
## ALLEQ:            1             1          1          1            4

ks      = c(       38,         33,      23,              18,     9,       4)
resdirs = c("CORRECT", "LAMBSAME", "QPAIR", "LAMBSAMEQPAIR", "QEQ", "ALLEQ")

for (m in seq_along(ks)) {
  k = ks[m]; 
  respath = function (x) sprintf("RES_HISSE_ST6NOULTRA_%s_1/%s", resdirs[m], x)

  chain = read.csv(respath("model.log"), sep="")
  nstates = 6
  parM = chain[,c("lambda_hid.1.",  "lambda_hid.2.",
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
  colnames(parM) = c("lambda_hid.A", "lambda_hid.B.",
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


  nburnin = 40
  pdf(respath("splom.pdf"), width = 65, height = 65)
  print(splom(parM[nburnin:nrow(parM),], pch = 20, main = sprintf("Posterior from %d samples, from idx = %d to the end.", nrow(parM)-nburnin, nburnin)))
  dev.off()

  pdf(respath("matplot.pdf"), width=11, height=11)
  layout(matrix(c(1,2),nrow=1), width=c(2.7,1))
  par(mar=c(5,4,4,0)) #No margin on the right side
  matplot(parM, type="l", col = 1:ncol(parM), lty = 1:ncol(parM), main="MCMC Convergence")
  par(mar=c(5,0,4,2))
  plot(c(0,1),type="n", axes=F, xlab="", ylab="")
  legend(1.0, 0.99, legend = colnames(parM), col = 1:ncol(parM), lty = 1:ncol(parM))
  dev.off()


  pdf(respath("marginals.pdf"), width=6, height = 10)
  print(bwplot(variable ~ value, data=melt(parM[nburnin:nrow(parM),]),
         panel = function(x, y, ...){
           panel.bwplot(x, y, ...)
           panel.text(x=colMeans(parM[nburnin:nrow(parM),]),
                      y=1:ncol(parM)+0.44, labels=sprintf("%.03f",colMeans(parM[nburnin:nrow(parM),])), cex=0.5)
         },
         main = "Posterior Marginals", pch = 18))
  dev.off()

  ## Marginal plot v2
  
  cairo_pdf(respath("marginals2.pdf"), width=6, height = 8.8, family="DejaVu Sans")
  parM2 = parM[nburnin:nrow(parM),c("q_obs.12","q_obs.13","q_obs.14","q_obs.15","q_obs.16",
                                    "q_obs.23","q_obs.24","q_obs.25","q_obs.26",
                                    "q_obs.34","q_obs.35","q_obs.36",
                                    "q_obs.45","q_obs.46",
                                    "q_obs.56",
                                    "q_hid.A", "lambda_hid.A", "lambda_hid.B.",
                                    "lambda_obs.1",  "lambda_obs.2",
                                    "lambda_obs.3",  "lambda_obs.4",
                                    "lambda_obs.5",  "lambda_obs.6")]
  names(parM2) = c("q: Africa<->Asia","q: Africa<->Europe","q: Africa<->N.Amer","q: Africa<->Oceania","q: Africa<->S.Amer",
                   "q: Asia<->Europe","q: Asia<->N.Amer","q: Asia<->Oceania","q: Asia<->S.Amer",
                   "q: Europe<->N.Amer","q: Europe<->Oceania","q: Europe<->S.Amer",
                   "q: N.Amer<->Oceania","q: N.Amer<->S.Amer",
                   "q: Oceania<->S.Amer",
                   "q: hidden state", "λ: hidden state A", "λ: hidden state B",
                   "λ: Africa","λ: Asia", "λ: Europe","λ: N.Amer", "λ: Oceania","λ: S.Amer")
  mp2dat = melt((parMmp2 = parM2))
  print(bwplot(variable ~ value, data=mp2dat,
         panel = function(x, y, ...){
           panel.bwplot(x, y, ...)
           panel.text(x=colMeans(parMmp2),
                      y=1:ncol(parMmp2)+0.44, labels=sprintf("%.03f",colMeans(parMmp2)), cex=0.5)
           panel.lines(x=c(-1,3), y=c(16.65,16.65), col="black", lty=3)
         },
         par.settings = list(box.rectangle = list(col="black"),
                             box.umbrella = list(col="black"),
                             plot.symbol = list(col="black", pch = 20)),
         main = "Posterior Marginals", pch = 20, col=1, do.out=F, cex=0.5, xlab = ""))
  dev.off()

  ## BIC score
  mapidx = which.max(chain$Posterior)
  bic = k*log(3585) -2*chain$Posterior[mapidx]
  sink(respath("bic.txt"))
  cat(sprintf("%.03f\n", bic))
  sink()

  sink(respath("bayesfactor.txt"))
  bayesfac = log(nrow(chain) - nburnin) - logSumExp(-chain$Likelihood[nburnin:nrow(chain)])
  cat(sprintf("%.03f\n", bayesfac))
  sink()

  ## Plot Lambda's
  pdf(respath("lambda_plot.pdf"))
  xlim = quantile((tmp <- c(as.matrix(parM[,c("lambda_obs.1",  "lambda_obs.2",
                                              "lambda_obs.3",  "lambda_obs.4",
                                              "lambda_obs.5",  "lambda_obs.6")]))), c(0.001, .999))

  multi.hist(as.matrix(parM[,c("lambda_obs.1",  "lambda_obs.2",
                               "lambda_obs.3",  "lambda_obs.4",
                               "lambda_obs.5",  "lambda_obs.6")]),
             xlim = xlim, main = c("Lambda[Africa]", "Lambda[Asia]", "Lambda[Europe]", "Lambda[N. America]", "Lambda[Oceania]", "Lambda[S. Amer]"))
  dev.off()
}
