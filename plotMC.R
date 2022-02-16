#!/usr/bin/env Rscript

GLOB_PARAM = list(plotdir  = 'plots',
                  chaindir = 'chain_data/RES_HISSE_ST6NOULTRA_{{namestub}}{{run}}/model{{chnno}}.log')
PROC_PARAM = list(
  namestub          = c('full-prior3',        'qpair-prior3',        'lambsameqpair-prior3',     'lambsame-prior3',     'qeq-prior3',          'full-withmu', 'qpair-withmu','lambsameqpair-withmu','lambsame-withmu','qeq-withmu','alleq-withmu'),
  nchains           = c(16,                   16,                    32,                         16,                    32,                    32,            32,             32,                    32,               32,          32           ),
  concats     = list( c('_2'),                c('_2'),               c('_PART1','_2'),           c('_PART1', '_2'),     c('_2'),               c('_2'),       c('_2'),        c('_2'),               c('_2'),          c('_2'),     c('_2')      ),
  burnin            = c(0.5,                  0.5,                   0.5,                        0.6,                   0.5,                   0.25,          0.25,           0.25,                  0.25,             0.25,        0.25         ),
  papername         = c('Model I (Inform.)',  'Model II (Inform.)',  'Model III (Inform.)',     'Model IV (Inform.)',   'Model IV (Inform.)',  'Model I',     'Model II',    'Model III',            'Model IV',       'Model V',   'Model VI'   )
)

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
  library(rmonad)
  library(functional)
  library(magrittr)
  library(lattice)
  library(latticeExtra)
  library(reshape2)
  library(data.table)
  library(mcmcr)
  library(psych)
})
ltheme = canonical.theme(color = FALSE)     ## built-in B&W theme 
ltheme$strip.background$col = 'transparent' ## change strip bg 
lattice.options(default.theme = ltheme)     ## set as default 

lift = function (ms) Map(esc, ms) %>>% {.}
print.mymcmc = function (L, ...)    cat(sprintf('An MCMC sample for \'%s\' consisted of %d independent chains\n', attr(L, 'namestub'), length(L)))
as_mymcmc    = function (L, spec) {
  class(L) = c('mymcmc', class(L))
  attr(L, 'namestub')  = spec[['namestub']]
  attr(L, 'colselect') = coldict[[spec[['namestub']]]]
  cat('--- ALL COLUMN NAMES ---\n');      print(colnames(L[[1]]))
  cat('--- SELECTED COLUMN NAMES ---\n'); print(attr(L, 'colselect'))
  L
}
get_spec     = function (proc_param_) {
  res = list()
  for (i in seq_along(proc_param_[[1]])) {
    res[[i]] = vector("list", length(proc_param_) + 1L)
    names(res[[i]]) = names(proc_param_)
    for (n in names(proc_param_))
      res[[i]][[n]] = proc_param_[[n]][[i]]
  }
  res
}
trim_burnin = function (chns, spec) {
  "We trim the chains to the same length because the chains were stopped with a time-budget.
   Except in a non-convergence situation like the unrestricted model, the number of iterations
   contained in each chains are very similar, often varying less than only a dozen of iterations
   or two."
  maxlen = max(sapply(chns, nrow))
  minlen = min(sapply(chns, nrow))
  burnin = round(minlen * spec[['burnin']])
  chns_lentrimmed = lapply(chns, function (M) M[1:minlen,])
  chns_burnintrimmed = lapply(chns_lentrimmed, function(M) M[-seq_len(burnin),])
  attr(chns_burnintrimmed, 'minlen') = minlen
  attr(chns_burnintrimmed, 'maxlen') = maxlen
  attr(chns_burnintrimmed, 'burnin') = burnin
  attr(chns_burnintrimmed, 'namestub') = spec[['namestub']]
  (attr(chns_burnintrimmed, 'colselect') = attr(chns, 'colselect')) -> colselect
  attr(chns_burnintrimmed, 'posterior_parameter_names') =  colselect[! colselect %in% c("Iteration", "Likelihood", "Prior")]
  class(chns_burnintrimmed) = class(chns)

  ## Print the length of each chain...
  cat(paste("----", spec$namestub,"----\n"))
  cat("The length of the original chain:\n")
  print(sapply(chns, nrow))
  cat(sprintf("Trimmed it to %d cycles\n", nrow(chns_burnintrimmed[[1]])))
  chns_burnintrimmed
}

source("varratio.R")
source("select_columns.R")
source("panel_lollipop.R")

process_model = function (spec) {
  read_chain = function (chnno, spec) {
     mymc = as_monad(spec[['concats']]) %>%
      loop( function (concat, chnno, namestub) {
          model_partial_sub = gsub(x=GLOB_PARAM$chaindir, pattern='\\{\\{namestub\\}\\}', replacement=gsub(x=toupper(namestub), pattern='-', replacement='_'))
          concat_partial_sub = gsub(x=model_partial_sub, pattern='\\{\\{run\\}\\}', replacement=toupper(concat))
          filename = gsub(x=concat_partial_sub, pattern='\\{\\{chnno\\}\\}', replacement=sprintf('%02d', chnno))
          as_monad(read.csv(filename, sep="", check.names=F))
      }, chnno = chnno, namestub = spec[['namestub']]) %*>% rbind
  }
  plot_stuff = function (chns, chns_burnintrimmed, spec) {
    plotdir = GLOB_PARAM$plotdir
    namestub = spec$namestub
    qpair = grepl("qpair", namestub, fixed=TRUE)

    source("plot_rhat.R")
    RhatPerColumn = GelmanRubinRhat(chns_burnintrimmed, attr(chns_burnintrimmed, 'colselect'))
    names(RhatPerColumn) = sapply(attr(chns_burnintrimmed, 'colselect'), function(col)humanise_colnames(attr(chns_burnintrimmed, 'colselect'),col,qpair=qpair))
    cairo_pdf(paste0(plotdir,"/",namestub,"-Rhat.pdf"), family="Liberation Serif", width=12,height=12)
    plot_rhat(RhatPerColumn)
    dev.off()


    ## ------------------------- Trace plot
    plts = list()
    for (parname in attr(chns_burnintrimmed, 'colselect')) {
      pl = list()
      if (parname %in% attr(chns_burnintrimmed, 'posterior_parameter_names')) {
        ylim = Reduce(
          function (minmax1, minmax2) {
            c(min(minmax1[1],minmax2[1]), max(minmax1[2],minmax2[2]))
          },
          x = lapply(seq_along(chns),
                     function(j)
                       sapply(list(min,max),
                              function (f) f(chns[[1]][[parname]][-seq_len(attr(chns_burnintrimmed, 'burnin'))]))),
          init=c(Inf,-Inf))
        panelfn = list( function (...) {
          panel.polygon(c(attr(chns_burnintrimmed, 'burnin'),attr(chns_burnintrimmed, 'burnin'),attr(chns_burnintrimmed, 'minlen'),attr(chns_burnintrimmed, 'minlen')),
                        c(ylim[1],ylim[2],ylim[2],ylim[1]),
                        col='#cccccc', border="#cccccc")
        } )
        for (j in seq_along(chns)) {
          pl[[j]] = xyplot(chns[[j]][[parname]]~seq_len(nrow(chns[[j]])),
                           xlim=c(1,attr(chns_burnintrimmed, 'maxlen')),
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
                             xlim=c(1,attr(chns_burnintrimmed, 'maxlen')),
                             ylim=ylim,
                             xlab = "Iteration",
                             ylab = "",
                             scales=list(x=list(rot=90)),
                             panel=function(...) {
                               panel.polygon(c(attr(chns_burnintrimmed, 'burnin'),attr(chns_burnintrimmed, 'burnin'),attr(chns_burnintrimmed, 'minlen'),attr(chns_burnintrimmed, 'minlen')),c(-1e6,1e6,1e6,-1e6),
                                             col='#cccccc',
                                             border="#cccccc")
                             })
        plts[[humanise_colnames(attr(chns_burnintrimmed, 'colselect'),parname,qpair=qpair)]] = Reduce(`+`, pl, basepolygon)
      }
    }

    plttrace = do.call(c, plts)
    cairo_pdf(paste0(plotdir,"/",namestub,"-traceplot.pdf"), width=20, height=20, family="Liberation Serif")
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
    for (parname in attr(chns_burnintrimmed, 'colselect')) {
      plts[[humanise_colnames(attr(chns_burnintrimmed, 'colselect'),parname,qpair=qpair)]] = bwplot(as.formula(sprintf("`%s`~chainno", parname)),
                                                                      horizontal=F,
                                                                      data = merged,
                                                                      pch = '.')
    }
    plttrace = do.call(c, plts)
    cairo_pdf(paste0(plotdir,"/",namestub,"-perchain-box.pdf"), width=23, height=14, family="Liberation Serif")
    print(plttrace)
    dev.off()



    ##--------------- Box plot of all-chain marginal distribution
    source('boxmarginal.R')
    boxmerged = as.list(merged)[attr(chns_burnintrimmed, 'colselect')[!(attr(chns_burnintrimmed, 'colselect') %in% c('Posterior','sig0'))]]
    names(boxmerged) = sapply(names(boxmerged), function(col) humanise_colnames(attr(chns_burnintrimmed, 'colselect'), col, qpair=qpair))
    cairo_pdf(paste0(plotdir,"/",namestub,"-marginal-box.pdf"), width=9, height=if (length(attr(chns_burnintrimmed, 'colselect')) > 30) 11.5 else 7, family="Liberation Serif")
    boxmarginal(boxmerged)
    dev.off()
    cat('Marginal_posteriors:\n')
    print(apply(as.data.frame(boxmerged), 2, function (x) quantile(x, probs=c(0.05,0.95))))

    ##--------------- Histogram of the diversification rates
    #xlim = quantile(do.call(c, as.list(merged)[sprintf('lambda_obs[%d]',1:6)]), prob=c(.001,.999))
    #cairo_pdf(paste0(plotdir,"/",namestub,"-diversification-rates-hist.pdf"), width=10, height=8)
    #par(mar=c(0,0,0,0))
    #multi.hist(do.call(cbind, as.list(merged)[sprintf('lambda_obs[%d]',1:6)]),
    #           xlim = xlim, main = c("Lambda[Africa]", "Lambda[Asia]", "Lambda[Europe]", "Lambda[N. America]", "Lambda[Oceania]", "Lambda[S. Amer]"))
    #dev.off()


    ##--------------- Compute the BIC
    #bic = kdict[[namestub]]*log(3585) -2*max(merged[['Posterior']])
    #sink(paste0(plotdir,"/",namestub,"-bic.txt"))
    #cat(sprintf("%.03f\n", bic))
    #sink()
    #cat(sprintf("BIC=%0.03f\n", bic))

    ##--------------- Box plot of all-chain marginal distribution but with inset zoom in.
    source("zoominbox.R")
    #cairo_pdf(paste0(plotdir,"/",namestub,"-marginal-box-zoomed.pdf"), width=9, height=8, family="Liberation Serif")
    #cairo_pdf(paste0(plotdir,"/",namestub,"-marginal-box-zoomed.pdf"), width=7, height=7, family="Liberation Serif")
    cairo_pdf(paste0(plotdir,"/",namestub,"-marginal-box-zoomed.pdf"), width=4.95, height=4.95, family="Libertinus Serif Semibold")
    zoominbox(namestub, boxmerged)
    dev.off()


    ##--------------- 95% C.I. of various stuff
    source('cireport.R')
    dfmerged = as.data.frame(merged)
    sink(paste0(plotdir,'/',namestub,'-ci.txt'))
    ci_report(dfmerged, attr(chns_burnintrimmed, 'colselect'), qpair=qpair)
    sink()

    'OK'
  }

  mcchaintag_this        = paste0('MCMCchains:',spec[['namestub']])
  mcchaintagtrimmed_this = paste0('MCMCchains_trimmed:',spec[['namestub']])
  m = spec[['nchains']] %>>% seq_len %>%
    loop(read_chain, spec=spec) %>>%
    as_mymcmc(spec=spec) %>%
    tag('MCMCchains') %>%
    tag(mcchaintag_this) %>>%
    trim_burnin(spec=spec) %>%
    tag('MCMCchains_trimmed') %>%
    tag(mcchaintagtrimmed_this)

  funnel(view(m, mcchaintag_this),
         view(m, mcchaintagtrimmed_this)        ) %*>%
    plot_stuff(spec=spec)
}


compute_allrhat = function (mods, mods_burnintrimmed) {
  'Compute the Rhat statistics for all models and all parameters'
  r = lapply(mods_burnintrimmed, function (chs) {
    rs = GelmanRubinRhat(chs, attr(chs, 'colselect'))
    names(rs) = sapply(attr(chs, 'colselect'),
                       function(col)
                         humanise_colnames(attr(chs, 'colselect'), col, qpair=grepl("qpair", attr(chs, 'namestub'), fixed=TRUE)))
    rs
  })
  names(r) = sapply(mods, function (M) attr(M, 'namestub'))
  r
}

plot_allrhat = function (allrhat) {
  plts = vector("list", length(allrhat)+2L)
  names(plts) = c(sapply(names(allrhat), function(N) PROC_PARAM[['papername']][PROC_PARAM[['namestub']]==N]), 'y.same', 'x.same')
  print(names(plts))
  for (i in seq_along(allrhat)) {
    notprop = which((!grepl('_prop', names(allrhat[[i]]))) & (!grepl('%', names(allrhat[[i]]))))
    namestub = names(allrhat)[i]
    plts[[i]] = xyplot(factor(names(allrhat[[i]])[notprop]) ~ allrhat[[i]][notprop],
           pch = 19, xlab = "Gelman-Rubin's Rhat", ylab = '', lwd = 1,
           xlim = c(0, 1.55),
           scales = list(y = list(rot = 0, cex=0.65)),
           panel = function(x, y, ...) {
             panel.lollipop(x, y, ...)
           },auto.key = list(columns = 6))
  }
  plts[['y.same']] = FALSE
  plts[['x.same']] = FALSE
  plttrace = do.call(c, plts[c(6,5,4,3,2,1,11,10,9,8,7,12,13)])
  cairo_pdf(paste0(GLOB_PARAM$plotdir,"/allrhats.pdf"), width=11, height=16.8, family="Libertinus Serif Semibold")
  print(plttrace)
  dev.off()
  'OK'
}

m1=  get_spec(PROC_PARAM) %>% as_monad %>% loop(process_model)
m = m1 %__%
  funnel(lift(views(m, 'MCMCchains')),
         lift(views(m, 'MCMCchains_trimmed'))) %*>%
  compute_allrhat %>% tag('EveryRhats') %>>% 
  plot_allrhat

# saveRDS(m, file='numerical_result_monad.rds')

## To get the MCMC chain, do
##
##    chns = view(m, 'MCMCchains:full-prior3') %>% esc
##
## or, for the burn-in-trimmed sample, use
##
##    chns_trimmed = view(m, 'MCMCchains_trimmed:full-withmu') %>% esc
##
## You can get all tags and explore the monad by using this command:
##
##    m %>% get_tag %>% unlist %>% unique
##
