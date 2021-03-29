plot_rhat = function (Rhat, ylim="auto") {
  par(mar=par("mar")*c(1,2.8,1,1))
  upper_Rhat_lim = 1.1
  lower_Rhat_lim = 1.0
  if (ylim == "auto") {
    mx = max(Rhat)
    mn = min(Rhat)
    delta  = mx - mn
    margin = 0.08*delta
    if (mx > 1.18)
      xlim  = c(1-margin, max(mx,upper_Rhat_lim)+margin)
    else
      xlim  = c(0.95,1.18)
  }
  plot(y=seq_along(Rhat), x=Rhat, ylim=c(1,length(Rhat)), xlim=xlim,
       ylab='', pch=20, las=1, yaxt = "n", xlab="Rhat")
  segments(y0=seq(1,length(Rhat)),
           x0=rep(1,length(Rhat)),
           y1=seq(1,length(Rhat)),
           x1=Rhat)
  abline(v=upper_Rhat_lim)
  abline(v=lower_Rhat_lim)
  axis(2, at=seq_along(Rhat), labels=names(Rhat), las=1)
}
