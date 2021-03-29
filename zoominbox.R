suppressPackageStartupMessages({
  library(shape)
})

zoominbox = function (stub, boxmerged) {
  if (stub == "lambsame-withmu") {
    ;;
  } else if (stub == "lambsameqpair-withmu") {
    ;;
  } else if (stub == "qpair-withmu") {
    ## Draw base box plot with grids.
    oldmgp = par('mgp')
    oldtck = par('tck')
    oldmar = par('mar')
    par(mar=oldmar*c(0.7,1.5,0.6,0.7), tck = .01, mgp=oldmgp*c(1,0.3,1))
    boxplot(boxmerged, horizontal=T, las=1, pch=20, outcex=0.6, boxwex=0.4, border = NA,
            axes=F)
    grid(26,26)
    title("Posterior Marginals")
    rect(xleft=-0.07, ybottom=10.5, xright=0.37, ytop=26, lty=2)
    boxplot(boxmerged, horizontal=T, las=1, pch=20, outcex=0.6, boxwex=0.4, outline=F,
            col="white", medlty=0, add=T, cex.axis=0.7)
    df = as.data.frame(boxmerged)
    colmns = colMeans(df)
    points(x=colmns, y=1:ncol(df), pch=20, cex=0.5)
    text(x=colmns, y=(1:ncol(df))+0.42,
         labels=sprintf("%.03f",colmns), cex=0.5)
    Arrows(x0=3/8, y0=17.6, x1 = 0.6, y1 = 19.0, lwd=1.4, arr.type="triangle")
    par(new=T, fig=c(0.45,0.965,0.52, 0.94), tck= .007, mgp=oldmgp*c(1, 0.2, 1), mar = c(1,0.6,0.5,0.3))
    boxplot(boxmerged[11:25], horizontal=T, las=1, pch=20, outcex=0.6, boxwex=0.4, outline=F,
            border = "#000000FF",
            col    = "#FFFFFFFF",
            xaxt   = 'n',
            yaxt   = 'n',
            medlty = 0,
            ylim   = c(0.0, 0.35))
    axis(2, at=1:15, labels=names(boxmerged)[11:25], las=1, cex.axis=0.7)
    axis(1, at=pretty(c(0,0.35)), labels=format(pretty(c(0,0.35),digits=2)), cex.axis=0.7)
    points(x=colmns[11:25], y=1:15, pch=20, cex=0.5)
    par(tck = oldtck, mgp=oldmgp, mar=oldmar)
  } else {
  }
}
