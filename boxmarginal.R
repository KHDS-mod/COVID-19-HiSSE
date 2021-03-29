boxmarginal = function (boxmerged, ...) {
  par(mar=par("mar")*c(1,2.8,1,1))
  boxplot(boxmerged, horizontal=T, las=1, pch=20, outcex=0.6, boxwex=0.4, outline=F,
          col="white", medlty=0, ...)
  df = as.data.frame(boxmerged)
  colmns = colMeans(df)
  points(x=colmns, y=1:ncol(df), pch=20, cex=0.5)
  text(x=colmns, y=(1:ncol(df))+0.36,
       labels=sprintf("%.03f",colmns), cex=0.6)
}
