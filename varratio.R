## In Gelman's Bayes book page 297.
GelmanRubinRhat = function (chns, columns) {
  out = double(length(columns))
  names(out) = columns
  n = nrow(chns[[1]])
  m = length(chns)
  for (cl in columns) {
    psibar._ = sapply(chns, function (chn) {
      sum(chn[,cl])
    })/n
    psibar.. = sum(sapply(chns, function (chn) mean(chn[,cl])))/m
    B = n / (m-1) * sum((psibar._ - psibar..)^2)
    W = mean(sapply(1:m, function (j) sum(sapply(1:n, function (i) {
      (chns[[j]][i,cl] - psibar._[j])^2 }))/(n-1)))
    varhat = ((n-1)/n)*W+B/n
    Rhat = sqrt(varhat/W)
    out[cl] = Rhat
  }
  out
}
