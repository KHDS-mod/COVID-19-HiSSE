suppressPackageStartupMessages({
  library("phytools")
  library("diversitree")
})

e2d = function (edgelen) edgelen * 365
d2e = function (days) days / 365
tiprootdist = function (tr) { (dist.nodes(tr))[seq_along(tr$tip.label),tr$edge[1,1]] }
ultrametrise = function (tr) {
  trd = tiprootdist(tr)
  maxdist = max(trd)
  for (i in seq_along(tr$tip.label)) {
    idx = which(tr$edge[,2] == i)
    if(is.na(maxdist - trd[i] + tr$edge.length[idx])) browser()
    tr$edge.length[idx] = maxdist - trd[i] + tr$edge.length[idx]
  }
  tr
}
make_edgelen2day = function (tr, total_days = 113) {
  ratio = total_days/max(tiprootdist(tr))
  function (edgelen) edgelen * ratio
}
