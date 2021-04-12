suppressPackageStartupMessages({
  library("phytools")
  library("diversitree")
})

source("utils.R")

st = read.csv2("RAW_PHYLOGENY/data_MUSSE_26042020.csv")
stord = sort(as.character(unique(st$Region)))

tr = read.tree("RAW_PHYLOGENY/nextstrain_ncov_global_timetree.nwk")
tr$edge.length = tr$edge.length * 365/30

trd = tiprootdist(tr)
tr = multi2di(tr, random=F)
tr$edge.length[which(tr$edge.length <= 0)] = d2e(1/24) * 365/30

## Now make it ultrametric.
tru = collapse.singles(tr)
tru = force.ultrametric(tru, method = "extend")
stopifnot(is.ultrametric(tru))

tru$tip.state = sapply(st$Region, function (x) which(stord == x))
names(tru$tip.state) = tru$tip.label


saveRDS(tru, "ALL_CLEANED_DAT_ULTRAMETRIC/tru.rds", version = 2)
write.nexus(tru, file="ALL_CLEANED_DAT_ULTRAMETRIC/tru.nex")
write.tree(tru, file="ALL_CLEANED_DAT_ULTRAMETRIC/tru.newick")
saveRDS(stord, "ALL_CLEANED_DAT_ULTRAMETRIC/stord.rds", version = 2)
write.csv(data.frame(Tip=tru$tip.label, State=tru$tip.state-1),
          file="ALL_CLEANED_DAT_ULTRAMETRIC/tipstate.csv", row.names=F, quote = F)


## Plot the amount of day extended versus tip number
##pdf("PLOTS/days_extended.pdf")
##hist(e2d((tiprootdist(tru) - trd) / 365 * 30), xlab = "Days",
##     main = "Histogram of number of days extended per tips" )
##dev.off()

##pdf("PLOTS/extended_tree.pdf", width = 210, height = 210)
##plot(tru, show.tip.label=F)
##dev.off()

