##
## MUST be run by latest CRAN versions of diversitree/ape/phytools
##

suppressPackageStartupMessages({
  library("phytools")
  library("diversitree")
})
source("utils.R")
tr = read.tree("DATA_MUSSE2/nextstrain_ncov_global_timetree.nwk")
tr$edge.length = tr$edge.length * 365/30
trd = tiprootdist(tr)
tr = multi2di(tr, random=F)

tr$edge.length[which(tr$edge.length <= 0)] = d2e(1/24) * 365/30

tru = collapse.singles(tr)
st = read.csv2("DATA_MUSSE2/data_MUSSE_26042020.csv")
stord = sort(as.character(unique(st$Region)))
tru$tip.state = sapply(st$Region, function (x) which(stord == x))
names(tru$tip.state) = tru$tip.label

saveRDS(tru, "ALL_CLEANED_DAT_ST6NOULTRA/tru.rds", version = 2)
saveRDS(stord, "ALL_CLEANED_DAT_ST6NOULTRA/stord.rds", version = 2)
write.nexus(tru, file="ALL_CLEANED_DAT_ST6NOULTRA/tru.nex")
write.tree(tru, file="ALL_CLEANED_DAT_ST6NOULTRA/tru.newick")
write.csv(data.frame(Tip=tru$tip.label, State=tru$tip.state-1),
          file="ALL_CLEANED_DAT_ST6NOULTRA/tipstate.csv", row.names=F, quote = F)

pdf("PLOTS/non-ultra-st6.pdf", width = 210, height = 210)
plot(tru, show.tip.label=F, tip.color = tru$tip.state)
dev.off()


## Plot the amount of day extended versus tip number

pdf("PLOTS/days_extended.pdf")
hist(e2d((tiprootdist(tru) - trd) / 365 * 30), xlab = "Days",
     main = "Histogram of number of days extended per tips" )
dev.off()


