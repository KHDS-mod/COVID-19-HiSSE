suppressPackageStartupMessages({
  library(phytools)
  library(plotrix)
  library(RColorBrewer)
})


## Change these to plot different models
character_file = "chain_data/RES_HISSE_ST6NOULTRA_QPAIR_WITHMU_2/marginal_character01.tree"
outfile        = "plots/qpair-withmu-simmap-fan-01.pdf"

## Whether to plot in new window or plot in PDF
write_pdf      = T

if (write_pdf) pdf(outfile, width = 50, height = 50)

sim2 = read.simmap(file=character_file, format="phylip")
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colv = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
## pie(rep(1,n), col=col_vector)
colsel = colv[c(48, 18, 52, 22, 26, 16,
                9 , 17, 15, 56, 25, 54)]
colsel = setNames(colsel, as.character(0:11))
plotSimmap(sim2, colsel, type = "fan", ftype="off", lwd=3, offset=1, split.vertical = T)

x = 3.52
y = 4.80

trans_ord = function (x) {
  dim(x) = c(length(x)/2,2)
  x = t(x)
  c(x)
}

labs = c("Africa A", "Asia A", "Europe A", "N.Amer. A", "Oceania A", "S.Amer. A",
         "Africa B", "Asia B", "Europe B", "N.Amer. B", "Oceania B", "S.Amer. B")
labst = trans_ord(labs)
y = y - (0:(length(colsel) - 1)) * 0.25
x = rep(x, length(y))
text(x + 0.2, y, labst, pos=4, cex=7.4, font=2)
colselt = trans_ord(colsel)
mapply(draw.circle, x=x, y=y, col=colselt, MoreArgs = list(nv=200, radius=0.08, border="white"))

if (write_pdf)  dev.off()

