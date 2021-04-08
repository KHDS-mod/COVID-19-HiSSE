suppressPackageStartupMessages({
  library(phytools)
  library(plotrix)
  library(RColorBrewer)
})

## Change these to plot different models
character_file = "chain_data/RES_HISSE_ST6NOULTRA_QPAIR_WITHMU_2/marginal_character01.tree"
outfile        = "plots/qpair-withmu-simmap-regionfan-01.pdf"

## Whether to plot in new window or plot in PDF
write_pdf      = T

if (write_pdf) pdf(outfile, width = 50, height = 50)
sim2 = read.simmap(file=character_file, format="phylip")

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colv = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
## pie(rep(1,n), col=colv)
colsel = colv[rep(c(24, 18, 46, 48, 60, 16),2)]
colsel = setNames(colsel, as.character(0:11))
plotSimmap(sim2, colsel, type = "fan", ftype="off", lwd=3, offset=1, split.vertical = T)
##  plotSimmap(sim2, colsel, type = "fan", fsize=0.7, lwd=3, offset = 1, split.vertical = T)

x = 2.85
y = 4.7

trans_ord = function (x) {
  dim(x) = c(length(x)/2,2)
  x = t(x)
  c(x)
}

labs = c("Africa", "Asia", "Europe", "N. America", "Oceania", "S. America")
y = y - (0:(length(colsel)/2 - 1)) * 0.35
x = rep(x, length(y))
text(x + 0.2, y, labs, pos=4, cex=10.3, font=2)
colselt = colsel[1:6]
mapply(draw.circle, x=x, y=y, col=colselt, MoreArgs = list(nv=200, radius=0.08, border="white"))

if (write_pdf)  dev.off()
