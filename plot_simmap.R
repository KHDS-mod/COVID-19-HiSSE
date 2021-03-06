suppressPackageStartupMessages({
  library(plotrix)
  library(phytools)
  library(RColorBrewer)
  library(RevGadgets)
})

## Change these to plot different models
character_file = "chain_data/RES_HISSE_ST6NOULTRA_QPAIR_WITHMU_2/marginal_character01.tree"
outfile        = "plots/qpair-withmu-simmap01.pdf"

## Whether to plot in new window or plot in PDF
write_pdf      = T

if (write_pdf) pdf(outfile, width = 100, height = 100)

sim2 = read.simmap(file=character_file, format="phylip")

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colv = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
## pie(rep(1,n), col=col_vector)
colsel = colv[c(48, 18, 52, 22, 26, 16,
                 9 , 17, 15, 56, 25, 54)]
colsel = setNames(colsel, as.character(0:11))
plotSimmap(sim2, colsel, fsize=0.7, lwd=1, split.vertical=TRUE, ftype="i")

# add legend
x = 0.05
y = 3000

trans_ord = function (x) {
  dim(x) = c(length(x)/2,2)
  x = t(x)
  c(x)
}

labs = c("Africa [A]", "Asia [A]", "Europe [A]", "N. Amer [A]", "Oceania [A]", "S Amer [A]",
         "Africa [B]", "Asia [B]", "Europe [B]", "N. Amer [B]", "Oceania [B]", "S Amer [B]")
labst = trans_ord(labs)
y = y - (0:(length(colsel) - 1)) * 50
x = rep(x, length(y))
text(x + 0.2, y, labst, pos=4, cex=8)
colselt = trans_ord(colsel)
mapply(draw.circle, x=x, y=y, col=colselt, MoreArgs = list(nv=200, radius=0.02, border="white"))
if (write_pdf)  dev.off()

## Any calls to RevGadgets::plot_ancestral_states fail with "Error in sidx[k]:eidx[k] : NA/NaN argument"
## inside read.beast. No idea why.

