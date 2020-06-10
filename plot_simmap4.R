suppressPackageStartupMessages({
  library(phytools)
  library(plotrix)
  library(RColorBrewer)
})

mods = c("FULL", "LAMBSAME", "QPAIR", "LAMBSAMEQPAIR", "QEQ", "ALLEQ")

for (m in mods) {
  respath = function (x) sprintf("RES_HISSE_ST6NOULTRA_%s_1/%s", m, x)
  character_file = respath("marginal_character.tree")
  astree_file = respath("ast.tree")
  write_pdf = TRUE

  if (write_pdf) pdf(respath("simmap_hiddenfan.pdf"), width = 50, height = 50)

  sim2 = read.simmap(file=character_file, format="phylip")

  n = 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colv = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  ## pie(rep(1,n), col=col_vector)
  colsel = colv[c(rep(18,6),rep(22,6))]
  colsel = setNames(colsel, as.character(0:11))
  plotSimmap(sim2, colsel, type = "fan", ftype="off", lwd=3, offset=1, split.vertical = T)
##  plotSimmap(sim2, colsel, type = "fan", fsize=0.7, lwd=5, offset = 1, split.vertical = T)

  x = 2.45
  y = 4.6

  trans_ord = function (x) {
    dim(x) = c(length(x)/2,2)
    x = t(x)
    c(x)
  }

  labs = c("Hidden = A", "Hidden = B")
  y = y - (0:(length(colsel)/6 - 1)) * .39
  x = rep(x, length(y))
  text(x + 0.2, y, labs, pos=4, cex=12.5, font=2)
  colselt = colsel[c(1,7)]
  mapply(draw.circle, x=x, y=y, col=colselt, MoreArgs = list(nv=200, radius=0.08, border="white"))

  if (write_pdf)  dev.off()
}
