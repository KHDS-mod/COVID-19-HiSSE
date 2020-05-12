suppressPackageStartupMessages({
  library(phytools)
  library(plotrix)
  library(RColorBrewer)
})


mods = c("CORRECT", "LAMBSAME", "QPAIR", "LAMBSAMEQPAIR", "QEQ", "ALLEQ")

for (m in mods) {
  respath = function (x) sprintf("RES_HISSE_ST6NOULTRA_%s_1/%s", m, x)
  character_file = respath("marginal_character.tree")
  astree_file = respath("ast.tree")
  write_pdf = TRUE

  if (write_pdf) pdf(respath("simmap_regionfan.pdf"), width = 50, height = 50)

  sim2 = read.simmap(file=character_file, format="phylip")

  n = 60
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
}
