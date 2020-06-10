options(width = 160)
paths = sapply(c(FULL="FULL", LAMBSAME="LAMBSAME", QPAIR="QPAIR",
                 LAMBSAMEQPAIR="LAMBSAMEQPAIR", QEQ="QEQ", ALLEQ="ALLEQ"),
               function (x) {
                 function(y) {
                   sprintf("RES_HISSE_ST6NOULTRA_%s_1/%s", x,y)
                 }
               })

cat("BIC:\n")
print(unlist(sapply(paths, function(p) {
  c(unname(read.csv(p("bic.txt"), header=F)))
})))


cat("Log of integrals in the Bayes Factors:\n")
print((bayesfac <- unlist(sapply(paths, function(p) {
  c(unname(read.csv(p("bayesfactor.txt"), header=F)))
}))))

## cat("Bayes Factors table: (log10 based)\n")
## print(outer(bayesfac, bayesfac, function (x,y) { (x-y)/log(10) }))

cat("Number of samples (including burn-in):\n")
print((bayesfac <- unlist(sapply(paths, function(p) {
  tsv = read.csv(p("model.log"), sep="")
  nrow(tsv)
}))))
