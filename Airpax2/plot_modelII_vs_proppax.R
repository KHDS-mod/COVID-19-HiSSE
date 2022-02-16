#---- This script plots a model's estimated transition/diversification rates
#     versus the PAX volume. Estimates can be entered here manually or
#     loaded from the MCMC chain into this format; and the PAX
#     data is read automatically from ./PropPaxVol.RData.
#
#     If you write it manually it will be have rounding error, of course!


#------------ START FILL IN THE BLANK
#---- Fill in Lambda's estimates
EST_LAMBDA = c(Afr=0.880, Asia=0.1562 , Eur=1.957 , NAm=1.567, Oc=0.707, SAm=1.233)

#---- Fill in q's estimates here
EST_QPAIR  = c(AfrAsia=0.006,  AfrEur=0.060,  AfrNAm=0.009,  AfrOc=0.035,  AfrSAm=0.056,
                               AsiaEur=0.103, AsiaNAm=0.083, AsiaOc=0.094, AsiaSAm=0.008,
                                              EurNAm=0.267,  EurOc=0.076,  EurSAm=0.061,
                                                             NAmOc=0.086,  NAmSAm=0.021,
                                                                           OcSAm=0.034)
#EST_QFULL  = c(            AfrAsia=,  AfrEur=,  AfrNAm=,  AfrOc=,  AfrSAm=,
#               AsiaAfr=,              AsiaEur=, AsiaNAm=, AsiaOc=, AsiaSAm=,
#               EurAfr=,    EurAsia=,            EurNAm=,  EurOc=,  EurSAm=,
#               NAmAfr=,    NAmAsia=,  NAmEur=,            NAmOc=,  NAmSAm=,
#               OcAfr=,     OcAsia=,   OcEur=,   OcNAm=,         ,  NAmSAm=,
#               SAmAfr=,    SAmAsia=,  SAmEur=,  SAmNAm=,  SAmOc=,          ))
#------------ END FILL IN THE BLANK

load('PropPaxVol.RData')

#---- Get PAX intracontinent
paxintra = l_PropPaxVol$'2020'$APR$vlambda_prop

#---- Get PAX intercontinent
if ('EST_QFULL' %in% ls()) {
  EST_QN = EST_QFULL
  Minter = "Minter"
  paxinter = double(30)
  k = 1L
  for (i in 1:6) {
    for (j in 1:6) {
      if (j != i) {
        paxinter[k] = l_PropPaxVol$'2020'$APR[[Minter]][i,j]
        k = k+1L
      }
    }
  }
} else {
  EST_QN = EST_QPAIR
  Minter = "Minter_sym"
  paxinter = t(l_PropPaxVol$'2020'$APR[[Minter]])[lower.tri(t(l_PropPaxVol$'2020'$APR[[Minter]]))]
}

## Normalize estimates to percentages
EST_LAMBDA = EST_LAMBDA / sum(EST_LAMBDA)
EST_QN     = EST_QN / sum(EST_QN)

c(paxintra, paxinter) -> paxD
c(EST_LAMBDA, EST_QN) -> estD

cat('Total Pearson correlation:\n'); cor.test(paxD, estD, method='pearson')
cat('Total Kendall\'s tau:\n');      cor.test(paxD, estD, method='kendall')
cat('Lamb correlation:\n');          cor.test(paxD[1:6], estD[1:6], method='pearson')
cat('Lamb Kendall\'s tau:\n');       cor.test(paxD[1:6], estD[1:6], method='kendall')
cat('q\'s correlation:\n');          cor.test(paxD[-(1:6)], estD[-(1:6)], method='pearson')
cat('q\'s Kendall\'s tau:\n');       cor.test(paxD[-(1:6)], estD[-(1:6)], method='kendall')

#---- Plot the correlation between PAX and the model
cairo_pdf("QLpaxProps.pdf", width=4.95, height=4.95, family="Libertinus Serif")
plot(x=paxD, y=estD, pch=19,xlab="Passenger Proportions",ylab="Rate Proportions",cex.lab=1.5,cex.axis=1.5,col=c(rep("red",6), rep("blue",15)),cex=1.5)
legend("topright",legend=c("Intercontinental/transitions","Intracontinental/diversifications"),col=c("blue","red"),pch=15,bty="n",cex=0.9)
dev.off()

