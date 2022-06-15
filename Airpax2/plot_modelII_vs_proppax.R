#---- This script plots a model's estimated transition/diversification rates
#     versus the PAX volume. Estimates can be entered here manually or
#     loaded from the MCMC chain into this format; and the PAX
#     data is read automatically from ./PropPaxVol.RData.
#
#     If you write it manually it will be have rounding error, of course!


#------------ START FILL IN THE BLANK
#---- Fill in %Lambda's estimates
EST_LAMBDA = c(Afr=0.111, Asia=0.198, Eur=0.248 , NAm=0.199, Oc=0.089, SAm=0.156)

#---- Fill in %q's estimates here
EST_QPAIR  = c(AfrAsia=0.00571178,  AfrEur=0.06026786,  AfrNAm=0.00935319,  AfrOc=0.03545135,  AfrSAm=0.05530632,
                                    AsiaEur=0.10263116, AsiaNAm=0.08346917, AsiaOc=0.09374223, AsiaSAm=0.00792283,
                                                        EurNAm=0.26731920,  EurOc=0.07647835,  EurSAm=0.06120832,
                                                                            NAmOc=0.08568067,  NAmSAm=0.02133913,
                                                                                               OcSAm=0.03411844)
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


c(paxintra, paxinter) -> paxD
c(EST_LAMBDA, EST_QN) -> estD

cat('#--------- With Asia as outlier ---------\n')
cat('#\tTotal Pearson correlation:\n'); cor.test(paxD, estD, method='pearson')                |> print()
cat('#\tTotal Kendall\'s tau:\n');      cor.test(paxD, estD, method='kendall')                |> print()
cat('#\tLamb correlation:\n');          cor.test(paxD[1:6], estD[1:6], method='pearson')      |> print()
cat('#\tLamb Kendall\'s tau:\n');       cor.test(paxD[1:6], estD[1:6], method='kendall')      |> print()
cat('#\tq\'s correlation:\n');          cor.test(paxD[-(1:6)], estD[-(1:6)], method='pearson')|> print()
cat('#\tq\'s Kendall\'s tau:\n');       cor.test(paxD[-(1:6)], estD[-(1:6)], method='kendall')|> print()

#---- Plot the correlation between PAX and the model
if (! ('EST_QFULL' %in% ls())) {
        cairo_pdf("QLpaxProps.pdf", width=4.5, height=4.5, family="Libertinus Serif")
        par(mar=c(3.1,3.1,0.4,0.2))
        plot(x=paxD, y=estD, pch=19,xlab="Passenger Proportions",ylab="Estimated Rate Proportions",cex.lab=1,cex.axis=1,col=c(rep("red",6), rep("blue",15)),cex=1.2)
        legend("bottomright",legend=c("Intercontinental/transitions","Intracontinental/diversifications"),col=c("blue","red"), pch=15,cex=0.92,bg='#00000000')
        outly_asia = which.max(paxD)
        text(x=paxD[outly_asia]-0.01, y=estD[outly_asia]+0.01, label='Asia')
        dev.off()
}

paxD_noasia = paxD[-outly_asia]
estD_noasia = estD[-outly_asia]

if (! ('EST_QFULL' %in% ls())) {
        cairo_pdf("QLpaxProps_nooutly.pdf", width=4.5, height=4.5, family="Libertinus Serif")
        par(mar=c(3.1,3.1,0.4,0.2))
        plot(x=paxD_noasia, y=estD_noasia, pch=19,xlab="Passenger Proportions",ylab="Estimated Rate Proportions",cex.lab=1,cex.axis=1,col=c(rep("red",5), rep("blue",15)),cex=1.2, asp=1)
        legend("bottomright",legend=c("Intercontinental/transitions","Intracontinental/diversifications"),col=c("blue","red"), pch=15,cex=0.92,bg='#00000000')
        dev.off()
}

## Lambda's correlation without Asia.
{
        cat('#--------- Without Asia as outlier ---------\n')
        cat('#\tTotal Pearson correlation:\n'); cor.test(paxD_noasia, estD_noasia, method='pearson')                |> print()
        cat('#\tTotal Kendall\'s tau:\n');      cor.test(paxD_noasia, estD_noasia, method='kendall')                |> print()
        cat('#\tLamb correlation:\n');          cor.test(paxD_noasia[1:5], estD_noasia[1:5], method='pearson')      |> print()
        cat('#\tLamb Kendall\'s tau:\n');       cor.test(paxD_noasia[1:5], estD_noasia[1:5], method='kendall')      |> print()
        ## These q's correlation will be the same as before removing Asia's lambda
        ##cat('#\tq\'s correlation:\n');          cor.test(paxD_noasia[-(1:5)], estD_noasia[-(1:5)], method='pearson')|> print()
        ##cat('#\tq\'s Kendall\'s tau:\n');       cor.test(paxD_noasia[-(1:5)], estD_noasia[-(1:5)], method='kendall')|> print()
}

months = names(l_PropPaxVol$'2020')[-1]
idx    = which(c(upper.tri(l_PropPaxVol$'2020'$JAN$Minter_sym)))
nm     = outer(rownames(l_PropPaxVol$'2020'$JAN$Minter_sym),
               colnames(l_PropPaxVol$'2020'$JAN$Minter_sym),
               function(x,y) paste0(x,'↔',y))[idx]
ts     = vector('list', length(idx))
names(ts) = nm
for (j in seq_along(idx)) {
        ts[[j]] = numeric(length(months))
        for (k in seq_along(months)) {
                ts[[j]][k] = l_PropPaxVol$'2020'[[months[k]]]$Minter_sym[idx[j]]
        }
}
ts = as.data.frame(ts)
colnames(ts) = nm
rownames(ts) = months
cairo_pdf('intercontinental_passenger_vol_changes.pdf',
          width = 6.0, height=4.2, family="Libertinus Serif")
par(mar=c(3.85,3.85,0.4,7.4))
plot(x=NA, y=NA, ylim=c(-0.06,0.37), xlim=c(0.5,12.5),
     xlab = 'Month in 2020', ylab = 'Passenger volume (%)',
     mgp = c(2.3, 0.2, 0))
for (i in seq_len(ncol(ts))) {
        lines(y=ts[,i], x=1:12, col=i, lty=i, lwd=1.6)
}
legend('right', colnames(ts), col=1:ncol(ts), lty=1:ncol(ts), lwd=1.6,
       inset=-0.30, xpd=T)
dev.off()

library('fishualize'); library('plot.matrix')

cairo_pdf('traffic2020.pdf', width=5.5, height=6.7)
#par(mfrow=c(5,3), mar=c(1.45,1.8,1.0,0.0), family='Libertinus Serif SemiBold')
layout(matrix(c(1, 2, 3,
                     4, 5, 6,
                     7, 8, 9,
                     10, 11, 12,
                     13, 13, 13), # and third plot
            nrow = 5,
            ncol = 3,
            byrow = TRUE),
            heights = c(3,3,3,3, 1))
par(mar=c(1.45,1.6,0.75,0.1), family='Libertinus Serif SemiBold')
for (i in seq_along(months)) {
        XX = l_PropPaxVol$'2020'[[months[i]]]$Minter
        diag(XX) = NA
        plot(XX, border=NA, col=fish(option='Lepomis_megalotis', n=20),
             key=NULL, main=months[i], cex.main=1.0, cex.axis=0.65,
             breaks=seq(from=0.0, to=0.16, length.out=20),
             axis.col=list(padj=-1, cex.axis=0.65))
}
options(digit=2)
#plot(x=NA,y=NA,bty='n',xaxt='n',yaxt='n',xlim=c(0,1),ylim=c(0,1))
par(mar=c(1.5,1.6,0.3,0.1), mgp=c(1,0.0,0.0))
plot(t(matrix(seq(from=0.0, to=0.16, length.out=20), ncol=1L)),
     breaks = seq(from=0.0, to=0.16, length.out=20),
     xlab='', ylab='',
     key=NULL, cex.axis=0.8, border=NA, col=fish(option='Lepomis_megalotis', n=20),
     main='', axis.row=list(labels='',col='#00000000',col.ticks='#00000000',cex=0.7),
     axis.col=list(cex=0.7,padj=-1,labels=sprintf('%.02f',seq(from=0.0, to=0.16, length.out=20))))
dev.off()
