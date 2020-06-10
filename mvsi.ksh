#!/bin/ksh93 -e

srcpath="/home/hckiang/src/wuhan-musse2"
sipath="/home/hckiang/src/ARTICLE/SI"
typeset -A namedict
namedict=( [FULL]=ModelI   [QPAIR]=ModelII [LAMBSAMEQPAIR]=ModelIII \
           [LAMBSAME]=ModelIV [QEQ]=ModelV    [ALLEQ]=ModelVI )
for m in "${!namedict[@]}"; do
        cp "${srcpath}/RES_HISSE_ST6NOULTRA_${m}_1/marginals2.pdf" "${sipath}/marginal_posterior_${namedict[${m}]}.pdf"
	cp "${srcpath}/RES_HISSE_ST6NOULTRA_${m}_1/lambda_plot.pdf" "${sipath}/lambda_posterior_${namedict[${m}]}.pdf"
        mlow="${namedict[${m}]//Model/model}"
        cp "${srcpath}/RES_HISSE_ST6NOULTRA_${m}_1/matplot.pdf" "${sipath}/MCMCchain_${mlow}.pdf"
done
