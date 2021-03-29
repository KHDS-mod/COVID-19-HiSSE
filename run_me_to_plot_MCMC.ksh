#!/bin/ksh93 -e

## --- Make file names of every chains
set -A TASKNAME alleq-withmu full-withmu lambsameqpair-withmu lambsame-withmu qeq-withmu qpair-withmu
set -A NCHAIN 32 32 32 32 32 32
for i in "${!TASKNAME[@]}"; do
    OUTDIR[i]='RES_HISSE_ST6NOULTRA_'${TASKNAME[i]/-/_}_'2'
    SRUN_OUTFILE_STUB[i]='stdout-'${TASKNAME}'-%02d.log'
    SRUN_ERRFILE_STUB[i]='stderr-'${TASKNAME}'-%02d.log'
done
typeset -u OUTDIR
PLOTDIR="plots"
LOCALDATAROOT="."

modellogfile() {
    printf "model%02d.log" "$1"
}

# source fetch_from_remote.ksh

plot_allchns() {
    for i in "${!TASKNAME[@]}"; do
	mkdir -p "${LOCALDATAROOT}/${OUTDIR[i]}"
	for j in {1.."${NCHAIN[0]}"}; do
	    FILELIST[j]="${LOCALDATAROOT}/${OUTDIR[i]}/"`modellogfile "${j}"`
	done
	mkdir -p "${PLOTDIR}"
	Rscript  ./plotMC.R "${FILELIST[@]}" "${PLOTDIR}" "${TASKNAME[i]}"
    done
}
best_bic() {
    for i in "${!TASKNAME[@]}"; do
	printf "${TASKNAME[i]}\t"
	cat "${PLOTDIR}"/"${TASKNAME[i]}-bic.txt"
    done
}


#print "Fetching CSV files from Sigma..."
#fetch_csv
print "Plotting MCMC chains"
plot_allchns
best_bic

exit
