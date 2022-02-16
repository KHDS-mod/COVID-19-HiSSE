#!/bin/ksh93 -e
#
#SBATCH -J hisse-lambsameqpair-prior3
#SBATCH -t 48:00:00
#SBATCH -n 32
#SBATCH -o stdout-lambsameqpair-prior3.log
#SBATCH -e stderr-lambsameqpair-prior3.log
#

NCHAIN=32
TASKNAME='lambsameqpair-prior3'
typeset -u OUTDIR
OUTDIR='RES_HISSE_ST6NOULTRA_'${TASKNAME/-/_}_'2'
SRUN_OUTFILE_STUB='stdout-'${TASKNAME}'-%02d.log'
SRUN_ERRFILE_STUB='stderr-'${TASKNAME}'-%02d.log'

mkdir -p "${OUTDIR}"
print "+" `date` >> "${OUTDIR}/timecheckin.txt"
for j in {1..$((${NCHAIN}-1))}; do
    CNO=`printf '%02d' ${j}`
    print \
	'chnname="'${CNO}'"; '   \
	'outdir="'${OUTDIR}'"; ' \
	'source("covid19-genome-hisse-'${TASKNAME}'.rb")'| \
	srun --exclusive -N1 -n1 -o `printf "${SRUN_OUTFILE_STUB}" ${CNO}` -e `printf "${SRUN_ERRFILE_STUB}" ${CNO}` rbnompi &
done
j=${NCHAIN}
CNO=`printf '%02d' ${j}`
print \
    'chnname="'${CNO}'"; '   \
    'outdir="'${OUTDIR}'"; ' \
    'source("covid19-genome-hisse-'${TASKNAME}'.rb")'| \
    srun --exclusive -N1 -n1 -o `printf "${SRUN_OUTFILE_STUB}" ${CNO}` -e `printf "${SRUN_ERRFILE_STUB}" ${CNO}` rbnompi
wait
print "-" `date` >> "${OUTDIR}/timecheckin.txt"
