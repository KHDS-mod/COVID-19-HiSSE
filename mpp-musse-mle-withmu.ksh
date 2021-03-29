#!/bin/ksh93 -e
#
#SBATCH -J musse-withmu
#SBATCH -t 50:00:00
#SBATCH -n 32
#SBATCH -o stdout-musse-withmu-outer.log
#SBATCH -e stderr-musse-withmu-outer.log
#

NINIT=32


SRUN_OUTFILE_STUB='stdout-musse-withmu-%02d.log'
SRUN_ERRFILE_STUB='stderr-musse-withmu-%02d.log'

for j in {1..${NINIT}}; do
    print 'initnum='${j}'; '\
          'set.seed('${j}'); '\
          'source("musse-mle-withmu.R")' |\
      srun -n 1 -o `printf $SRUN_OUTFILE_STUB $j` -e `printf $SRUN_ERRFILE_STUB $j` R --vanilla &
done

wait
