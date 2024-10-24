#!/usr/bin/env bash
#SBATCH --job-name=PSMC
#SBATCH -c 1
#SBATCH --time 12:00:00
#SBATCH --mem-per-cpu 1G
#SBATCH -o /home/crq857/projects/Chapter2/slurmout/Shanxi1Wolf.out
#SBATCH -e /home/crq857/projects/Chapter2/slurmout/Shanxi1Wolf.err

#convert bam file to psmcfa

DIR=/projects/mjolnir1/people/crq857/Chapter2/04_Demographichistory/PSMC/fq_files

module load psmc

OUTDIR=/projects/mjolnir1/people/crq857/Chapter2/04_Demographichistory/PSMC/psmcfa

fq2psmcfa -q20 ${DIR}/Shanxi1Wolf.CanFam31.realigned.bam.fq.gz > ${OUTDIR}/Shanxi1Wolf.psmcfa
