#!/usr/bin/env bash
#SBATCH --job-name=PSMC
#SBATCH -c 1
#SBATCH --time 12:00:00
#SBATCH --mem-per-cpu 5G
#SBATCH -o /home/crq857/projects/Chapter2/slurmout/PSMC_psmc_rerun_testl.out
#SBATCH -e /home/crq857/projects/Chapter2/slurmout/PSMC_psmc_rerun_testlc.err

DIR=/projects/mjolnir1/people/crq857/Chapter2/04_Demographichistory/PSMC/psmcfa
LOW=Chinese_CAN6.CanFam31_12x.psmcfa
HIGH=Chinese_CAN6.CanFam31.realigned.bam_autosomes.bam.psmcfa
OUTDIR=/projects/mjolnir1/people/crq857/Chapter2/04_Demographichistory/PSMC/plot_Shanxi

module load psmc

psmc -p "1+1+1+1+25*2+4+6" -r5 -t15 -N25 -o ${OUTDIR}/Chinese_CAN6.CanFam31_12x_1_1_1_1.psmc ${DIR}/Chinese_CAN6.CanFam31_12x.psmcfa
