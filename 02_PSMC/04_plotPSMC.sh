#!/usr/bin/env bash
#SBATCH --job-name=PSMC
#SBATCH -c 1
#SBATCH --time 03-10:00
#SBATCH --mem-per-cpu 1G
#SBATCH -o /home/crq857/projects/Chapter2/slurmout/PSMC_plot_recheck.out
#SBATCH -e /home/crq857/projects/Chapter2/slurmout/PSMC_plot_recheck.err

module load psmc

## Shanxi 

psmc_plot.pl -u 4.5e-09 -g 4.4 -R -M "sample1=0.12" Shanxi1Wolf1_1_1_1_u4.5e9g4.4_FNR0.12  Shanxi1Wolf1_1_1_1.psmc
