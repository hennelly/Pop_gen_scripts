#!/usr/bin/env bash
#SBATCH --job-name=psmc
#SBATCH -c 1
#SBATCH --time 2-14:00:00
#SBATCH --mem-per-cpu 5G
#SBATCH -o /home/crq857/projects/Chapter2/slurmout_Oct/downsample.out
#SBATCH -e /home/crq857/projects/Chapter2/slurmout_Oct/downsample.err

module load samtools 

DIR=/projects/mjolnir1/people/crq857/Chapter2/bams_auto
OUT=/projects/mjolnir1/people/crq857/Chapter2/bam_downsample

#Inner Mongolia is 20.04, so 12.35/20.04 is 0.61
samtools view -bs 42.61 ${DIR}/Chinese_CAN6.CanFam31.realigned.bam_autosomes.bam > ${OUT}/Chinese_CAN6.CanFam31_12x.realigned.bam_autosomes.bam
