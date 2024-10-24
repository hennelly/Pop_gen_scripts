#!/usr/bin/env bash
#SBATCH --job-name=keepauto
#SBATCH -c 1
#SBATCH --time 1-12:00:00
#SBATCH --mem-per-cpu 5G
#SBATCH --array=1-97
#SBATCH -o /home/crq857/projects/Chapter2/slurmout/keepauto_%A_%a.out
#SBATCH -e /home/crq857/projects/Chapter2/slurmout/keepauto_%A_%a.err

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
BAM=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/crq857/projects/Chapter2/files/listbam_autosomes.txt | cut -f1)
NAME=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/crq857/projects/Chapter2/files/listbam_autosomes.txt | cut -f2)

OUTDIR=/projects/mjolnir1/people/crq857/Chapter2/bams_autosomes
BEDIN=/home/crq857/projects/Chapter2/files/chrlist_autosomes.bed

module load bedtools

#bedtools intersect -abam ${BAM} -b ${BEDIN} > ${OUTDIR}/${NAME}_autosomes.bam

module load samtools

samtools index ${OUTDIR}/${NAME}_autosomes.bam
