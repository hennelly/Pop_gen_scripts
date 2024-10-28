#!/usr/bin/env bash
#SBATCH --job-name=selectTSHR
#SBATCH -c 1
#SBATCH --time 01:10:00
#SBATCH --mem-per-cpu 1G
#SBATCH -o /home/crq857/projects/TSHRpaper/slurmout/aPortugueseWolf_%A_%a.out
#SBATCH -e /home/crq857/projects/TSHRpaper/slurmout/aPortugueseWolf_%A_%a.err
#SBATCH --array=1-56

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
DIR=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /projects/mjolnir1/people/crq857/SarabiaData/Bams/list_Jul1 | cut -f1)
BAM=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /projects/mjolnir1/people/crq857/SarabiaData/Bams/list_Jul1 | cut -f2)

echo ${BAM}

BED=/home/crq857/projects/reproductivegenes/scripts_Dec12/TSHR_deletion.bed

module load samtools 

samtools index ${DIR}/${BAM}

samtools view -bh -L ${BED} -o /projects/mjolnir1/people/crq857/TSHRpaper/InsertionAnalysis/Bams/${BAM}_TSHR.bam ${DIR}/${BAM}

samtools index /projects/mjolnir1/people/crq857/TSHRpaper/InsertionAnalysis/Bams/${BAM}_TSHR.bam



scp -r crq857@mjolnirgate.unicph.domain:/projects/mjolnir1/people/crq857/TSHRpaper/InsertionAnalysis/Bams/*.bam* ~/Desktop
