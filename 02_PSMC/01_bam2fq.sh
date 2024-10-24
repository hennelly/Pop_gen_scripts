#!/usr/bin/env bash
#SBATCH --job-name=PSMC
#SBATCH -c 1
#SBATCH --time 1-12:00:00
#SBATCH --mem-per-cpu 1G
#SBATCH --array=1
#SBATCH -o /home/crq857/projects/Chapter2/slurmout/PSMC_De22_test_bam2fq_%A_%a.out
#SBATCH -e /home/crq857/projects/Chapter2/slurmout/PSMC_De22_test_bam2fq_%A_%a.err

## Bam files were only for the autosomes
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
BAM=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/crq857/projects/Chapter2/list_Dec/psmc_lc_list_Shanxi.txt | cut -f1)
DIR=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/crq857/projects/Chapter2/list_Dec/psmc_lc_list_Shanxi.txt | cut -f2)
MIND=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/crq857/projects/Chapter2/list_Dec/psmc_lc_list_Shanxi.txt | cut -f3) #minimum depth of 0.5 of the average depth per sample
MAXD=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/crq857/projects/Chapter2/list_Dec/psmc_lc_list_Shanxi.txt | cut -f4) #maximum depth of 
OUTDIR=/projects/mjolnir1/people/crq857/Chapter2/04_Demographichistory/PSMC/fq_files
REF=/projects/mjolnir1/people/crq857/Chapter2/ref/canFam31.fasta

module load perl
module load gsl/2.5
module load samtools
module load bcftools

bcftools mpileup -q 20 -Q 20 -C 50 -f ${REF} ${DIR}/${BAM}  | bcftools call -c | vcfutils.pl vcf2fq -d ${MIND} -D ${MAXD} | gzip > ${OUTDIR}/${BAM}.fq.gz 

