#!/usr/bin/env bash
#SBATCH --job-name=pulloutgenes
#SBATCH -c 1
#SBATCH --time 01:30:00
#SBATCH --mem-per-cpu 1G
#SBATCH -o /home/crq857/projects/TSHRpaper/slurmout_Oct/pulloutgenes_filter.out
#SBATCH -e /home/crq857/projects/TSHRpaper/slurmout_Oct/pulloutgenes_filter.err

module load perl
module load bcftools/1.19

VCF=/projects/mjolnir1/data/user_owned_folders/Shyam_sharing/refpanel_1715g/1715g_WildSled.SNP.INDEL.chrAll.vcf.gz
OUTDIR=/projects/mjolnir1/people/crq857/TSHR_Oct2024/01_Reprogenes_scan
BED=/home/crq857/projects/reproductivegenes/1_dataset/Oct23_genes_sorted.bed
OUTVCF=/projects/mjolnir1/people/crq857/TSHR_Oct2024/01_Reprogenes_scan/1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23.vcf

bcftools view --regions-file ${BED} ${VCF} > ${OUTVCF}


module load perl
module load vcftools

vcftools --vcf ${OUTVCF} --minDP 3 --max-missing 0.8 --non-ref-ac-any 1 --out /projects/mjolnir1/people/crq857/TSHR_Oct2024/01_Reprogenes_scan/1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8 --recode --recode-INFO-all
