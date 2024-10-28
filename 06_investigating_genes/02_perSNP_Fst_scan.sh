#!/usr/bin/env bash
#SBATCH --job-name=fst
#SBATCH -c 1
#SBATCH --time 1:00:00
#SBATCH --mem-per-cpu 1G
#SBATCH -o /home/crq857/projects/TSHRpaper/slurmout_Oct/01_fstscan.out
#SBATCH -e /home/crq857/projects/TSHRpaper/slurmout_Oct/01_fstscan.err

############################
#### STEP ONE: RUN FST ####
############################

VCF=/projects/mjolnir1/people/crq857/TSHR_Oct2024/01_Reprogenes_scan/1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8.recode.vcf
dogpop=/home/crq857/projects/reproductivegenes/1_dataset/keepdogs_Oct242024.txt
wolfpop=/home/crq857/projects/reproductivegenes/1_dataset/keepwolves_Sept212024.txt
OUT=/projects/mjolnir1/people/crq857/TSHR_Oct2024/01_Reprogenes_scan/01_Fst_results/1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8_fst

module load perl
module load vcftools

vcftools --vcf ${VCF} \
--weir-fst-pop ${dogpop} \
--weir-fst-pop ${wolfpop} \
--out ${OUT}

#Outputs 1241518 per-SNP fst with -nans

####################
### REMOVE -NAN ####
####################
DIR=/projects/mjolnir1/people/crq857/TSHR_Oct2024/01_Reprogenes_scan/01_Fst_results

sed '/n/d' ${DIR}/1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8_fst.weir.fst >  ${DIR}/1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8_noNAN.weir.fst

scp -r crq857@mjolnirgate.unicph.domain:/projects/mjolnir1/people/crq857/TSHR_Oct2024/01_Reprogenes_scan/01_Fst_results/Autosomes_Canis_modern_v9_candidategenes_minDP3_maxmiss0.8_fst_noNAN_header.weir.fst.txt ~/Desktop
