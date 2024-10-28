#!/usr/bin/env bash
#SBATCH --job-name=snpeff
#SBATCH -c 1
#SBATCH --time 1:00:00
#SBATCH --mem-per-cpu 1G
#SBATCH -o /home/crq857/projects/TSHRpaper/slurmout_Oct/01_snpeff_final.out
#SBATCH -e /home/crq857/projects/TSHRpaper/slurmout_Oct/01_snpeff_final.err

##################################################
#### STEP FIVE: SELECT SNPS for SNPEFF ON VCF ####
##################################################

#pull out SNPs with fst information

VCF=/projects/mjolnir1/people/crq857/TSHR_Oct2024/01_Reprogenes_scan/1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8.recode.vcf
BED=/home/crq857/projects/TSHRpaper/files_Oct2024/top_584_variants_Fst.txt
OUT=/projects/mjolnir1/people/crq857/TSHR_Oct2024/01_Reprogenes_scan/1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8_584topvariants.vcf
module load perl
module load vcftools

vcftools --vcf ${VCF} --positions ${BED} --recode --out ${OUT}

#pull out top 58 SNPs with fst information

VCF=/projects/mjolnir1/people/crq857/TSHR_Oct2024/01_Reprogenes_scan/1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8_584topvariants.vcf.recode.vcf
BED=/home/crq857/projects/TSHRpaper/files_Oct2024/top_58_variants_Fst.txt
OUT=/projects/mjolnir1/people/crq857/TSHR_Oct2024/01_Reprogenes_scan/1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8_58_topvariants_for_plotting.vcf
module load perl
module load vcftools

vcftools --vcf ${VCF} --positions ${BED} --recode --out ${OUT}


############################
#### RUN  SNPEFF ON VCF ####
############################

# pull out top 0.001% of SNPs to annotate

conda activate /projects/mjolnir1/apps/conda/snpeff-5.2

VCF=/projects/mjolnir1/people/crq857/TSHR_Oct2024/01_Reprogenes_scan/1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8_584topvariants.vcf.recode.vcf
CONFIG=/home/crq857/projects/reproductivegenes/scripts/snpEff_copy.config
OUT=/projects/mjolnir1/people/crq857/ReproGenes/vcffiles_Dec12/1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8_584topvariants_final.vcf

java -Xmx8g -jar /projects/mjolnir1/apps/snpeff-5.1d/snpEff.jar -c ${CONFIG} -v CanFam3.1.99 ${VCF} > ${OUT}


## take first couple columns

grep -v "#" 1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8_584topvariants_final.vcf > 1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8_584topvariants_final_nohash.vcf
cut -f 1-10 1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8_584topvariants_final_nohash.vcf  > 1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8_584topvariants_final_nohash_final.vcf

scp -r crq857@mjolnirgate.unicph.domain:/projects/mjolnir1/people/crq857/ReproGenes/vcffiles_Dec12/1715g_WildSled.SNP.INDEL.chrAll_candidategenes_Oct23_minDP3_maxmiss0.8_584topvariants_final_nohash_final.vcf ~/Desktop
