#!/usr/bin/env bash
#SBATCH --job-name=intersectvcf
#SBATCH -c 1
#SBATCH --time 2:00:00
#SBATCH --mem-per-cpu 5G
#SBATCH -o /home/crq857/projects/TSHRpaper/slurmout/vcf2plink2vcf.out
#SBATCH -e /home/crq857/projects/TSHRpaper/slurmout/vcf2plink2vcf.err

conda activate /projects/mjolnir1/apps/conda/plink-1.90b6.21

# Convert to plink to plot better in R
OUTVCF_filtered=/projects/mjolnir1/people/crq857/ReproGenes/vcffiles_Dec12/TSHR_1715g_WildSled_SNP_INDEL_candidategenes_Dec12_min3DP_maxmissing0.8_snpeff_May92024_SNPEFF_highfst.vcf.recode.vcf
OUTBED=/projects/mjolnir1/people/crq857/ReproGenes/vcffiles_Dec12/TSHR_1715g_WildSled_SNP_INDEL_candidategenes_Dec12_min3DP_maxmissing0.8_snpeff_May92024_SNPEFF_highfst_bed
OUTPLINK=/projects/mjolnir1/people/crq857/ReproGenes/vcffiles_Dec12/TSHR_1715g_WildSled_SNP_INDEL_candidategenes_Dec12_min3DP_maxmissing0.8_snpeff_May92024_SNPEFF_highfst_plink
OUTVCF_final=/projects/mjolnir1/people/crq857/ReproGenes/vcffiles_Dec12/TSHR_1715g_WildSled_SNP_INDEL_candidategenes_Dec12_min3DP_maxmissing0.8_snpeff_May92024_SNPEFF_highfst_final

plink --vcf ${OUTVCF_filtered} --make-bed --const-fid --dog --out ${OUTBED} #to sort the chromosomes
plink --bfile ${OUTBED} --recode --out ${OUTPLINK} --const-fid --dog
plink --file  ${OUTPLINK} --recode vcf --out ${OUTVCF_final} --const-fid --dog


##### Convert vcf to coding ######
sed -i 's:0/0:0:g' 1715g_WildSled.SNP.INDEL.chrAll_minDP3_maxmissing0.8_noindels_bialleleic_minQ30_dogswolves_chr8_TSHRgene_finalcoded.vcf
sed -i 's:0/1:1:g' 1715g_WildSled.SNP.INDEL.chrAll_minDP3_maxmissing0.8_noindels_bialleleic_minQ30_dogswolves_chr8_TSHRgene_finalcoded.vcf
sed -i 's:1/0:1:g' 1715g_WildSled.SNP.INDEL.chrAll_minDP3_maxmissing0.8_noindels_bialleleic_minQ30_dogswolves_chr8_TSHRgene_finalcoded.vcf
sed -i 's:1/1:2:g' 1715g_WildSled.SNP.INDEL.chrAll_minDP3_maxmissing0.8_noindels_bialleleic_minQ30_dogswolves_chr8_TSHRgene_finalcoded.vcf
sed -i 's:./.:3:g' 1715g_WildSled.SNP.INDEL.chrAll_minDP3_maxmissing0.8_noindels_bialleleic_minQ30_dogswolves_chr8_TSHRgene_finalcoded.vcf

scp -r crq857@mjolnirgate.unicph.domain:/projects/mjolnir1/people/crq857/TSHR_Oct2024/01_Reprogenes_scan/chr8_Canis_modern_v9_TSHR.recode.vcf ~/Desktop


## PLOTTING IN R 

dat_hap <- read.csv("Oct26_Figure1_TSHR_Haplotype_phased.csv", header=TRUE)

library (ggplot2)
library(dplyr)
library (tidyverse)

## Organize the dataset
data_long <- gather(dat_hap, Sample, Allele, Affenpinscher01_33500_A:Wolf95_B, factor_key=TRUE)

#Make my variables in factor format
data_long$POS <- as.factor(data_long$POS)
data_long$Allele <- as.factor(data_long$Allele)

#Plot with geom_tile entire TSHR gene
p1 <- ggplot(data_long,aes(x=POS,y=Sample,fill=Allele))+
      geom_tile() 
p2 <- p1 + scale_fill_manual(values = c("lightblue", "red", "black", "black", "black")) 

p3 <- p2 + scale_x_discrete(breaks=c("53394526","53413008"))

p3


+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_x_discrete(expand=c(0,0),
                   breaks=c("53394526","53413008"))
p2 + scale_x_discrete(position = "top") 





ggsave("Oct26_TSHR.tiff", width=10,height=3) 



