#!/usr/bin/env bash
#SBATCH --job-name=depth
#SBATCH -c 1
#SBATCH --time 12:00:00
#SBATCH --mem-per-cpu 1G
#SBATCH -o /home/crq857/slurmout/depthTSHR.out
#SBATCH -e /home/crq857/slurmout/depthTSHR.err

module load vcftools

VCF=/projects/mjolnir1/people/crq857/ReproGenes/vcffiles/1715g_WildSled_SNP_INDEL_TSHR_allsites_Dec15_minDP3_allcanids_village.vcf
OUT=/projects/mjolnir1/people/crq857/ReproGenes/depth/TSHR_per_genotype_village

vcftools --vcf ${VCF} --geno-depth --out ${OUT}


## Plotting in R

library (ggplot2)
library(dplyr)
library (tidyverse)
library(reshape2)
install.packages("reshape2")
library(plotly)
install.packages("plotly")

dat <- read.csv("TSHR_depth_per_site_Sept26.csv", header=TRUE )

#Convert wide format to long format:
data_long <- gather(dat2, Sample, Depth, AlaskanWolf:YorkshireTerrier76_14941, factor_key=TRUE)

#Pull out region of interest
dat1 <- dat[dat[,2] > 53394526,]
dat2 <- dat1[dat1[,2] < 53413008,]

#Set maximum depth to 50x 
data_long$Depth = ifelse(data_long$Depth > 50, 50, data_long$Depth)

#Make my variables in factor format
data_long$POS <- as.factor(data_long$POS)

#Plot
q <- ggplot(data_long, aes(x=POS,y=Sample,fill=Depth)) +
  geom_raster() + scale_fill_viridis_c()
p <- q + theme(axis.text.x = element_text(angle = 45, hjust = 1))  + theme(axis.text.y = element_text(size = 10))  
r <- p + scale_x_discrete(expand=c(0,0),
                   breaks=c("53394527","53413007")) + scale_y_discrete(expand=c(0,0),
                   breaks=c( "Dingo01","NewGuineaSingingDog17")) 

                   ggsave("Sept26_2024_TSHR_insertion_allcanids.tiff", width=5,height=7)
