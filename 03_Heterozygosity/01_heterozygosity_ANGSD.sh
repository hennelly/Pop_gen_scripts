#!/usr/bin/env bash
#SBATCH --job-name=het
#SBATCH -c 1
#SBATCH --time 1-12:00:00
#SBATCH --mem-per-cpu 10G
#SBATCH -o /home/crq857/projects/Chapter2/slurmout/het_Dec26_sw.out
#SBATCH -e /home/crq857/projects/Chapter2/slurmout/het_Dec26_sw.err

BAM=wSierraMorena_autosomes.bam
ANC=/projects/mjolnir1/people/crq857/Chapter2/04_Demographichistory/ANGSD/AndeanFox_mapped2canfam31.fa.gz
REF=/projects/mjolnir1/people/crq857/Chapter2/ref/canFam31.fasta

## You need to make the ancestral fasta to fun this script in ANGSD. This would be the outgroup bam file, which you can convert to the fa.gz file. I'll try 
# find my code to convert it! The script should be within the program ANGSD to generate the ancestral fasta. 
module load angsd 
module load samtools 

samtools index /projects/mjolnir1/people/crq857/Chapter2/bams_auto//wSierraMorena_autosomes.bam

samtools faidx ${ANC}
#this makes the saf.idx file in ANGSD:
angsd -i /projects/mjolnir1/people/crq857/Chapter2/bams_auto/${BAM} -ref ${REF} -remove_bads 1 -only_proper_pairs 1 -C 50 -setMinDepth 4 -doCounts 1 -GL 1 -P 4 -minMapQ 20 -minQ 20 -anc ${ANC} -doSaf 1 -out /projects/mjolnir1/people/crq857/Chapter2/04_Demographichistory/ANGSD_Dec26/${BAM}_het
