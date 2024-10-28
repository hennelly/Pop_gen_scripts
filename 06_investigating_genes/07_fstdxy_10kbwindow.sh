#!/usr/bin/env bash
#SBATCH --job-name=dxypifst
#SBATCH -c 1
#SBATCH --time 02:00:00
#SBATCH --mem-per-cpu 5G
#SBATCH -o /home/crq857/projects/TSHRpaper/slurmout_Oct/xxTSHR_Oct25_%A_%a.out
#SBATCH -e /home/crq857/projects/TSHRpaper/slurmout_Oct/xxTSHR_Oct25_%A_%a.err
#SBATCH --array=2-38

#30,31,24,29

module load python/2.7.18
module load numpy
module load perl
module load vcftools 
module load bcftools

CHR=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /home/crq857/projects/Chapter2/files/chrlist.txt | cut -f1)

#conda activate /projects/mjolnir1/apps/conda/numpy-1.21.2
DIR=/projects/mjolnir1/people/crq857/TSHR_Oct2024/01_Reprogenes_scan
DIRGENO=/projects/mjolnir1/people/crq857/TSHR_Oct2024/04_dxyfst/genofiles

## Convert to Geno
python3 /projects/mjolnir1/people/crq857/bin/genomics_general/VCF_processing/parseVCF.py -i ${DIR}/1715g_WildSled.SNP.INDEL.chrAll_minDP3_maxmissing0.8_noindels_bialleleic_minQ30_dogswolves_${CHR}.recode.vcf --skipIndels | gzip > ${DIRGENO}/1715g_WildSled.SNP.INDEL.chrAll_minDP3_maxmissing0.8_noindels_bialleleic_minQ30_dogswolves_${CHR}_geno.gz

POPFILE=/home/crq857/projects/TSHRpaper/files_Oct2024/dxyfst_dogswolves.popfile
OUTDIR=/projects/mjolnir1/people/crq857/TSHR_Oct2024/04_dxyfst/results

### Run Fst analysis
python3 /projects/mjolnir1/people/crq857/bin/genomics_general/popgenWindows.py -g ${DIRGENO}/1715g_WildSled.SNP.INDEL.chrAll_minDP3_maxmissing0.8_noindels_bialleleic_minQ30_dogswolves_${CHR}_geno.gz -f phased -o ${OUTDIR}/1715g_WildSled.SNP.INDEL.chrAll_minDP3_maxmissing0.8_noindels_bialleleic_minQ30_dogswolves_${CHR}_final_Oct25_2024.csv -m 50 -w 10000 -s 2500 -p dog -p wolf --popsFile ${POPFILE}
