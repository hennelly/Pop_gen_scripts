
To infer the historical demography of wolves in Eurasia, we used the Pairwise
847 Sequential Markovian Coalescent (PSMC; Li and Durbin 2011). PSMC uses a
848 coalescent approach to estimate the history of change in effective population sizes over
849 time. We only included the autosomal sequences of each gray wolf individual. We
850 converted each bam file to a fasta-like consensus sequence by first using the mpileup
851 command with SAMtools and subsequently using BCFtools view –c to call variants and
852 vcfultils.pl vcf2fq to convert the vcf file to fastq format (Danecek et al. 2021). We
853 excluded any reads that were less than 20 for minimum mapping quality and minimum
854 base quality (-q 20 -Q 20) and excluded reads with excessive mismatches (-C 50). We
855 also removed sites with more than double or less than a third of the average depth of
coverage for each sample. We tested different combinations of parameters to infer the
857 PSMC, which were “psmc -N25 -t15 -r5 -p 4+25*2+4+6”, “psmc -N25 -t15
858 -r5 -p 2+2+25*2+4+6” and “psmc -N25 -t15 -r5 -p 1+1+1+1+25*2+4+6”
859 following previous studies on gray wolves and recent recommendations for inferring
860 PSMC trajectories (Hilgers et al. 2024). Some samples showed a false peak when using
861 the parameter “psmc -N25 -t15 -r5 -p 4+25*2+4+6”, therefore, we used the
862 parameter “psmc -N25 -t15 -r5 -p 1+1+1+1+25*2+4+6” for our final analyses
863 (Figure S25).
864
865 To account for our selected low coverage genomes (15-20x), we estimated the false
866 negative rate (FNR) by a downsampling high coverage gray wolf genome that were in
867 the same geographic region to the specific depth of the low coverage genome. To
868 determine the best FNR, we visually estimated the best correspondence between the
869 PSMC plots with the high coverage (>20x) regional wolf genome and downsampled
870 genomes with various FNR corrections (Figure S26). We then applied the best
871 estimated FNR to the low coverage gray wolf genomes to infer their demographic
history. We used a mutation rate of 4.5 x 10-9 872 (Koch et al. 2019) and a generation time
873 of 4.4 years (Mech and Barber-Meyer 2017).
