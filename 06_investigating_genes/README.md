## Scripts for investigating structural variations and genes of interests 

### Investigating per-SNPs differences within selected genes 
 - Extract genes from VCF and filter: 01_extract_genes.sh
      - For my analyses, I kept indels in my VCF to better detect small structural variants, and converted any individual's site that was under 3x as missing (--minDP 3)
 - Per-SNP scan: 02_perSNP_Fst_scan.sh
      - I then looked for which individual sites were in the top 0.01% of all the SNPs I inferred Fst for. I took the top 1% or 0.1% or so to run SNPEff on.
 - Running SNPeff to infer functional impact: 03_runSNPeff.sh
 - You can also explore further interesting regions by doing a window-based fst and dxy scan : 07_fstdxy_10kbwindow.sh
### Investigating presence of structural variants with per-site depth 

- inferring the per-site depth is helpful for detecting deletions or the starting points of an inversion, in which the reference genome has a certain sequence, and the resequenced individuals mapping to the reference genome do not have reads with those sequences.
   - 04_persite_depth.sh
     - I inferred per-site depth for areas of interest, like from the Fst and SNpeff analyses where there could be interesting variants. Otherwise, if you know the approximate area of the inversion breakpoints, you could select this region and plot per-site depth to better detect the starting and end points. 
     
### Visualizing reads spanning the structural variants with IGV
   - IGV takes a BAM input and the reference genome (.fa file) to visualize reads.
   - To extact the region of interest, you can use: 05_extractBAM_for_IGV.sh
   - Now you can use IGV to see exactly where and aspects of the reads that correspond to the structural variant.
          - For example, with a deletion, often reads are "soft clipped" (if you click on a read, it will tell you if its soft-clipped) where part of the reads maps, and part of it doesn't. In addition, reads may have large insert sizes, which will span the entire deletion as part of the read maps to one side of the reference genome before the deletion, and the other read maps to the other side of the deletion on the reference genome.
          - I suspect for inversion, you will have something similar. Where part of the read will map to the reference genome right before the starting point of the inversion, and then the reads will be clipped with an abrupt zereo depth after that point. 

## Visualizing the gene in R 
- This script will plot the gene or region by color coating the reference and alternative alleles: 06_visualize_the_gene.sh
    - This can be useful to visually see difference in haplotypes 
