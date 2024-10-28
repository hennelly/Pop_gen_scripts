## Scripts for investigating structural variations and genes of interests 

### Investigating per-SNPs differences within selected genes 
 - Extract genes from VCF and filter: 01_extract_genes.sh
      - For my analyses, I kept indels in my VCF to better detect small structural variants, and converted any individual's site that was under 3x as missing (--minDP 3)
 - Per-SNP scan: 02_perSNP_Fst_scan.sh
      - I then looked for which individual sites were in the top 0.01% of all the SNPs I inferred Fst for. I took the top 1% or 0.1% or so to run SNPEff on.
 - Running SNPeff to infer functional impact:  
### Investigating presence of structural variants with per-site depth 
- 
