# env_genomics
Scripts used in: "New resources for environmental genomics in the wild Mediterranean grass B. distachyon", https://doi.org/10.1101/2023.06.01.543285

## pop_genomics

**LEA_PCA.R**

Create Structure plots, and calculate optimal k value with Evanno.

**estimate_pop_size.sh**

Estimate population size with Relate.

**Calculate_pi.sh**

Calculate pi values using pixy.

**Tajimas_d_from_pixy_pi.R**

Calculate Tajima's D values from pixy.

**bpp_CS_RH.sh**

Run bpp from bam files and plot a tree.

**gemma_loop.sh**

Loop for GEAs using the software GEMMA.

**regions_of_significance_gemma_loop_FDR_wrp.R**

Extract regions of significance based on p-values above the False Discovery Rate threshold.

**bedtools_intersect_loop_FDR_wrp.sh**

Extract the significant genes from the v3.2 annotation file of B. distachyon, based on the significant regions.

**gene_lists_gene_search_FDR_wrp.R**

Use the significant genes and create a table of environmental variables (columns) and significant genes (rows) with 0s (absence) and 1s (presence) as values.

**UpSet.R**

Create an Upset plot using the gene list table.

**gene_age_statistics.R**

Gene age statistics based on significant genes for each variable. Boxplots for gene ages.
