#!/bin/bash

vcftools --gzvcf ../0_data/Yann_filt_no_missing_noHet.vcf.gz --plink --out ../0_data/plinkfile

Rscript plinkfile_old_to_new_wrp.R

plink1.9 --file ../0_data/plinkfile --make-bed --out ../0_data/binary --noweb

Rscript binary_fam_old_to_new_wrp.R

bash gemma_loop_wrp.sh

mkdir ../2_output/gemma_output

mv ./output/* ../2_output/gemma_output/

#FDR
Rscript regions_of_significance_gemma_loop_FDR_wrp.R

bash ./bedtools_intersect_loop_FDR_wrp.sh

Rscript gene_lists_gene_search_FDR_wrp.R

Rscript Loop_for_plots_GWAS_wrp.R

Rscript UpSet_wrp.R
