# Find regions of significance based on p-values

library(rehh)
library(hdi)
library(tidyverse)

setwd("/media/nikos/LaCie/Nikos_data/GEA/GEA_Yann_noHe/2_output/gemma_output/")

# Open and read the file list including all the variable names
f_short <- list.files(full.names = F, pattern = "assoc")
options(scipen = 100, digits = 4)
# Start the loop
for (i in 1:length(f_short)) {
  
  # Create a filename
  gwscanfile <- paste(f_short[i], sep = "")
  
  # Import GEMMA results
  gwscan1 <- read.table(gwscanfile, as.is = "rs", header = TRUE)
  
  # exclude rows with NaN
  gwscan <- gwscan1[complete.cases(gwscan1), ]

  # Threshold
  FDR <- function(pvals, FDR){
    pvalss <- sort(pvals, decreasing = F)
    m = length(pvalss)
    cutoffs <- ((1:m)/m)*FDR
    logicvec <- pvalss <= cutoffs
    postrue <- which(logicvec)
    k <- max(c(postrue,0))
    cutoff <- (((0:m)/m)*FDR)[k + 1]
    return(cutoff)
  }
  
  fdrlog10 <- -log10(FDR(gwscan$p_lrt, 0.005))
  
  BC <- -log10(0.05/nrow(gwscan))
  
  # Convert the p-values to the -log10 scale
  gwscan <- transform(gwscan, p_lrt = -log10(p_lrt))
  
  # Convert table to rehh format
  gwscan_rehh <- cbind.data.frame(gwscan$chr, gwscan$ps, gwscan$p_lrt)
  colnames(gwscan_rehh) <- c("CHR", "POSITION", "LOGPVALUE")
  
  # Delineate candidate regions
  significant_regions_gwscan = calc_candidate_regions(
    gwscan_rehh,
    threshold = fdrlog10,
    pval = TRUE,
    ignore_sign = FALSE,
    window_size = 5000,
    overlap = 2500,
    right = TRUE,
    min_n_mrk = 4,
    min_n_extr_mrk = 2, # change this to two if you want at least 2 SNPs to be significant for a certain area to be extracted
    min_perc_extr_mrk = 0,
    join_neighbors = TRUE)
  
  # Create folder if it does not exist
  my_folder1 <- "../FDR"
  if (!file.exists(my_folder1)) {
    dir.create(my_folder1)
  }
  my_folder2 <- "../FDR/rehh_output"
  if (!file.exists(my_folder2)) {
    dir.create(my_folder2)
  }
  my_folder3 <- "../FDR/Bd_chr"
  if (!file.exists(my_folder3)) {
    dir.create(my_folder3)
  }
  my_folder4 <- "../FDR/bedtools_intersect"
  if (!file.exists(my_folder4)) {
    dir.create(my_folder4)
  }
  my_folder5 <- "../FDR/gene_list_per_variable"
  if (!file.exists(my_folder5)) {
    dir.create(my_folder5)
  }
  my_folder6 <- "../QQ_MH_plots"
  if (!file.exists(my_folder6)) {
    dir.create(my_folder6)
  }
  # Write table
  write.table(significant_regions_gwscan, paste("../FDR/rehh_output/", f_short[i], ".FDR.sig.txt", sep = ""))
  bedtools_file <- significant_regions_gwscan %>%
    select(chrom = CHR, start = START,end = END)
  bedtools_file$chrom <- sub("^","Bd",bedtools_file$chrom)
  # write table as bed file
  write.table(bedtools_file, paste("../FDR/Bd_chr/", f_short[i], ".FDR.sig.bed", sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")
}
write.table(f_short, paste("../../0_data/file_list.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
