# Order phenotype to 332 GWAS binary.fam

library(tidyverse)

setwd("/media/nikos/LaCie/Nikos_data/Lineage_Clade_vcf/noHe_Analysis/Clade_A_East/0_data")

# File name to rename
filename <- 'binary.fam'

# Check if file exists.
if (file.exists(filename)) {
  # Rename file name
  file.rename(filename,'binary_old.fam')
}else{
  print('File Not found :')
}

options(scipen = 100, digits = 4)

# import data
originalfam <- read.table('binary_old.fam')
originalfam$V6 <- -9
phenofile <- read.table('Bdis332_SNPs_Yann_final_noHet_binary.fam')

filtered_phenofile <- phenofile %>%
  filter(phenofile$V1 %in% originalfam$V1)



#list_of_missing <- originalfam %>%
#  filter(!(originalfam$V1 %in% phenofile$V1))
#list_of_missing <- list_of_missing[c(1, 5, 6)]
#colnames(list_of_missing) <- c("V1", "V2", "V3")

# fill the phenofile until it has the same number of rows
#phenofile <- rbind(phenofile, list_of_missing)

# match my_data according to originalfam names order
ordered_phenofile <- filtered_phenofile[match(originalfam$V1, filtered_phenofile$V1),]

famfile <- ordered_phenofile

write.table(famfile, "binary.fam", sep = " ", row.names = F, col.names = F, quote = F)

varnames_tbl <- read.table('Bdis332_SNPs_Yann_final_noHet_binary_varNames.fam', header = T)

varnames <- cbind(colnames(varnames_tbl[6:length(varnames_tbl)]), 1:(length(famfile) - 5))

write.table(varnames, "binary_varNames.fam", sep = ",", row.names = F, col.names = F, quote = F)

