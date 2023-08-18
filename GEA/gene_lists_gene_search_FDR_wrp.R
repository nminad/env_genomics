# 1. get gene names for each variable
# 2. Find specific genes based on gene ID
# 3. Write a presence/absence table for genes and accessions
# Nikos Minadakis
# 9/9/2022 - 3/2/2023

library(tidyverse)
library(readr)
library(stringr)
setwd("/2_output/FDR/bedtools_intersect")

# 1. convert and only get the gene names

## insert files
f <- list.files(full.names = T)
f_short <- list.files(full.names = F)

for (i in 1:length(f)) {
  
  gene_list <- read_delim(f[i], col_names = F)
  
  ## give names to columns
  colnames(gene_list) <- c("chr", "start", "end", "chr2", "source", "type", "pos", "un1", "un2", "un3", "un4", "ID")
  
  ## if tibble 0x0 move to the next variable
  k = 0
  for (b in 1:nrow(gene_list)) {
    if (nrow(gene_list) != 0) {
      k <- k + 1
    }
  }
  
  if (k == 0) {
    next
  }
  
  
  ## only keep the ones interesting
  gene_list <- gene_list %>%
    filter(type == "gene") %>%
    select("ID")
  
  ## if there are no genes move to the next variable
  if (nrow(gene_list) == 0) {
    next
  }
  
  ## split by ";"
  split_equal <- as_tibble(gene_list %>%
    str_split(";", simplify = T))
  
  ## transform the vector to rows of one column and only select the ones including "Name"
  split_equal_intermediate = pivot_longer(split_equal, cols = 1:ncol(split_equal)) %>%
    filter(grepl("Name", value)) %>%
    select(value)
  
  ## split by ","
  split_gene_intermediate <- as_tibble(split_equal_intermediate %>%
    str_split(",", simplify = T))
  
  ## pivot longer and only get rows that include "ID" and the column "value"
  split_gene_intermediate_2 <- split_gene_intermediate %>%
    pivot_longer(cols = 1:ncol(split_gene_intermediate)) %>%
    filter(!grepl("ID", value)) %>%
    select(value)
  
  ## split by "="
  spl_gene_new <- as_tibble(split_gene_intermediate_2 %>%
    str_split("=", simplify = T))
  
  spl_gene_new2 <- spl_gene_new[,2:ncol(spl_gene_new)]
  
  spl_gene_names <- NA
  
  for (j in 1:ncol(spl_gene_new2)) {
    spl_gene_names[j] <- substr(spl_gene_new2[1,j], 1, 12)
  }
  spl_gene_names <- rbind.data.frame(spl_gene_names)
  
  if (length(spl_gene_names) == 1) {
    write.table(spl_gene_names[,1],
                paste("../gene_list_per_variable/spl_gene_names_",f_short[i], sep = ""),
                col.names = F, row.names = F, quote = F)
  } else {
    write.table(spl_gene_names[,2:length(spl_gene_names)],
              paste("../gene_list_per_variable/spl_gene_names_",f_short[i], sep = ""),
              col.names = F, row.names = F, quote = F)
  }
}


# 2. find specific genes based on gene ID

setwd('../2_output/FDR0.005/gene_list_per_variable/')

# we use the gene_list per variable as input

## insert files
f <- list.files(full.names = T)
f_short <- list.files(full.names = F)

## get a vector with all gene names combined
gene_names <- data.frame(matrix(0, ncol = 0, nrow = 0))
`%notin%` <- Negate(`%in%`)
for (i in 1:length(f)) {
  gene_list <- read_delim(f[i], col_names = F, delim = " ")
  gene_names_new <- pivot_longer(gene_list, cols = 1:length(gene_list)) %>%
    select(value)
  gene_names <- rbind(gene_names, gene_names_new)
}
gene_names <- unique(gene_names)

## make a dataframe with the right dimensions
## length of all gene names, length(f)
common_genes_df <- data.frame(matrix(0, ncol = length(f), nrow = nrow(gene_names)))

## combine everything
common_genes_df <- cbind(gene_names,common_genes_df)
gene_count <- rep(0, length(common_genes_df))
common_genes_df[nrow(common_genes_df) + 1,] <- gene_count

for (l in 1:length(f)) {
  gene_list <- read_delim(f[l], col_names = F, delim = " ")
  gene_names_new <- pivot_longer(gene_list, cols = 1:length(gene_list)) %>%
    select(value)
  counter_row <- 0
  for (k in 1:nrow(common_genes_df)) {
    if (common_genes_df[k, 1] %in% gene_names_new$value) {
      common_genes_df[k, l + 1] <- 1
      counter_row <- counter_row + 1
    } else {
      common_genes_df[k, l + 1] <- 0
    }
  }
  common_genes_df[nrow(common_genes_df),l + 1] <- counter_row
}

# create contig names and put them as colnames in my matrix
var_names_list = 0
for (i in 1:length(f_short)) {
  var_names_list[i] <- gsub(".*ID_", "contig_", f_short[i])
}
for (i in 1:length(f_short)) {
  var_names_list[i] <- gsub("_output.*", "", var_names_list[i])
}

colnames(common_genes_df) <- c("genes",var_names_list)
common_genes_df$genes[nrow(common_genes_df)] <- "count"

# 3. Write presence/absence table

write.table(common_genes_df, "../complete_gene_list_FDR.csv", row.names = F, col.names = T, quote = F, sep = ",")



  
  
  
  
