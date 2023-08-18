## UpSetR

library(tidyverse)
library(UpSetR)

setwd("/media/nikos/LaCie/Nikos_data/GEA/GEA_Yann_noHe/2_output/FDR/bedtools_intersect")

# I have to remove the rows that do not have 0s or 1s

common_genes_df <- read.csv("../complete_gene_list_FDR.csv")
test_upset <- common_genes_df[c(1:nrow(common_genes_df) - 1), c(2:length(common_genes_df))]

pdf("../UpSet_plot.pdf", width = 12, height = 6)
upset(test_upset, nsets = 200, line.size = 1, nintersects = NA,
      point.size = 1.5, text.scale = .6, mb.ratio = c(0.7, 0.3),
      sets.bar.color = "dodgerblue3", set_size.show = T, order.by = "freq")
dev.off()

# specific set
my_set <- test_upset %>%
  select(c("spl_gene_names_XtX_statistics_top0.1perc_bedtools_intersect_wb_out.txt",
           "spl_gene_names_ULMM_bio_12",
           "spl_gene_names_ULMM_bio_16",
           "spl_gene_names_ULMM_bio_13",
           "spl_gene_names_ULMM_prec_spring",
           "spl_gene_names_ULMM_aridity_spring",
           "spl_gene_names_ULMM_prec_winter",
           "spl_gene_names_ULMM_bio_18"))
pdf("../UpSet_morethan100genes.pdf", width = 18, height = 6)
upset(my_set, nsets = length(my_set), line.size = 1, point.size = 3, sets.bar.color = "dodgerblue3", set_size.show = T,
      text.scale = 2.8, order.by = "freq", nintersects = 34, mb.ratio = c(0.4,0.6))
dev.off()

