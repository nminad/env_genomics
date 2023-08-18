# Create sNMF plot, Evanno, PCA
# Nikos Minadakis - April 2023

library(SNPRelate)
library(LEA)
library(scales)
library(tidyverse)

setwd("/media/nikos/LaCie/Nikos_data/GEA/GEA_Yann_noHe/2_output/Admixture_LEA/0_data/")

# Only use to load a project --------------------
project = load.snmfProject(("Bdis332_Yann_final_noHet.gds.snmfProject"))

# Start here normally -------------
filename <- "0_data/Bdis332_SNPs_Yann_final_noHet.vcf.gz"

# convert VCF to SNP GDS

snpgdsVCF2GDS(filename, "Bdis332_Yann_final_noHet.gds", method = "biallelic.only", ignore.chr.prefix = "Bd")

snpgdsSummary("Bdis332_Yann_final_noHet.gds")

genofile <- openfn.gds("Bdis332_Yann_final_noHet.gds")

samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# LD prune
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.4, maf = 0.05)
snpset.id <- unlist(unname(snpset))
str(snpset)

# Get .map file with the SNPs that passed LD and MAF

snpgdsGDS2PED(genofile, "Bdis332_Yann_final_noHet.gds", snp.id = snpset.id)

ped2geno(input.file = "Bdis332_Yann_final_noHet.gds.ped")

write.table(samp.id, "acc_names.txt", quote = F, row.names = F, col.names = F)

# Run analysis
project = snmf(paste("Bdis332_Yann_final_noHet.gds", "geno", sep = "."),
               K = 2:20, 
               entropy = TRUE, 
               repetitions = 10,
               project = "new", 
               CPU = 4, ploidy = 2)

## CONTINUE HERE -------------
# Save plots to file

genofile <- openfn.gds("Bdis332_Yann_final_noHet.gds")

samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# plot cross-entropy criterion of all runs of the project
png('crossentropy_k2-20.png')
plot(project, cex = 1.2, pch = 19)
dev.off()

pdf('crossentropy_k2-20.pdf')
plot(project, cex = 1.2, pch = 19)
dev.off()

# Keep the Q.matrix for k = 5
best = which.min(cross.entropy(project, K = 5))
Q.matrix_k5 <- as.matrix((Q(project, K = 5, run = best)))
rownames(Q.matrix_k5) <- samp.id
write.table(Q.matrix_k5, "Q.matrix_k5.tbl", quote = F)

# PLOTS

my.colors <- c( "#96127d", "#5DC863FF", "#FDE725FF", "#3B528BFF", "#21908c")
show_col(my.colors)

# get the cross-entropy value for each run 
ce <- cross.entropy(project, K = 5)

# select the run with the lowest cross-entropy value
best <- which.min(ce)

pdf("/media/nikos/LaCie/Nikos_data/GEA/GEA_Yann_noHe/LEA_admixture_barchart.pdf")
barchart(project, K = 5, run = best,
         border = F, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las = 1,
     cex.axis = .4)
dev.off()



ce = cross.entropy(project, K = 5)
best = which.min(ce)
Q.matrix <- as.matrix(Q(project, K = 5, run = best))
rownames(Q.matrix) <- samp.id
colnames(Q.matrix) <- c("1", "2", "3", "4", "5")

# keep the ancestry Q shortage
# create a new table with the names and order of Q.matrix
Q.matrix.bp <- as.data.frame(cbind(row.names(Q.matrix), bp$order))
colnames(Q.matrix.bp) <- c("id", "order")
Q.matrix.bp <- Q.matrix.bp[order(as.numeric(as.character(Q.matrix.bp$order))),]
Q.matrix.ordered <- Q.matrix[match(row.names(Q.matrix), Q.matrix.bp$id),]

pdf("LEA_admixture_barplot.pdf", width = 16, height = 8)
barplot(t(Q.matrix.ordered), col = my.colors, las = 2, cex.names = 0.2, border = F)
dev.off()

pdf('LEA_k2-10_barplots.pdf')
par(mfrow = c(9,1),
    mar = c(0.2,1,0.2,1),
    oma = c(3,0,0,0), 
    lwd = 0.5)
for (k in 2:10) {
  ce = cross.entropy(project, K = k)
  best = which.min(ce)
  Q.matrix <- as.matrix(Q(project, K = k, run = best))
  rownames(Q.matrix) <- samp.id
  
  # keep the ancestry Q shortage
  # create a new table with the names and order of Q.matrix
  Q.matrix.bp <- as.data.frame(cbind(row.names(Q.matrix), bp$order))
  colnames(Q.matrix.bp) <- c("id", "order")
  Q.matrix.bp <- Q.matrix.bp[order(as.numeric(as.character(Q.matrix.bp$order))),]
  Q.matrix.ordered <- Q.matrix[match(row.names(Q.matrix), Q.matrix.bp$id),]
  if (k == 10) {
    barplot(t(Q.matrix.ordered), col = rainbow(k), las = 2, cex.names = 0.1)
  }
  else {
    barplot(t(Q.matrix.ordered), col = rainbow(k), las = 2, axisnames = F)
  }
}
dev.off()


#Evanno

# show the project
show(project)
# summary of the project  
summary(project)

# Plot cross entropy and delta cross entropy
library(ggplot2)
library(reshape2)
library(magrittr)

ce = cross.entropy(project, 2)

for (k in 3:20){
  ce <- cbind(ce, cross.entropy(project, k))
}

colnames(ce) <- as.factor(c(2:20))

ce_m <- melt(ce)
ce_m$Var2 <- as.factor(ce_m$Var2)

ggplot(ce_m) +
  geom_boxplot(aes(x=Var2, y=value)) +
  theme_bw(base_size=14, base_family = "Arial") +
  xlab("Ancestral populations") + ylab("Cross entropy")


# Evanno et al. 2005 for cross entropy rather than likelihoods

# Mean differences between successive CEs, CE'(K) = CE(K) - CE(K-1)
ce_diff <- data.frame(
  '2' = rep(NA, nrow(ce))
)

for (i in 2:ncol(ce)){
  diff = ce[,i] - ce[,i-1]
  ce_diff <- cbind(ce_diff, diff)
}

colnames(ce_diff) <- c(2:20)


# Second order differences CE''(K) = CE'(K+1) - CE'(K)
ce_diff2 <- data.frame(
  '3' = abs(ce_diff[,3] - ce_diff[,2])
)

for (i in 3:(ncol(ce_diff)-1)){
  diff = abs(ce_diff[,i+1] - ce_diff[,i])
  ce_diff2 <- cbind(ce_diff2, diff)
}

colnames(ce_diff2) <- c(3:19)


# Finally, delta K = mean |CE''(K)| / sd CE(K)
ce_pruned = ce[,c(2:18)]


dK = c()
for (k in 1:17){
  a = mean(ce_diff2[,k])  / sd(ce_pruned[,k])
  dK = c(dK, a)
}



dK_df <- data.frame(
  'K' = c(2:20),
  'value' = c(NA,dK, NA)
)

dK_df$step <- rep('deltaK', nrow(dK_df))

# Plot

comb <- ce_m[,c(2,3)]
comb$step <- rep('CE', nrow(comb))
colnames(comb) <- c('K', 'value', 'step')

ce_diff_m <- melt(ce_diff)[,c(2,3)]
ce_diff_m$step <- rep('CE_diff1', nrow(ce_diff))
colnames(ce_diff_m) <- c('K', 'value', 'step')


ce_diff2_m <- melt(ce_diff2)
ce_diff2_m$step <- rep('CE_diff2', nrow(ce_diff2))
colnames(ce_diff2_m) <- c('K', 'value', 'step')


comb <- rbind(comb, ce_diff_m, ce_diff2_m, dK_df)

pdf("Evanno.pdf", height = 8, width = 16)
ggplot(comb) +
  geom_boxplot(aes(x=K, y=value)) +
  theme_bw(base_size=14) +
  xlab("K") + ylab("")+
  facet_wrap(~step, scales='free_y') +
  theme_bw(base_size=12)
dev.off()

my_pca <- snpgdsPCA(genofile, snp.id = snpset.id)

closefn.gds(genofile) # do not forget to close the gds file

plot(my_pca$eigenvect[,1],my_pca$eigenvect[,2])
plot(my_pca$eigenvect[,2],my_pca$eigenvect[,3])

# calculate percentages
pc.percent <- my_pca$varprop*100
percentages <- round(pc.percent, 2)
percentages

my_pca_df <- data.frame(sample.id = my_pca$sample.id, # make data frame
                        PC1 = my_pca$eigenvect[,1], # 1st PC
                        PC2 = my_pca$eigenvect[,2], # 2nd PC
                        PC3 = my_pca$eigenvect[,3], # 3rd PC
                        PC4 = my_pca$eigenvect[,4], # 4th PC
                        PC5 = my_pca$eigenvect[,5], # 5th PC
                        PC6 = my_pca$eigenvect[,6], # 6th PC
                        PC7 = my_pca$eigenvect[,7], # 7th PC
                        stringsAsFactors = FALSE)

write.csv(my_pca_df, "pca_table_Bdis332_2023.csv")

info_clades <- read.csv("./CLADE_SOURCE_INFO.csv", header = T)

new_ordered <- as.data.frame(info_clades[match(my_pca_df$sample.id, info_clades$name),])
#new_ordered <- subset(new_ordered, (!is.na(new_ordered[,2])))

my_pca_df$clade <- new_ordered$cluster
my.colors <- c( "#96127d", "#21908c", "#5DC863FF", "#3B528BFF", "#FDE725FF")
show_col(my.colors)

my_pca_df$source <- new_ordered$source
my.colors.source <- c( "red", "green", "black", "blue", "azure4")
show_col(my.colors.source)

pc1_clade <- ggplot(my_pca_df, aes(x=PC1, y=PC2, label=sample.id, col = clade)) +
  geom_point(size=2, alpha = 0.8) +
  xlab(paste("PC1: ", percentages[1], "%", sep = "")) +
  ylab(paste("PC2: ", percentages[2], "%", sep = ""))

pdf("PC1-PC2_clade.pdf")
pc1_clade + scale_color_manual(values = my.colors) +
  theme_test()
dev.off()

pc1_source <- ggplot(my_pca_df, aes(x=PC1, y=PC2, label=sample.id, col = source)) +
  geom_point(size=2, alpha = 0.8) +
  xlab(paste("PC1: ", percentages[1], "%", sep = "")) +
  ylab(paste("PC2: ", percentages[2], "%", sep = ""))

pdf("PC1-PC2_source.pdf")
pc1_source + scale_color_manual(values = my.colors.source) +
  theme_test()
dev.off()


# rest of PCAs

pc2 <- ggplot(my_pca_df, aes(x=PC2, y=PC3, label=sample.id, col = clade)) +
  geom_point(size=2) +
  xlab(paste("PC2: ", percentages[2], "%", sep = "")) +
  ylab(paste("PC3: ", percentages[3], "%", sep = ""))

pdf("PC2-PC3_Lineage_A.pdf")
pc2 + scale_color_manual(values = my.colors) +
  theme_test()
dev.off()

pc3 <- ggplot(my_pca_df, aes(x=PC3, y=PC4, label=sample.id, col = clade)) +
  geom_point(size=2) +
  xlab(paste("PC3: ", percentages[3], "%", sep = "")) +
  ylab(paste("PC4: ", percentages[4], "%", sep = ""))

pdf("PC3-PC4_Lineage_A.pdf")
pc3 + scale_color_manual(values = my.colors) +
  theme_test()
dev.off()



