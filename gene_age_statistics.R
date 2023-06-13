# Extract gene ages for each list
library(tidyverse)
library(nortest)

setwd("/media/nikos/LaCie/Nikos_data/GEA/GEA_Yann_noHe/2_output/FDR0.005/Gene_age_estimates/")

options(scipen = 999)

my_folder <- "./2_output/"
if (!file.exists(my_folder)) {
  dir.create(my_folder)
}

# table with all gene info
all_genes <- read.delim("../../../0_data/All_genes_age_matrix.txt")
head(all_genes)
colnames(all_genes) <- c("chr", "pos1", "pos2", "Mean_age", "Number_of_SNPs")
all_vars <- all_genes[,1:4]

# loop for all files in the folder

## insert files
f <- list.files("../bedtools_intersect", full.names = T)
f_short <- list.files("../bedtools_intersect", full.names = F)

for (y in 1:length(f)) {
  my_var <- read.delim(f[y], header = FALSE)
  head(my_var)
  colnames(my_var) <- c("chr", "start", "end", "chr2", "source", "type", "pos1", "pos2", "un1", "un2", "un3", "ID")
  my_var <- my_var %>%
    filter(type == "gene") %>%
    select("chr", "pos1", "pos2")
  mycolname <- paste(f_short[y], "_age", sep = "")
  all_genes[[mycolname]] <- NA # create new column with the current file name
  position_counter <- 1
  for (i in 1:nrow(my_var)) {
    for (j in position_counter:nrow(all_genes)) {
      if ((my_var$chr[i] == all_genes$chr[j]) & (my_var$pos1[i] == all_genes$pos1[j])) {
        all_genes[[mycolname]][j] <- all_genes$Mean_age[j]
        position_counter <- j + 1
        next
      } else {
        all_genes[[mycolname]][j] <- NA
      }
    }
  }
}

### simplify names
colnames(all_genes) = gsub("wc2.1_30s_", "",colnames(all_genes))
colnames(all_genes) = gsub("ULMM_", "",colnames(all_genes))
colnames(all_genes) = gsub("_output.assoc.txt_bedtools_intersect_FDR_wb_out.txt_age", "",colnames(all_genes))
colnames(all_genes) = gsub("_et0_", "_",colnames(all_genes))
colnames(all_genes) = gsub("spring", "March-June",colnames(all_genes))
colnames(all_genes) = gsub("winter", "Nov-Feb",colnames(all_genes))
colnames(all_genes) = gsub("_statistics_top0.1perc_bedtools_intersect_wb_out.txt_age", "",colnames(all_genes))

summary(all_genes)

write.table(all_genes, "./2_output/all_genes_and_variables.tbl", row.names = F, col.names = T, quote = F, sep = ",")
write.table(summary(all_genes), "./2_output/all_genes_and_variables_summary.tbl", row.names = T, col.names = T, quote = F, sep = ",")

ggplot(stack(all_genes[,c(4,6:length(all_genes))]), aes(x = ind, y = values)) +
  geom_boxplot()

biovars <- select(all_genes, c(4,contains("bio")))
ggplot(stack(biovars), aes(x = ind, y = values)) +
  geom_boxplot()

monthvars <- select(all_genes, c(4,contains(c("Feb", "June", "elevation"))))
ggplot(stack(monthvars), aes(x = ind, y = values)) +
  geom_boxplot()

# get number of genes per variable
colSums(!is.na(all_genes))

# only keep genes with more than 4 SNPs
trustworthy_dataset_SNPs <- filter(all_genes, Number_of_SNPs > 4)

# only keep variables with more than 14 genes
trustworthy_dataset <- select(trustworthy_dataset_SNPs, c(4,6:length(trustworthy_dataset_SNPs)))
trustworthy_dataset <- trustworthy_dataset[,colSums(!is.na(trustworthy_dataset)) > 14]

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

# boxplots only with trustworthy dataset
myboxplot <- ggplot(stack(trustworthy_dataset), aes(x = ind, y = values)) +
  geom_boxplot() +
  scale_y_continuous(breaks = round(seq(0, max(trustworthy_dataset$Mean_age, na.rm = T), by = 10000),1)) +
  xlab("") + ylab("gene age in years") +
  stat_summary(fun.data = give.n, geom = "text", hjust = 0.5, position = position_dodge(0.6)) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
myboxplot

pdf("./2_output/boxplot_morethan31genes.pdf", width = 12, height = 5)
myboxplot
dev.off()




# check normality

shapiro.test(trustworthy_dataset$Mean_age[0:5000])
ad.test(trustworthy_dataset$Mean_age)

normality_data <- trustworthy_dataset %>%
  select(!Mean_age) %>%
  summarise_all(.funs = funs(statistic = shapiro.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))

# long format
data_long <- tidyr::gather(trustworthy_dataset, variable, value, 1:19)
data_long <- drop_na(data_long)

# one-way ANOVA
one.way <- aov(value ~ variable, data = data_long)
summary(one.way)

# Tukey's Honestly-Significant Difference
TukeyHSD(one.way)

# START HERE

library(tidyverse)
library(nortest)

setwd("/media/nikos/LaCie/Nikos_data/GEA/GEA_Yann_noHe/2_output/FDR_threshold/Gene_age_estimates/")
all_genes <- read.table("./2_output/all_genes_and_variables.tbl", sep = ",", header = T)


ggplot(stack(all_genes[,c(4,6:length(all_genes))]), aes(x = ind, y = values)) +
  geom_boxplot()

biovars <- select(all_genes, c(4,contains("bio")))
ggplot(stack(biovars), aes(x = ind, y = values)) +
  geom_boxplot()

monthvars <- select(all_genes, c(4,contains(c("Feb", "June", "elevation"))))
ggplot(stack(monthvars), aes(x = ind, y = values)) +
  geom_boxplot()

# get number of genes per variable
colSums(!is.na(all_genes))

# only keep variables with more than 31 genes
trustworthy_dataset <- select(all_genes, c(4,7:39))
trustworthy_dataset <- trustworthy_dataset[,colSums(!is.na(trustworthy_dataset)) > 31]

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

# boxplots only with trustworthy dataset
myboxplot <- ggplot(stack(trustworthy_dataset), aes(x = ind, y = values)) +
  geom_boxplot() +
  scale_y_continuous(breaks = round(seq(0, max(trustworthy_dataset$Mean_age, na.rm = T), by = 10000),1)) +
  xlab("") + ylab("gene age in years") +
  stat_summary(fun.data = give.n, geom = "text", hjust = 0.5, position = position_dodge(0.6)) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
myboxplot

pdf("./2_output/boxplot_morethan31genes.pdf", width = 12, height = 5)
myboxplot
dev.off()




# check normality

shapiro.test(trustworthy_dataset$Mean_age[0:5000])
ad.test(trustworthy_dataset$Mean_age)

normality_data <- trustworthy_dataset %>%
  select(!Mean_age) %>%
  summarise_all(.funs = funs(statistic = shapiro.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))

# long format
data_long <- tidyr::gather(trustworthy_dataset, variable, value, 1:19)
data_long <- drop_na(data_long)

# one-way ANOVA
one.way <- aov(value ~ variable, data = data_long)
summary(one.way)

# Tukey's Honestly-Significant Difference
TukeyHSD(one.way)

# FILTER for <5 SNPs per gene and REPEAT

library(tidyverse)
library(nortest)
library(datarium)
library(rstatix)

setwd("/media/nikos/LaCie/Nikos_data/GEA/GEA_Yann_noHe/2_output/FDR_threshold/Gene_age_estimates/")

all_genes <- read.table("./2_output/all_genes_and_variables.tbl", sep = ",", header = T)

all_genes <- all_genes %>%
  filter(Number_of_SNPs > 4)

ggplot(stack(all_genes[,c(4,6:length(all_genes))]), aes(x = ind, y = values)) +
  geom_boxplot()

biovars <- select(all_genes, c(4,contains("bio")))
ggplot(stack(biovars), aes(x = ind, y = values)) +
  geom_boxplot()

monthvars <- select(all_genes, c(4,contains(c("Feb", "June", "elevation"))))
ggplot(stack(monthvars), aes(x = ind, y = values)) +
  geom_boxplot()

# get number of genes per variable
colSums(!is.na(all_genes))

# only keep variables with more than 14 genes
trustworthy_dataset <- select(all_genes, c(4,7:39))
trustworthy_dataset <- trustworthy_dataset[,colSums(!is.na(trustworthy_dataset)) > 14]

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

# boxplots only with trustworthy dataset
myboxplot <- ggplot(stack(trustworthy_dataset), aes(x = reorder(ind, values, na.rm = TRUE, median), y = values)) +
  geom_boxplot() +
  scale_y_continuous(breaks = round(seq(0, max(trustworthy_dataset$Mean_age, na.rm = T), by = 10000),1)) +
  xlab("") + ylab("gene age in years") +
  stat_summary(fun.data = give.n, geom = "text", hjust = 0.5, position = position_dodge(0.6)) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
myboxplot

pdf("./2_output/boxplot_morethan14genes_morethan4SNPspergene_withXtX_ordered.pdf", width = 12, height = 5)
myboxplot
dev.off()

# check normality

# Build the linear model
library(ggpubr)
model <- lm(values ~ ind, data = stack(trustworthy_dataset))
ggqqplot(residuals(model))
shapiro_test(residuals(model))

shapiro.test(trustworthy_dataset$Mean_age[0:5000])
ad.test(trustworthy_dataset$Mean_age)

normality_data <- trustworthy_dataset %>%
  select(!Mean_age) %>%
  summarise_all(.funs = funs(statistic = shapiro.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))

# long format
data_long <- tidyr::gather(trustworthy_dataset, variable, value, 1:length(trustworthy_dataset))
data_long <- drop_na(data_long)

# data not normal. Kruskal-Wallis
library(FSA)
library(rcompanion)
library(robustbase)
my_test <- kruskal.test(value ~ variable, data = data_long)
my_test

# Dunn test for multiple comparisons (Benjamini-Hochberg correction for the p-values of multiple comparisons)
PT <- dunnTest(value ~ variable, data = data_long, method = "bh")
PT

colMeans(trustworthy_dataset, na.rm = T)
colMedians(as.matrix(trustworthy_dataset), na.rm = T)

#load DescTools library
library(DescTools)

#perform Dunnett's Test
DunnettTest(x = data_long$value, g = data_long$variable, control = "XtX")

# Wilcoxon
library(stats)

p1 <- wilcox.test(trustworthy_dataset$`aridity_March-June`, trustworthy_dataset$Mean_age)$p.value
p2 <- wilcox.test(trustworthy_dataset$`aridity_Nov-Feb`, trustworthy_dataset$Mean_age)$p.value
p3 <- wilcox.test(trustworthy_dataset$bio_1, trustworthy_dataset$Mean_age)$p.value
p4 <- wilcox.test(trustworthy_dataset$bio_10, trustworthy_dataset$Mean_age)$p.value
p5 <- wilcox.test(trustworthy_dataset$bio_11, trustworthy_dataset$Mean_age)$p.value
p6 <- wilcox.test(trustworthy_dataset$bio_12, trustworthy_dataset$Mean_age)$p.value
p7 <- wilcox.test(trustworthy_dataset$bio_13, trustworthy_dataset$Mean_age)$p.value
p8 <- wilcox.test(trustworthy_dataset$bio_14, trustworthy_dataset$Mean_age)$p.value
#p9 <- wilcox.test(trustworthy_dataset$bio_15, trustworthy_dataset$Mean_age)$p.value
p10 <- wilcox.test(trustworthy_dataset$bio_16, trustworthy_dataset$Mean_age)$p.value
#p11 <- wilcox.test(trustworthy_dataset$bio_17, trustworthy_dataset$Mean_age)$p.value
p12 <- wilcox.test(trustworthy_dataset$bio_18, trustworthy_dataset$Mean_age)$p.value
#p13 <- wilcox.test(trustworthy_dataset$bio_2, trustworthy_dataset$Mean_age)$p.value
#p14 <- wilcox.test(trustworthy_dataset$bio_4, trustworthy_dataset$Mean_age)$p.value
p15 <- wilcox.test(trustworthy_dataset$bio_5, trustworthy_dataset$Mean_age)$p.value
p16 <- wilcox.test(trustworthy_dataset$bio_9, trustworthy_dataset$Mean_age)$p.value
p17 <- wilcox.test(trustworthy_dataset$elevation, trustworthy_dataset$Mean_age)$p.value
p18 <- wilcox.test(trustworthy_dataset$`prec_March-June`, trustworthy_dataset$Mean_age)$p.value
p19 <- wilcox.test(trustworthy_dataset$`prec_Nov-Feb`, trustworthy_dataset$Mean_age)$p.value
p20 <- wilcox.test(trustworthy_dataset$`srad_March-June`, trustworthy_dataset$Mean_age)$p.value
#p21 <- wilcox.test(trustworthy_dataset$`srad_Nov-Feb`, trustworthy_dataset$Mean_age)$p.value
p22 <- wilcox.test(trustworthy_dataset$`tavg_March-June`, trustworthy_dataset$Mean_age)$p.value
#p23 <- wilcox.test(trustworthy_dataset$`tavg_Nov-Feb`, trustworthy_dataset$Mean_age)$p.value
p24 <- wilcox.test(trustworthy_dataset$`tmax_March-June`, trustworthy_dataset$Mean_age)$p.value
#p25 <- wilcox.test(trustworthy_dataset$`tmax_Nov-Feb`, trustworthy_dataset$Mean_age)$p.value
p26 <- wilcox.test(trustworthy_dataset$`tmin_March-June`, trustworthy_dataset$Mean_age)$p.value
p27 <- wilcox.test(trustworthy_dataset$XtX, trustworthy_dataset$Mean_age)$p.value


wilcox_bonferroni_results <- as.data.frame(p.adjust(c(p1,p2,p3,p4,p5,p6,p7,p8,p10,p12,p15,p16,p17,p18,p19,
                                                      p20,p22,p24,p26,p27), method = "bonferroni"))
rownames(wilcox_bonferroni_results) <- colnames(trustworthy_dataset)[2:length(trustworthy_dataset)]
write.table(wilcox_bonferroni_results, "./2_output/wilcox_bonferroni_results.tbl", row.names = T, col.names = T, quote = F, sep = ",")









