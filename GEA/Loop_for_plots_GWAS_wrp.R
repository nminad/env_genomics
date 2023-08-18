library(hdi)
library(tidyverse)
library(cowplot)

setwd("../2_output/gemma_output/")

# Open and read the file list including all the variable names
f <- list.files(full.names = T, pattern = "assoc")
f_short <- list.files(full.names = F, pattern = "assoc")

# Start the loop
for (i in 1:length(f)) {

  gwscan1 <- read.table(f[i], as.is = "rs", header = TRUE)

  # exclude rows with NaN
  gwscan <- gwscan1[complete.cases(gwscan1), ]

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
  
  n <- length(gwscan$p_lrt)
  
  # Compute the negative log10(p-values), and sort them from largest to smallest.
  y <- rev(sort(-log10(gwscan$p_lrt)))
  
  # Create the q-q plot.
  QQ <- ggplot(data.frame(x = -log10((1:n)/n),y = y),aes(x = x,y = y)) +
    geom_abline(intercept = 0,slope = 1,color = "magenta") +
    geom_point(color = "dodgerblue",shape = 20,size = 2) +
    labs(x = "Expected -log10 p-value",
         y = "Observed -log10 p-value") +
    geom_hline(yintercept = BC) +
    geom_hline(yintercept = fdrlog10, linetype = "dotted") +
    ggtitle(f_short[i]) +
    theme(axis.line = element_blank())
  
  
  
  n <- nrow(gwscan)
  gwscan <- cbind(gwscan,marker = 1:n)
  
  # Convert the p-values to the -log10 scale
  gwscan <- transform(gwscan,p_lrt = -log10(p_lrt))
  
  # Add column "odd.chr" to the table, and find the positions of the chromosomes along the x-axis
  gwscan <- transform(gwscan,odd.chr = (chr %% 2) == 1)
  x.chr  <- tapply(gwscan$marker,gwscan$chr,mean)
  
  # Create the genome-wide scan
  MH <- ggplot(gwscan,aes(x = marker,y = p_lrt,color = odd.chr)) +
    geom_point(size = 1.5, shape = 20) +
    scale_x_continuous(breaks = x.chr,labels = 1:5) +
    scale_y_continuous(limits = c(0, max(gwscan$p_lrt) + 1), expand = c(0,0)) +
    scale_color_manual(values = c("dodgerblue3","azure4"),guide = "none") +
    labs(x = "",y = "-log10 p-value") +
    geom_hline(yintercept = BC) +
    geom_hline(yintercept = fdrlog10, linetype = "dotted") +
    theme_cowplot() +
    ggtitle(f_short[i]) +
    theme(axis.line.y = element_line(lineend = "butt"),
          axis.line = element_blank(),
          axis.ticks.x = element_line(),
          plot.title = element_text(hjust = 0.5))
  #pdf(paste("/2_output/QQ_MH_plots/QQ_plot_",f_short[i],".pdf",sep = ""), width = 12, height = 4)
  #print(QQ)
  #dev.off()
  #pdf(paste("/2_output/QQ_MH_plots/MH_plot_",f_short[i],".pdf",sep = ""), width = 12, height = 4)
  #print(MH)
  #dev.off()
  
  png(paste("/2_output/QQ_MH_plots/QQ_plot_",f_short[i],".png",sep = ""), width = 12, height = 4, units = "in", res = 150)
  print(QQ)
  dev.off()
  png(paste("/2_output/QQ_MH_plots/MH_plot_",f_short[i],".png",sep = ""), width = 12, height = 4, units = "in", res = 150)
  print(MH)
  dev.off()
}

# for png: The units in which height and width are given. Can be px (pixels, the default), in (inches), cm or mm.

