

## Run bpp

## make a file with the path to the bam files: List_of_paths_to_bam.txt

## index the bam files that were not previously indexed:
samtools index -b SRR4162910_sort.bam
samtools index -b SRR4162906_sort.bam
samtools index -b SRR4162905_sort.bam
samtools index -b SRR4162891_sort.bam
samtools index -b D13_FDSW202357929-1r_H377MDSXY_L4_sort.bam
samtools index -b D60_merged.bam
samtools index -b D16_FDSW202357932-1r_H377MDSXY_L4_sort.bam
samtools index -b D69_FDSW202357985-1r_H375VDSXY_L3_sort.bam


## run bpp
./bpp_from_Christoph/scripts/consensus_from_bam.py -L List_of_paths_to_bam.txt -a bpp_from_Christoph/w1000_s100000.200intervals.sorted.bed

## re-naming fastas
mv Cm4.consensi_di_bam.fasta Cm4.fasta
mv Ren22.consensi_di_bam.fasta Ren22.fasta
mv Msa27.consensi_di_bam.fasta Msa27.fasta
mv Lb13.consensi_di_bam.fasta Lb13.fasta
mv San12.consensi_di_bam.fasta San12.fasta
mv Cb23.consensi_di_bam.fasta Cb23.fasta
mv Mca12.consensi_di_bam.fasta Mca12.fasta
mv Tso18.consensi_di_bam.fasta Tso18.fasta
mv 1a_32_12.consensi_di_bam.fasta 1a_32_12.fasta
mv 1c_25_14.consensi_di_bam.fasta 1c_25_14.fasta
mv 2_14_13.consensi_di_bam.fasta 2_14_13.fasta
mv 4_52_6.consensi_di_bam.fasta 4_52_6.fasta
mv BdTR2B.consensi_di_bam.fasta BdTR2B.fasta
mv Bd30-1.consensi_di_bam.fasta Bd30-1.fasta
mv BdTR3C.consensi_di_bam.fasta BdTR3C.fasta
mv BdTR7a.consensi_di_bam.fasta BdTR7a.fasta
mv Luc1.consensi_di_bam.fasta Luc1.fasta
mv Bd21.consensi_di_bam.fasta Bd21.fasta
mv Cro24.consensi_di_bam.fasta Cro24.fasta
mv ABR6.consensi_di_bam.fasta ABR6.fasta
mv D69_FDSW202357985-1r_H375VDSXY_L3_sort.consensi_di_bam.fasta D69.fasta
mv SRR4162910_sort.consensi_di_bam.fasta SRR4162910.fasta
mv SRR4162906_sort.consensi_di_bam.fasta SRR4162906.fasta
mv SRR4162905_sort.consensi_di_bam.fasta SRR4162905.fasta
mv SRR4162891_sort.consensi_di_bam.fasta SRR4162891.fasta
mv D13_FDSW202357929-1r_H377MDSXY_L4_sort.consensi_di_bam.fasta D13.fasta
mv D60_merged.consensi_di_bam.fasta D60.fasta
mv D16_FDSW202357932-1r_H377MDSXY_L4_sort.consensi_di_bam.fasta D16.fasta


find /media/HDD2/Robert/bpp -type f -name "*.fasta" > paths_to_fastas.txt
/media/HDD2/Robert/bpp/bpp_from_Christoph/scripts/fasta2supermatrix.py paths_to_fastas.txt phylyp

./fasta2supermatrix_2.py paths_to_fastas.txt phylyp

## make a supermatrix file (supermatrix_new_b_2.ph)

## run bpp
/media/HDD2/Robert/bpp/bpp-4.6.2-linux-x86_64/bin/bpp --cfile A00_re_run.bpp.ctl

mv FigTree.tre FigTree_run1.tre

## rum it a secound time
/media/HDD2/Robert/bpp/bpp-4.6.2-linux-x86_64/bin/bpp --cfile A00_re_run.bpp_2.ctl

mv FigTree.tre FigTree_run2.tre









## get bpp results
## In R


library(bppr)
library(bayesAB)
library(ape)
library(coda)


mcmc <- read.table('/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/bpp_analysis/Results/run1.mcmc', header=T)



myTree <- ape::read.tree(text='(((A_East: 0.000090, A_Italia: 0.000090)A_EastA_Italia[&height_95%_HPD={0.00004500, 0.00012900}, theta=0.0023758]: 0.000221, (B_East: 0.000160, B_West: 0.000160)B_EastB_West[&height_95%_HPD={0.00012300, 0.00019700}, theta=0.0012283]: 0.000151)A_EastA_ItaliaB_EastB_West[&height_95%_HPD={0.00023200, 0.00039900}, theta=0.0093400]: 0.000376, C_Italia: 0.000687)A_EastA_ItaliaB_EastB_WestC_Italia[&height_95%_HPD={0.00057200, 0.00081000}, theta=0.0242490];')

# Calibrate time
mcmc_cal <- msc2time.r(mcmc, u.mean = 7e-9, u.sd = 2.213594e-09, g.mean = 1, g.sd = 0.1)
mcmc_s <- mcmc.summary(mcmc_cal, prob = 0.95)

class(mcmc_cal)
write.table(mcmc_cal, file = "/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/bpp_analysis/Results/run1.mcmc_cal")


# Plot tree
pdf(file = "/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/bpp_analysis/Results/FigTree_run1_plot.pdf")      
mcmc2densitree(myTree, mcmc_cal, time.name="t_", thin = 0.01, alpha=0.01, col="blue")
title(xlab="Divergence time (years)")
dev.off()


mcmc <- read.table('/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/bpp_analysis/Results/run2.mcmc', header=T)

myTree <- ape::read.tree(text='(((A_East: 0.000090, A_Italia: 0.000090)A_EastA_Italia[&height_95%_HPD={0.00004700, 0.00013000}, theta=0.0023974]: 0.000224, (B_East: 0.000159, B_West: 0.000159)B_EastB_West[&height_95%_HPD={0.00012100, 0.00019700}, theta=0.0012538]: 0.000154)A_EastA_ItaliaB_EastB_West[&height_95%_HPD={0.00022800, 0.00040100}, theta=0.0092063]: 0.000370, C_Italia: 0.000683)A_EastA_ItaliaB_EastB_WestC_Italia[&height_95%_HPD={0.00057200, 0.00079000}, theta=0.0242791];')

# Calibrate time
mcmc_cal <- msc2time.r(mcmc, u.mean = 7e-9, u.sd = 2.213594e-09, g.mean = 1, g.sd = 0.1)
mcmc_s <- mcmc.summary(mcmc_cal, prob = 0.95)

write.table(mcmc_cal, file = "/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/bpp_analysis/Results/run2.mcmc_cal")


# Plot tree
pdf(file = "/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/bpp_analysis/Results/FigTree_run2_plot.pdf")      
mcmc2densitree(myTree, mcmc_cal, time.name="t_", thin = 0.01, alpha=0.01, col="blue")
title(xlab="Divergence time (years)")
dev.off()







