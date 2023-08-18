#!/bin/bash

while read line; do
	BIOCLIM=$(echo $line | cut -d',' -f1)
	bedtools intersect -a ../2_output/FDR/Bd_chr/${BIOCLIM}.FDR.sig.bed -b ~/Documents/Project/Bdistachyon/v3.2/annotation/Bdistachyon_556_v3.2.gene.gff3 -wb | tee ../2_output/FDR/bedtools_intersect/${BIOCLIM}_bedtools_intersect_FDR_wb_out.txt
done < ../0_data/file_list.txt

# -wb Reporting the original B feature
# Similarly, one can force bedtools intersect to report the original “B” feature when an overlap is found. If just -wb is used, the overlapping portion of A will be reported followed # by the original “B”.
