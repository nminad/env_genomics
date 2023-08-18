#!/bin/bash

while read line; do
	NUM=$(echo $line | cut -d',' -f2)
	CODENAME=$(echo $line | cut -d',' -f3)
	/home/nikos/software/gemma-0.98.5-linux-static-AMD64 -bfile ../0_data/binary \
	-k ../0_data/output/structure.cXX.txt -lmm 4 -n $NUM -maf 0.05 -o ULMM_${CODENAME}_output
done < phenotype_names.csv



