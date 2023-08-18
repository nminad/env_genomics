
## Calculate pi using pixy


## for the genes of interest 
pixy --stats pi --vcf /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis_332_Bd1_pixy.vcf.gz --populations /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis332_clades_popfile.txt --output_folder /media/HDD2/Robert/Pixy --output_prefix Bdis_332_Bd1_genes_pi_pixy_output --n_cores 1 --bed_file Bdistachyon_314_v3.2.gene.Bd1.bed

pixy --stats pi --vcf /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis_332_Bd2_pixy.vcf.gz --populations /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis332_clades_popfile.txt --output_folder /media/HDD2/Robert/Pixy --output_prefix Bdis_332_Bd2_genes_pi_pixy_output --n_cores 1 --bed_file Bdistachyon_314_v3.2.gene.Bd2.bed

pixy --stats pi --vcf /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis_332_Bd3_pixy.vcf.gz --populations /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis332_clades_popfile.txt --output_folder /media/HDD2/Robert/Pixy --output_prefix Bdis_332_Bd3_genes_pi_pixy_output --n_cores 1 --bed_file Bdistachyon_314_v3.2.gene.Bd3.bed

pixy --stats pi --vcf /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis_332_Bd4_pixy.vcf.gz --populations /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis332_clades_popfile.txt --output_folder /media/HDD2/Robert/Pixy --output_prefix Bdis_332_Bd4_genes_pi_pixy_output --n_cores 1 --bed_file Bdistachyon_314_v3.2.gene.Bd4.bed

pixy --stats pi --vcf /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis_332_Bd5_pixy.vcf.gz --populations /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis332_clades_popfile.txt --output_folder /media/HDD2/Robert/Pixy --output_prefix Bdis_332_Bd5_genes_pi_pixy_output --n_cores 1 --bed_file Bdistachyon_314_v3.2.gene.Bd5.bed


## genome wide pi:
pixy --stats pi --vcf /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis_332_Bd1_pixy.vcf.gz --populations /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis332_clades_popfile.txt --output_folder /media/HDD2/Robert/Pixy --output_prefix Bdis_332_Bd1_genome_wide_5000_windows_pi_pixy_output --n_cores 1 --window_size 5000

pixy --stats pi --vcf /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis_332_Bd2_pixy.vcf.gz --populations /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis332_clades_popfile.txt --output_folder /media/HDD2/Robert/Pixy --output_prefix Bdis_332_Bd2_genome_wide_5000_windows_pi_pixy_output --n_cores 1 --window_size 5000

pixy --stats pi --vcf /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis_332_Bd3_pixy.vcf.gz --populations /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis332_clades_popfile.txt --output_folder /media/HDD2/Robert/Pixy --output_prefix Bdis_332_Bd3_genome_wide_5000_windows_pi_pixy_output --n_cores 1 --window_size 5000

pixy --stats pi --vcf /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis_332_Bd4_pixy.vcf.gz --populations /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis332_clades_popfile.txt --output_folder /media/HDD2/Robert/Pixy --output_prefix Bdis_332_Bd4_genome_wide_5000_windows_pi_pixy_output --n_cores 1 --window_size 5000

pixy --stats pi --vcf /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis_332_Bd5_pixy.vcf.gz --populations /media/HDD2/Lars_temp/WholeGenomeGVCF/pixy/Bdis332_clades_popfile.txt --output_folder /media/HDD2/Robert/Pixy --output_prefix Bdis_332_Bd5_genome_wide_5000_windows_pi_pixy_output --n_cores 1 --window_size 5000




