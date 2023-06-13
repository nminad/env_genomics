#!/bin/bash
#SBATCH --mem=25G
#SBATCH -p serial
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=8
# Walltime format hh:mm:ss
#SBATCH --time=36:00:00
# Output and error files
#SBATCH -o filter.%J.out
#SBATCH -e filter.%J.err
# **** Put all #SBATCH directives above this line! ****
# **** Otherwise they will not be in effective! ****
#
# **** Actual commands start here ****
# Load modules here (safety measure)
module load gencore/1
module load gencore_variant_detection/1.0
module load gcc/4.9.3
source activate ANGSD
#./relate_v1.1.7_x86_64_static/bin/RelateFileFormats --mode GenerateSNPAnnotations --haps haploid_Bd1.haps --sample Bd1_haploid.sample --poplabels samples.poplabels -o Bd1_annotated --mut Output_Bd1.mut 
#rm Bd1_annotated.annot
#mv Bd1_annotated.mut Output_Bd1.mut 

#./relate_v1.1.7_x86_64_static/bin/RelateFileFormats --mode GenerateSNPAnnotations --haps haploid_Bd2.haps --sample Bd2_haploid.sample --poplabels samples.poplabels -o Bd2_annotated --mut Output_Bd2.mut 
#rm Bd2_annotated.annot
#mv Bd2_annotated.mut Output_Bd2.mut 

#./relate_v1.1.7_x86_64_static/bin/RelateFileFormats --mode GenerateSNPAnnotations --haps haploid_Bd3.haps --sample Bd3_haploid.sample --poplabels samples.poplabels -o Bd3_annotated --mut Output_Bd3.mut 
#rm Bd3_annotated.annot
#mv Bd3_annotated.mut Output_Bd3.mut 

#./relate_v1.1.7_x86_64_static/bin/RelateFileFormats --mode GenerateSNPAnnotations --haps haploid_Bd5.haps --sample Bd5_haploid.sample --poplabels samples.poplabels -o Bd5_annotated --mut Output_Bd5.mut 
#rm Bd5_annotated.annot
#mv Bd5_annotated.mut Output_Bd5.mut 

#rename Output_Bd Output_chr Output_Bd*
../relate_v1.1.7_x86_64_static/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i Output --first_chr 1 --last_chr 5 -m 7e-9 --pop_of_interest A_Balkan,A_East,A_Italy,B_East,B_West --poplabels ./samples.poplabels --threads 8 --years_per_gen 1 --threshold 0 --seed 1234 -o Eff_pop_size

