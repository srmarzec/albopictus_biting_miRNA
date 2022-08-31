#!/bin/bash
#SBATCH --job-name=quantifier_FULL --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sm3679@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script runs the quantifier module of miRDeep2 #
#-----------------------------------------------------------------------------#

#- Go into Conda environment where quantifier.pl is  --------------------------------#

source activate conda-env

cd /home/sm3679/albopictus_biting_miRNA/quantified_miRNA_dir/

#- Set variables ----------------------------------------------------------------#

reads=/home/sm3679/albopictus_biting_miRNA/miRDeep_dir/processed_reads.fa
mature_miRNA=/home/sm3679/albopictus_biting_miRNA/known_miRNAs/FULL_Mature_noSpace.fa
precursor_miRNA=/home/sm3679/albopictus_biting_miRNA/known_miRNAs/FULL_Precursors_noSpace.fa
config_file=/home/sm3679/albopictus_biting_miRNA/miRDeep_dir/config_file.txt


#- RUN quantifier.pl ----------------------------------------------------------------#

quantifier.pl -p ${precursor_miRNA} -m ${mature_miRNA} -r ${reads} -c ${config_file} -d -k

#- FIN -----------------------------------------------------------------------#
