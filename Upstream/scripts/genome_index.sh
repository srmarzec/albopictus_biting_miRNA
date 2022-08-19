#!/bin/bash
#SBATCH --job-name=genome_index --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=[]@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script runs makes a genome index with bowtie #
#-----------------------------------------------------------------------------#

module load bowtie/1.3.1

#- Set variables ----------------------------------------------------------------#

ref_genome=/home/sm3679/albopictus_biting_miRNA/albopictus_genome/aedes_albopictus_AalbF3.fa
index_path_prefix=/home/sm3679/albopictus_biting_miRNA/albopictus_genome/genome_index/Aal


#- RUN bowtie ----------------------------------------------------------------#

bowtie-build ${ref_genome} ${index_path_prefix}


#- Unload module----------------#
module unload bowtie/1.3.1

#- FIN -----------------------------------------------------------------------#
