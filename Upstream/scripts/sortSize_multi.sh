#!/bin/bash
#SBATCH --job-name=sortSize_multi --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=[]@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script runs the python script trimANDsizeSort.py on multiple files #
#-----------------------------------------------------------------------------#

module load python3/3.6.8

#- Set variables ----------------------------------------------------------------#

input_dir=/home/sm3679/albopictus_biting_miRNA/filtered_dir


#- RUN python script ----------------------------------------------------------------#

files=(${input_dir}/*_flt.fastq)
for file in ${files[@]}
do
python /home/sm3679/albopictus_biting_miRNA/scripts/python_scripts/trimANDsizeSort.py ${file} 
done

#- FIN -----------------------------------------------------------------------#
