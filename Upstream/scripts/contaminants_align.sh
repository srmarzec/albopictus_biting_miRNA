#!/bin/bash
#SBATCH --job-name=contaminants_align --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sm3679@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This runs bowtie2 for contaminant mapping#
#---------------------------------------------------------------------#
 
#- Load Module--------------------------------------------#
module load bowtie2/2.4.4
 
#- DEFINE FILE LOCATIONS--------------------------------------------#
 
index_path=/home/sm3679/albopictus_biting_miRNA/albopictus_genome/contaminant_index/contaminant
input_dir=/home/sm3679/albopictus_biting_miRNA/trim_dir
filtered_dir=/home/sm3679/albopictus_biting_miRNA/filtered_dir
sam_dir=/home/sm3679/albopictus_biting_miRNA/sam_dir
 
#- RUN command ----------------#
files=(${input_dir}/*_cln.fastq.gz)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _cln.fastq.gz`

bowtie2 -x ${index_path} -U ${input_dir}/${base}_cln.fastq.gz -S ${sam_dir}/${base}_flt.cln.sam --un-gz ${filtered_dir}/${base}_flt.fastq.gz

done 

#- Unload module----------------#
module unload bowtie2/2.4.4

#- FIN -----------------------#
