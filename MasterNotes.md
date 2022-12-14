# Workflow details for the *Aedes albopictus* biting miRNA analysis

## Upstream 
Information about sequence reads (raw read count, read alignment rate, etc.) for the data can be found [here](https://docs.google.com/spreadsheets/d/1F3cfOhkYX_w3hRzcgdPXem40MCuqEjhe6jSWc0WScSM/edit?usp=sharing)

### Data Accession
Data was generated by collaborator...

The raw reads are available in NCBI’s short read archive (SRA) under accession number XXXXYYYY

### Preprocessing and Quality Control

Trimmomatic (version 0.39) was used to trim sequence reads based on quality ([script](https://github.com/srmarzec/albopictus_biting_miRNA/blob/main/Upstream/scripts/trim.sh))

We removed specific adapter sequences for the Illumina small RNA 3' adapter, which can be found [here](https://github.com/srmarzec/albopictus_biting_miRNA/blob/main/Upstream/smRNA_NexFlex_adapters.fa) after contacting the sequencing center for the specific adapters used. 

FastQC (v0.11.9) was used for quality control visualization ([script](https://github.com/srmarzec/albopictus_biting_miRNA/blob/main/Upstream/scripts/fastqc.sh))

Preliminary fastqc showed a poor "per sequence base content" for the first few bases. Therefore, we used headcrop four in the begginning, but there is no command to crop the four at the end which could be any base pairs. Note that SE settings were used. 

From the fastqc files, you can see that the per base sequences quality improved and the adapter content was removed. There are still flags for the following with either warnings or fails: Per base sequence content, Per sequence GC content, Sequence Length Distribution, Sequence Duplication Levels, Overrepresented sequences, Per Tile Sequence Quality. It is under a general consensus, however, that these flags will not significantly affect our analysis and that there might be a biological reason behind them. 

#### Cleaning out other small RNAs
Aim: to remove tRNA and other contaminates.

All tRNAs and rRNAs features from the gff file for AalbF3 were collected and we made a subset fasta with those sequences as a "contaminants" file.

We indexed this contaminants file using Bowtie2 (v2.4.4) with this [script](https://github.com/srmarzec/albopictus_biting_miRNA/blob/main/Upstream/scripts/contaminants_index.sh).

Ran alignment and put non-aligned reads into filtered fasta files using this [script](https://github.com/srmarzec/albopictus_biting_miRNA/blob/main/Upstream/scripts/contaminants_align.sh).

The output files contain all the sequences that did NOT align with our contanminants file i.e. files with sequences that were not tRNAs or rRNAs (presumably mostly miRNAs are left).

#### Size sorting

We wanted to only keep reads that fell within 18 - 24 bases which is what we consider the size of miRNAs. This was done with a custom [python script](https://github.com/srmarzec/albopictus_biting_miRNA/blob/main/Upstream/scripts/python_scripts/trimANDsizeSort.py) for all the files with this [script](https://github.com/srmarzec/albopictus_biting_miRNA/blob/main/Upstream/scripts/sortSize_multi.sh).

### miRDeep2 
#### Index with Bowtie

Mapping was done using the *Aedes albopictus* reference genome (GCA_018104305.1) found on NCBI

The genome was indexed for miRDeep2 using Bowtie (v1.3.1) with this [script](https://github.com/srmarzec/albopictus_biting_miRNA/blob/main/Upstream/scripts/genome_index.sh). Note, for miRDeep2 to run properly, there can be no whitespace in the headers of the reference genome fasta so I had to remove/replace the whitespaces before indexing.

#### Run miRDeep2 to count reads for each miRNA

miRDeep2 (miRDeep2.0.1.3) was installed and run in a conda virtual environment which was created following this [markdown](https://github.com/srmarzec/Culex_Biting_miRNA/blob/main/misc/Conda_VirtualEnvironment.md).

miRDeep2 was [run](https://github.com/srmarzec/albopictus_biting_miRNA/blob/main/Upstream/scripts/miRDeep2.sh) and we looked for novel miRNAs identified in our samples and added these to the input list of known miRNAs for sake of quantification (mapping reads to miRNAs).

We then mapped reads to the known miRNAs with this [script](https://github.com/srmarzec/albopictus_biting_miRNA/blob/main/Upstream/scripts/quantifier_FULL.sh). This produced a count matrix we used for downstream analysis. 

## Downstream

All downstream analysis done in R (v4.0.2)

### DESeq
Using DESeq2 (v1.30.1) ([script](https://github.com/srmarzec/albopictus_biting_miRNA/blob/main/Downstream/DESeq.R))

As a result of performing differential expression analysis on our miRNAs, we obtained no differentially abundant miRNAs.
