# This script runs DESeq2 (specifically starting with a count matrix file)

##Install DESeq2
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install('EnhancedVolcano')



#Load Libraries
library(DESeq2); library(ggplot2); library(EnhancedVolcano); library(tidyverse)


# Set working directory to source file location

# Read in the count data matrix
countData <- read.csv('../data/counts_miRNA_DESeq_set.csv', header = TRUE)
head(countData)

# Removing the NB9 column because this is an obvious outlier with its lack of reads
countData <- countData %>%
  select(-NB9)

#Create the metaData table (this could alternatively be made externally and read in)
sampleNames <- colnames(countData)[2:8] #Here I am calling the column names of the Count matrix, minus the first column which is the name for the miRNA annotations
sampleConditions <- substr(sampleNames, 1, 2) #to get conditions I'm pulling the second letter from the sample names, which is either BI (biting) or NB (non-biting)

#metaData <- data.frame(sampleName = sampleNames,
#                          condition = sampleConditions)
metaData <- read.csv("../data/miRNA_metadata.csv", header = T)
str(metaData)
metaData$condition <- as.factor(metaData$condition)
metaData$cage <- as.factor(metaData$cage)
metaData$day <- as.factor(metaData$day)
metaData <- metaData %>%
  filter(sampleName != "NB9")

# Make the DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = metaData, 
                              design= ~ condition, tidy = TRUE)
dds

#DESeq recommends a pre-filtering step to reduce memory size and increase speed. They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Relevel the condition for what to compare to. Likely not important in our case with two conditions, but you would want everything compared to the control. (Default is first condition alphabetically)
dds$condition <- relevel(dds$condition, ref = "NB")

# A bit of quality control
# Look at the distribution of the count data across the samples
librarySizes <- colSums(counts(dds))

barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, ylim = c(0,1.6e+07),
        main="Barplot of library sizes")

logcounts <- log2(counts(dds) + 1)
head(logcounts)

#Is there any difference between per gene counts for each of the sample groups?
statusCol <- as.numeric(factor(dds$condition)) + 1  # make a colour vector

boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=statusCol)

#Adding median log counts
abline(h=median(as.matrix(logcounts)), col="blue")

# Transform normalized counts using the rlog function
rld <- rlog(dds, blind=TRUE)

# PCA plot with other PCs besides 1 & 2
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

pca_df <- cbind(metaData, pca$x)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
ggplot(pca_df) + 
  geom_point(aes(x=PC6, y=PC7, color = condition), size = 3) +
  xlab(paste0("PC6: ",round(percentVar[6] * 100),"% variance")) +
  ylab(paste0("PC7: ",round(percentVar[7] * 100),"% variance"))

# look at loadings of PCA
head(pca$rotation)

# Differential Expression Analysis
#Run DESeq
dds <- DESeq(dds)
res <- results(dds)
res

# Apparently, it is common to shrink the log fold change estimate for better visualization and ranking of genes. Often, lowly expressed genes tend to have relatively high levels of variability so shrinking can reduce this.
resultsNames(dds)
res_LFC <- lfcShrink(dds, coef="condition_BI_vs_NB", type="apeglm")
res_LFC



# Volcano plot of the lfcShrink data. For various setting try: `browseVignettes("EnhancedVolcano")`
# Note for selecting fold change cutoff: log2foldchange 0.58 is equal to a 1.5 fold change
EnhancedVolcano(res_LFC,
                lab = rownames(res_LFC),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = NA,
                title = NULL,
                subtitle = NULL,
                xlim = c(-3.5, 3.5),
                ylim = c(0,15),
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 5.0,
                caption = NULL)


# Full result print out
full_res <- res_LFC %>%
  data.frame() %>%
  rownames_to_column(var="miRNA") %>% 
  as.data.frame()

# Write out a table of these significant differentially expressed genes
write.csv(dplyr::select(full_res, miRNA, log2FoldChange, padj), 
          file="../output/BvNB_LFCshrink_padj_FULL.csv", row.names = F)
