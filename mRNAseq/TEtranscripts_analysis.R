####################################################################################################
# Functions for the analysis of RNASeq data (or data processed with DESeq2)
# source("rnaSeqAnalysisFuncs.R") to get the functions in the environment
# Adapted by Miguel Almeida from a script by Jon Price. 
####################################################################################################

#######################################
#Load libraries  ---- 
#######################################

library(ggplot2)
library(pheatmap)
library(reshape2)
library(genefilter) 
library(ggrepel)
library(eulerr)
library(DESeq2)
library(UpSetR)
library(GenomicFeatures)
library(biomaRt)
library(tidyverse)
library(tximport)
library(ashr)



#######################################
#Load Counts  ---- 
#######################################


# Choose the files 
TEcounts = read_table("../pathToYour/TETranscriptsOutput.cntTable")

# manually add names to columns

colnames(TEcounts)[1] <- c("TE") 
  
# tidy up file names (will depend on your own naming)
#colnames(TEcounts)[2:7] <- str_replace(colnames(TEcounts)[2:7],
  #                                     '.+5_STARTEtranscripts//(.+)', '\\1')
#colnames(TEcounts)

# to remove everything after replicate number
#colnames(TEcounts)[2:7] <- str_replace(colnames(TEcounts)[2:7],
 #                                      '(.+)Aligned.sortedByCoord.out.bam.+', '\\1')
  
colnames(TEcounts)


###############################################
# make sample info  and TE counts matrix ---- #
###############################################

sampleInfoTE = as.data.frame(colnames(TEcounts[2:7]))
#sampleInfoTE$group = sub("_.*", "", sampleInfoTE$`colnames(TEcounts[2:7])`)
sampleInfoTE$group = sub("_r.*", "", sampleInfoTE$`colnames(TEcounts[2:7])`)
sampleInfoTE$group = sub(".*_", "", sampleInfoTE$group)

colnames(sampleInfoTE) = c("sample", "group")
row.names(sampleInfoTE) = sampleInfoTE$sample

TEcounts_matrix <- as.data.frame(TEcounts)
rownames(TEcounts_matrix) <- TEcounts_matrix$TE
head(TEcounts_matrix)

TEcounts_matrix <- TEcounts_matrix[,2:7]
head(TEcounts_matrix)

### for GEO sub raw counts
write.table(TEcounts, file = "rawCounts_mRNA_GenesAndTEs.txt",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


#######################################
# Run DESeq2  ---- 
#######################################

ddsObjTE.raw <- DESeqDataSetFromMatrix(countData = TEcounts_matrix,
                                       colData = sampleInfoTE,
                                       design = ~ group)

ddsObjTE <- DESeq(ddsObjTE.raw)

plotDispEsts(ddsObjTE)

rlogObjTE = rlog(ddsObjTE)

#######################################
# save normalised and not normalised counts  ---- 
#######################################
normTECounts = assay(rlogObjTE)

notNormTECounts = counts(ddsObjTE, normalized = F)


## save as tables only with TE expression
normTECounts
nrow(normTECounts)
normTECounts_save <-normTECounts[grep(":",rownames(normTECounts)),]
nrow(normTECounts_save)

View(normTECounts_save)

write.table(normTECounts_save,
            "NormTEcounts.txt",sep="\t",
            col.names=T, quote=F, row.names=T)



write.table(notNormTECounts,
            "NotNormTEcounts.txt",sep="\t",
            col.names=T, quote=F, row.names=T)

#######################################
# PCA  ---- 
#######################################

plotPCAEx(rlogObjTE,PCx = 1,PCy = 2,cond = "group",ntop = 500, T)
plotPCAEx(rlogObjTE,PCx = 2,PCy = 3,cond = "group",ntop = 500, T)


