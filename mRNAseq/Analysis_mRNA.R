####################################################################################################
# Functions for the analysis of RNASeq data (or data processed with DESeq2)
# source("rnaSeqAnalysisFuncs.R") to get the functions in the environment
# Adapted by Miguel Almeida from a script by Jon Price. 
####################################################################################################

#######################################
#Load libraries  ---- 
#######################################

library(ggplot2)
library(reshape2)
library(genefilter) 
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

# Get the tx2gene table
tx2gene = read.table("../PathToYourTx2geneFile/tx2gene.tsv", header = T)

# Choose the files 
files = list.files(pattern = ".*sf", 
                   path =  "../pathToYour/SalmonResults/", 
                   recursive = T,
                   full.names = T)
# add names
names(files) = list.files(path =  "../pathToYour/SalmonResults/")


# Remove a sample if nessecary 
#files = files[-12]

# scaled using the average transcript length over samples and then the library size (lengthScaledTPM)
#nice 
tpm <- tximport(files, type = "salmon",
                tx2gene = tx2gene, 
                countsFromAbundance = "lengthScaledTPM")


# This is the equivelant of raw counts
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                countsFromAbundance = "no")


#######################################
# make sample info  ---- 
#######################################

sampleInfo = as.data.frame(colnames(txi$counts))
sampleInfo$group <- sub(".*_ov_.*", "ov", sub(".*_tes_.*", "tes", sampleInfo$`colnames(txi$counts)`))

colnames(sampleInfo) = c("sample", "group")
row.names(sampleInfo) = sampleInfo$sample



#######################################
# Make annotation   ---- 
#######################################

# See what annotations are available
listMarts()

# Set ENSEMBL as the annotation and see what species are there
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
head(datasets, 50)

# choose O. niloticus, as an example, change dataset according to your needs.
ensembl = useDataset("oniloticus_gene_ensembl",mart=ensembl)

# Filters are the info we have 
filters = listFilters(ensembl)
filters[1:5,]

# attributes are the information we want to get 
attributes = listAttributes(ensembl)
attributes[1:50,]


# Get the annotation 
annotation = getBM(attributes=c('ensembl_gene_id',
                                'description' , 'gene_biotype', 
                                'external_gene_name','external_gene_source',
                                'entrezgene_accession',"chromosome_name","start_position"), 
                   filters = 'ensembl_gene_id', 
                   values = row.names(txi$abundance), 
                   mart = ensembl, uniqueRows = T)

# See if any have duplicate entries 
annotation[duplicated(annotation$ensembl_gene_id),]
annotation
dim(annotation)


# isolate protein coding genes
protein_codingGenes = annotation[annotation$gene_biotype == "protein_coding",]
nrow(protein_codingGenes)


# If you want to add other annotation sets:
#t1 = read.table("../info/morc1/exampleGeneSet.txt", header = F, sep = "\t")
#protein_codingGenes$example = 0
#protein_codingGenes[protein_codingGenes$ensembl_gene_id %in% t1$V1,]$example = 1




#######################################
# Run DESeq2  ---- 
#######################################

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleInfo,
                                       design = ~ group)
ddsObj.raw = ddsObj.raw[protein_codingGenes$ensembl_gene_id,]

ddsObj <- DESeq(ddsObj.raw)

plotDispEsts(ddsObj)

rlogObj = rlog(ddsObj)

#######################################
# save normalised and not normalised counts  ---- 
#######################################
normCounts = assay(rlogObj)
normCountsAnnot = cbind(protein_codingGenes, normCounts )

notNormCounts = counts(ddsObj, normalized = F)
notNormCountsAnnot = cbind(protein_codingGenes, notNormCounts  )


#######################################
# PCA  ---- 
#######################################

plotPCAEx(rlogObj,PCx = 1,PCy = 2,cond = "group",ntop = 500, T)
plotPCAEx(rlogObj,PCx = 2,PCy = 3,cond = "group",ntop = 500, T)



### TABLES OF COUNTS FOR ALL GENES

### get table with gene names , chromosome, normalized reads

#View(normCountsAnnot)
normCounts_parsed <- normCountsAnnot[,c(1,4,7,9:14)]
length(normCounts_parsed$ensembl_gene_id)
#View(normCounts_parsed)


write.table(normCounts_parsed, file = "NormCounts.txt", quote = FALSE)


### Get a file with the annotations

AnnotationGenes <- normCountsAnnot[,c(1,2,3,4,7,8)]

write.table(AnnotationGenes, file = "AnnotationGenes_ensembl.txt",
            quote = FALSE, row.names = FALSE)




normCounts_parsed

normCounts_parsed2 <- normCounts_parsed

normCounts_parsed2 <- normCounts_parsed2[c(1,4:9)]

normCounts_parsed2m <- melt(normCounts_parsed2)
normCounts_parsed2m
#View(normCounts_parsed2m)


normCounts_parsed2m <- mutate(normCounts_parsed2m, tissue_terms = str_extract(variable, "(ov|tes)"))
normCounts_parsed2m$species = "[YourSpecies/sampleName]"

write.table(normCounts_parsed2m, file = "NormCounts_AllGenes_melted.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


write.table(tpm$abundance, file = "LengthScaled_TPM.txt", sep = "\t", col.names = TRUE, quote = FALSE)

