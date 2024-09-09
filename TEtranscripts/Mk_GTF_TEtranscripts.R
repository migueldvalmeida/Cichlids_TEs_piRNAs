library(tidyverse)

### Created by Alexandra Dallaire, received 28/09/2021. Modified slightly by Miguel Almeida.

#Make GTF file to input in TEtranscripts

######################################################################
############### TE GTF for A. calliptera  ############################
######################################################################

#Make a TE GTF file, to be transformed into a gff file. It will be the input for TEtranscripts (Molly Hammell; https://github.com/mhammell-laboratory/TEtranscripts)
annotTE <- read.table("./pathToYourRepeatMaskerOutput/*.toplevel.fa.out", header = FALSE, fill=TRUE)
View(annotTE)
annotTE <- na.omit(annotTE)


TE1 <- subset(annotTE, select=c(V5))
TE1$V2 <- "Acal_rmsk"
TE1$V3 <- "exon"
TE1$V4 <- annotTE$V6
TE1$V6 <- annotTE$V7
TE1$V7 <- "."
TE1$V8 <- annotTE$V9 
TE1$V9 <- "."

View(TE1)

#Print all attribute types in V11. Make sure I don't delete composite names with grepl
att <- data.frame(unique(annotTE$V11))

#Remove unknown repeats, simple repeats, rRNAs, scRNAs, snRNAs, srpRNAs, and tRNAs. 
TE1$V10 = annotTE$V10
TE1$V11 = annotTE$V11
TE2 <- dplyr::filter(TE1, !grepl("Unknown|Low_complexity|Simple_repeat|Satellite|tRNA|rRNA|snRNA", V11, ignore.case=TRUE))
kept <- data.frame(unique(TE2$V11))
#Remove 'C's from orientation col
TE2$V8 <- gsub("\\C.*", "-", TE2$V8)

#Play with the last column. Aim =
#gene_id "NINJA_I"; transcript_id "NINJA_I-[uniqueIDnumber]"; family_id "Pao"; class_id "LTR";

TE3 <- TE2
string1 = "gene_id \""
TE3$V12 = paste(string1, TE3$V10, sep="")
string2 = "\"; transcript_id \""
TE3$V13 = paste(TE3$V12, string2, sep="")
TE3$V14 <- paste(TE3$V13, TE3$V10, sep="")
#Give transcripts unique IDs. THESE ARE NOT THE SAME IDs as in the initial RM output.
TE3$newnumber = seq.int(nrow(TE3))
TE3$V15 <- paste(TE3$V14, "-", TE3$newnumber, sep="")
string3 = "\"; family_id \""
TE3$V16 <- paste(TE3$V15, string3, sep="")

TE3$V12 <- NULL
TE3$V13 <- NULL
TE3$V14 <- NULL
TE3$newnumber <- NULL

#Split the family and class
library(stringr)
TE4 <- data.frame(str_split_fixed(TE3$V11, "/", 2))
#If TE family is unknown, replace by Class
TE4$X2[TE4$X2==""]<-"Unknown"

# Add family and class to TE3 df
TE3$V17 <- paste(TE3$V16, TE4$X2, sep="")
string4 = "\"; class_id \""
TE3$V18 <- paste(TE3$V17, string4, sep="")
TE3$V19 <- paste(TE3$V18, TE4$X1, sep="")
string5 = "\" ;"
TE3$V20 <- paste(TE3$V19, string5, sep="")

TE3$V10 <- NULL
TE3$V11 <- NULL
TE3$V15 <- NULL
TE3$V16 <- NULL
TE3$V17 <- NULL
TE3$V18 <- NULL
TE3$V19 <- NULL

out_dir <- getwd()
write.table(TE3, paste(out_dir,"YourFileName.TEtranscripts.gtf",sep="/"), col.names=F, quote=F, sep="\t", row.names=F)






