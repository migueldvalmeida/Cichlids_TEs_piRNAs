## libraries needed:
library(RColorBrewer)
library(Gviz)
library(GenomicFeatures)
library(data.table)


gen <- 'fAstCal1.2'


## gene coordinates
mygenes <- fread("genes_to_plot.txt") # a tab separated file with different columns: "chromosome", "start", "end" and "name"
gene <- mygenes[1] ## change row according to which gene on the table you want to plot

chr <- gene[["chromosome"]]
mygene_start <- gene[["start"]]
mygene_end <- gene[["end"]]
mygene <- gene[["name"]]

# create track with gene annotations:
txdb <- makeTxDbFromGFF("Astatotilapia_calliptera.fAstCal1.2.104.gtf", format="gtf") # gene annotation for Astatotilapia calliptera obtained from Ensembl.
options(ucscChromosomeNames=FALSE)
grtrack <- GeneRegionTrack(
  txdb,
  genome=gen,
  chromosome=chr,
  name='fAstCal1.2 genes',
  collapseTranscripts=FALSE,
  shape="smallArrow",
  stacking="squish",
  # start=mygene_start,
  # end=mygene_end,
  showId = TRUE,
  fill="#000000",
  col= NULL,
  col.line=NULL
)

txdb2 <- makeTxDbFromGFF("Astatotilapia_calliptera.fAstCal1.2.TEtranscripts_PiosLibrary_v3.2.gtf", format="gtf") # curated TE annotation for A. calliptera. 
options(ucscChromosomeNames=FALSE)
TEtrack <- GeneRegionTrack(
  txdb2,
  genome=gen,
  chromosome=chr,
  name='TEs',
  collapseTranscripts=FALSE,
  shape="box",
  stacking="dense",
  # start=mygene_start,
  # end=mygene_end,
  showId = TRUE,
  fill="#000000",
  col= NULL,
  col.line=NULL
)


bw1 <- "bigWigFile1.bw"
bw2 <- "bigWigFile2.bw"
bw3 <- "bigWigFile3.bw"
bw4 <- "bigWigFile4.bw"

name1 <- "Name Bigwig file 1"
name2 <- "Name Bigwig file 2"
name3 <- "Name Bigwig file 3"
name4 <- "Name Bigwig file 4"

bw1Track <- DataTrack(
  range = bw1,
  genome = gen,
  chromosome = chr,
  name = name1,
  type = "histogram",
  col.histogram="#3f007d",
  fill="#3f007d",
  ylim=c(0,4))



bw2Track <- DataTrack(
  range = bw2,
  genome = gen,
  chromosome = chr,
  name = name2,
  type = "histogram",
  col.histogram="#9e9ac8",
  fill="#9e9ac8",
  ylim=c(0,4))




bw3Track <- DataTrack(
  range = bw3,
  genome = gen,
  chromosome = chr,
  name = name3,
  type = "histogram",
  col.histogram="#0570b0",
  fill="#0570b0",
  ylim=c(0,4))


bw4Track <- DataTrack(
  range = bw4,
  genome = gen,
  chromosome = chr,
  name = name4,
  type = "histogram",
  col.histogram="#006d2c",
  fill="#006d2c",
  ylim=c(0,4))



## adds position in chromosome
axisTrack <- GenomeAxisTrack()


track_list <- list(
  bw1Track,
  bw2Track,
  bw3Track,
  bw4Track,
  axisTrack,
  grtrack,
  TEtrack
)

plotTracks(track_list,
           chromosome=chr,
           from = mygene_start - 5000,  # adjust number according to desired window range.
           to = mygene_end + 5000,      # adjust number according to desired window range.
           background.title="white",
           col.title="black",
           col.axis="black")




