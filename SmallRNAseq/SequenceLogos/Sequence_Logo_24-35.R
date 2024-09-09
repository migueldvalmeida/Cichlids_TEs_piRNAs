library(tidyverse)
library(reshape2)
library(ggseqlogo)
library(phylotools)

###############################################
### Load Tables and make data tidy ---- #######
###############################################


BigFasta_AcaTes <- phylotools::read.fasta(file = "./Aca_tes_Allreps.fasta", clean_name = TRUE)
#View(BigFasta_AcaTes)

# make the Ts into Us
BigFasta_AcaTes$seq.text <- gsub("T", "U", BigFasta_AcaTes$seq.text)


###############################################
### Plot ---- #################################
###############################################

plot1 <- ggplot () +
  geom_logo(BigFasta_AcaTes$seq.text, method = 'prob', seq_type='rna') +
  annotate('rect', xmin = 9.5, xmax = 10.5, ymin = -0.05, ymax = 1.05, alpha = .1, col='black', fill='yellow') +
  theme_logo()

pdf("Aca_Tes_seqLogo_prob.pdf", width = 10.79, height = 3.65)
plot(plot1)
dev.off()


plot2 <- ggplot () +
  geom_logo(BigFasta_AcaTes$seq.text, method = 'bits', seq_type='rna') +
  theme_logo()

pdf("Aca_Tes_seqLogo_bits.pdf", width = 10.79, height = 3.65)
plot(plot2)
dev.off()

## print with 10.79 X 3.65 inches

######################################################
### Load Tables and make data tidy Aca_ov ---- #######
######################################################


BigFasta_AcaOv <- phylotools::read.fasta(file = "./Aca_ov_Allreps.fasta", clean_name = TRUE)
#View(BigFasta_AcaOv)

# make the Ts into Us
BigFasta_AcaOv$seq.text <- gsub("T", "U", BigFasta_AcaOv$seq.text)

######################################################
### Plot Aca_ov ---- #################################
######################################################


plot3 <- ggplot () +
  geom_logo(BigFasta_AcaOv$seq.text, method = 'prob', seq_type='rna') +
  annotate('rect', xmin = 9.5, xmax = 10.5, ymin = -0.05, ymax = 1.05, alpha = .1, col='black', fill='yellow') +
  theme_logo()

pdf("Aca_Ov_seqLogo_prob.pdf", width = 10.79, height = 3.65)
plot(plot3)
dev.off()


plot4 <- ggplot () +
  geom_logo(BigFasta_AcaOv$seq.text, method = 'bits', seq_type='rna') +
  theme_logo()

pdf("Aca_Ov_seqLogo_bits.pdf", width = 10.79, height = 3.65)
plot(plot4)
dev.off()



#######################################################
### Load Tables and make data tidy Aca_mus ---- #######
#######################################################


BigFasta_AcaMus <- phylotools::read.fasta(file = "./Aca_mus_rep2.fasta", clean_name = TRUE)
#View(BigFasta_AcaMus)

# make the Ts into Us
BigFasta_AcaMus$seq.text <- gsub("T", "U", BigFasta_AcaMus$seq.text)


######################################################
### Plot Aca_mus ---- ################################
######################################################


plot5 <- ggplot () +
  geom_logo(BigFasta_AcaMus$seq.text, method = 'prob', seq_type='rna') +
  annotate('rect', xmin = 9.5, xmax = 10.5, ymin = -0.05, ymax = 1.05, alpha = .1, col='black', fill='yellow') +
  theme_logo()

pdf("Aca_Mus_seqLogo_prob.pdf", width = 10.79, height = 3.65)
plot(plot5)
dev.off()


plot6 <- ggplot () +
  geom_logo(BigFasta_AcaMus$seq.text, method = 'bits', seq_type='rna') +
  theme_logo()

pdf("Aca_Mus_seqLogo_bits.pdf", width = 10.79, height = 3.65)
plot(plot6)
dev.off()

#######################################################
### Load Tables and make data tidy Mze_tes ---- #######
#######################################################


BigFasta_MzeTes <- phylotools::read.fasta(file = "./Mze_tes_Allreps.fasta", clean_name = TRUE)
#View(BigFasta_MzeTes)

# make the Ts into Us
BigFasta_MzeTes$seq.text <- gsub("T", "U", BigFasta_MzeTes$seq.text)



######################################################
### Plot Mze_tes ---- ################################
######################################################


plot7 <- ggplot () +
  geom_logo(BigFasta_MzeTes$seq.text, method = 'prob', seq_type='rna') +
  annotate('rect', xmin = 9.5, xmax = 10.5, ymin = -0.05, ymax = 1.05, alpha = .1, col='black', fill='yellow') +
  theme_logo()

pdf("Mze_Tes_seqLogo_prob.pdf", width = 10.79, height = 3.65)
plot(plot7)
dev.off()


plot8 <- ggplot () +
  geom_logo(BigFasta_MzeTes$seq.text, method = 'bits', seq_type='rna') +
  theme_logo()

pdf("Mze_Tes_seqLogo_bits.pdf", width = 10.79, height = 3.65)
plot(plot8)
dev.off()


#######################################################
### Load Tables and make data tidy Tma_tes ---- #######
#######################################################


BigFasta_TmaTes <- phylotools::read.fasta(file = "./Tma_tes_Allreps.fasta", clean_name = TRUE)
#View(BigFasta_TmaTes)

# make the Ts into Us
BigFasta_TmaTes$seq.text <- gsub("T", "U", BigFasta_TmaTes$seq.text)



######################################################
### Plot Tma_tes ---- ################################
######################################################


plot9 <- ggplot () +
  geom_logo(BigFasta_TmaTes$seq.text, method = 'prob', seq_type='rna') +
  annotate('rect', xmin = 9.5, xmax = 10.5, ymin = -0.05, ymax = 1.05, alpha = .1, col='black', fill='yellow') +
  theme_logo()

pdf("Tma_Tes_seqLogo_prob.pdf", width = 10.79, height = 3.65)
plot(plot9)
dev.off()


plot10 <- ggplot () +
  geom_logo(BigFasta_TmaTes$seq.text, method = 'bits', seq_type='rna') +
  theme_logo()

pdf("Tma_Tes_seqLogo_bits.pdf", width = 10.79, height = 3.65)
plot(plot10)
dev.off()



