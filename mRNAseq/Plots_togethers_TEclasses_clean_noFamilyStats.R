library(tidyverse)
library(pheatmap)
library(reshape2)
library(ggrepel)
library(UpSetR)
library(RColorBrewer)
library(ggpubr)
library(patchwork)
library(ggplotify)
library(gplots)




##########################################################################################
########  Load Relevant files and boxplots per family and class usign ALL COUNTS ----  ###
##########################################################################################

### First for A. CALLIPTERA  ----

Aca_normTEcounts <- read.table(file = "Acal_NormTEcounts.txt", sep = "\t")
Aca_normTEcounts
#View(Aca_normTEcounts)
nrow(Aca_normTEcounts)

normTECounts2_Aca <- as.data.frame(Aca_normTEcounts)
normTECounts2_Aca$TE <- rownames(normTECounts2_Aca)

rownames(normTECounts2_Aca) <- NULL


normTECounts3_Aca <- normTECounts2_Aca[,c(9,1:8)]
#View(normTECounts3_Aca)
nrow(normTECounts3_Aca)
normTECounts3m_Aca <- melt(normTECounts3_Aca)
#View(normTECounts3m_Aca)
nrow(normTECounts3m_Aca)


normTECounts3m_Aca$TEclass <-  sub(".*:", "", normTECounts3m_Aca$TE)
normTECounts3m_Aca$TEfam <- sub("^[^:]*:([^:]*).*", "\\1", normTECounts3m_Aca$TE)
#View(normTECounts3m_Aca)

normTECounts3m_Aca$tissue <-  sub("^[^_]*_([^_]*).*", "\\1", normTECounts3m_Aca$variable)
#View(normTECounts3m_Aca)

nrow(normTECounts3m_Aca)
normTECounts3m_Aca_final <- filter(normTECounts3m_Aca, !TEclass == "Retroposon" & !TEclass == "SINE?")
nrow(normTECounts3m_Aca_final)


unique(normTECounts3m_Aca_final$TEclass)

normTECounts3m_Aca_final$TEclass <- ordered(normTECounts3m_Aca_final$TEclass,
                                        levels=c("SINE","RC","LTR","LINE","DNA"))



plotTECLASSES_Aca <- ggplot(normTECounts3m_Aca_final, aes(x = TEclass, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,18)) + 
  labs(title = NULL, y = "Log-Normalized counts", x = "TE class") + 
  scale_fill_brewer(palette = "Oranges")+
  theme_bw()

plotTECLASSES_Aca2 <- plotTECLASSES_Aca + coord_flip()
plotTECLASSES_Aca2

# stats
fortestClassesDNA_Aca <- normTECounts3m_Aca_final[grepl("DNA", normTECounts3m_Aca_final$TEclass),]
nrow(fortestClassesDNA_Aca)

fortestClassesLINE_Aca <- normTECounts3m_Aca_final[grepl("LINE", normTECounts3m_Aca_final$TEclass),]
nrow(fortestClassesLINE_Aca)

fortestClassesLTR_Aca <- normTECounts3m_Aca_final[grepl("LTR", normTECounts3m_Aca_final$TEclass),]
nrow(fortestClassesLTR_Aca)

fortestClassesRC_Aca <- normTECounts3m_Aca_final[grepl("RC", normTECounts3m_Aca_final$TEclass),]
nrow(fortestClassesRC_Aca)

fortestClassesSINE_Aca <- normTECounts3m_Aca_final[grepl("SINE", normTECounts3m_Aca_final$TEclass),]
nrow(fortestClassesSINE_Aca)


testDNAAca <- compare_means(value ~ tissue, data = fortestClassesDNA_Aca, p.adjust.method = "BH")
View(testDNAAca)

testLINEAca <- compare_means(value ~ tissue, data = fortestClassesLINE_Aca, p.adjust.method = "BH")
View(testLINEAca)

testLTRAca <- compare_means(value ~ tissue, data = fortestClassesLTR_Aca, p.adjust.method = "BH")
View(testLTRAca)

testRCAca <- compare_means(value ~ tissue, data = fortestClassesRC_Aca, p.adjust.method = "BH")
View(testRCAca)

testSINEAca <- compare_means(value ~ tissue, data = fortestClassesSINE_Aca, p.adjust.method = "BH")
View(testSINEAca)



### PER FAMILY -> DNA Acalliptera ---- 

normTEcounts_perfamALL_DNAm_Aca <- filter(normTECounts3m_Aca_final, TEclass == "DNA")
#View(normTEcounts_perfamALL_DNAm_Aca)

levels(factor(normTEcounts_perfamALL_DNAm_Aca$TEfam))
normTEcounts_perfamALL_DNAm_Aca$TEfam <- fct_rev(normTEcounts_perfamALL_DNAm_Aca$TEfam)
levels(factor(normTEcounts_perfamALL_DNAm_Aca$TEfam))

plotTEFamDNA_Aca <- ggplot(normTEcounts_perfamALL_DNAm_Aca, aes(x = TEfam, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.5) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,16)) + 
  labs(x = NULL, y = "Log-Normalized counts", title = "DNA transposons") + 
  #labs(title = "Brood size at 25°C", y = "Brood Size", x = NULL) + 
  scale_fill_brewer(palette = "Oranges")+
  theme_bw()

plotTEFamDNA_Aca2 <- plotTEFamDNA_Aca + coord_flip()
plotTEFamDNA_Aca2


# Retrotransposons A calliptera ----

normTEcounts_perfamALL_retro1_Aca <- filter(normTECounts3m_Aca_final, TEclass == "LINE")
#View(normTEcounts_perfamALL_retro1_Aca)

normTEcounts_perfamALL_retro2_Aca <- filter(normTECounts3m_Aca_final, TEclass == "LTR")
#View(normTEcounts_perfamALL_retro2_Aca)

normTEcounts_perfamALL_retro3_Aca <- filter(normTECounts3m_Aca_final, TEclass == "Retroposon")
#View(normTEcounts_perfamALL_retro3_Aca)

normTEcounts_perfamALL_retro4_Aca <- filter(normTECounts3m_Aca_final, TEclass == "SINE")
#View(normTEcounts_perfamALL_retro4_Aca)

normTEcounts_perfamALL_retroTEs_Aca <- rbind(normTEcounts_perfamALL_retro1_Aca,normTEcounts_perfamALL_retro2_Aca,
                                         normTEcounts_perfamALL_retro3_Aca,normTEcounts_perfamALL_retro4_Aca)

#View(normTEcounts_perfamALL_retroTEs_Aca)


levels(factor(normTEcounts_perfamALL_retroTEs_Aca$TEfam))
normTEcounts_perfamALL_retroTEs_Aca$TEfam <- fct_rev(normTEcounts_perfamALL_retroTEs_Aca$TEfam)
levels(factor(normTEcounts_perfamALL_retroTEs_Aca$TEfam))

plotTEFamRetro_Aca <- ggplot(normTEcounts_perfamALL_retroTEs_Aca, aes(x = TEfam, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.5) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,16)) + 
  labs(x = NULL, y = "Log-Normalized counts", title = "Retrotransposons") + 
  #labs(title = "Brood size at 25°C", y = "Brood Size", x = NULL) + 
  scale_fill_brewer(palette = "Oranges")+
  theme_bw()

plotTEFamRetro_Aca2 <- plotTEFamRetro_Aca + coord_flip()
plotTEFamRetro_Aca2



testRetroFAMsAca <- compare_means(value ~ tissue, data = normTEcounts_perfamALL_retroTEs_Aca, group.by = "TEfam", p.adjust.method = "BH")
View(testRetroFAMsAca)

testDNAFAMsAca <- compare_means(value ~ tissue, data = normTEcounts_perfamALL_DNAm_Aca, group.by = "TEfam", p.adjust.method = "BH")
View(testDNAFAMsAca)



### Then for A. CALLIPTERA  with counts obtained when a curated TE annotation is used  ----

AcaPio_normTEcounts <- read.table(file = "Malawi_NormTEcountsPio.txt", sep = "\t")
AcaPio_normTEcounts
#View(AcaPio_normTEcounts)
nrow(AcaPio_normTEcounts)

AcaPio_normTEcounts <- AcaPio_normTEcounts[,c(1:2,5,10:14)]

normTECounts2_AcaPio <- as.data.frame(AcaPio_normTEcounts)
normTECounts2_AcaPio$TE <- rownames(normTECounts2_AcaPio)

rownames(normTECounts2_AcaPio) <- NULL


normTECounts3_AcaPio <- normTECounts2_AcaPio[,c(9,1:8)]
#View(normTECounts3_AcaPio)
nrow(normTECounts3_AcaPio)
normTECounts3m_AcaPio <- melt(normTECounts3_AcaPio)
#View(normTECounts3m_AcaPio)
nrow(normTECounts3m_AcaPio)


normTECounts3m_AcaPio$TEclass <-  sub(".*:", "", normTECounts3m_AcaPio$TE)
normTECounts3m_AcaPio$TEfam <- sub("^[^:]*:([^:]*).*", "\\1", normTECounts3m_AcaPio$TE)
#View(normTECounts3m_AcaPio)

normTECounts3m_AcaPio$tissue <-  sub("^[^_]*_([^_]*).*", "\\1", normTECounts3m_AcaPio$variable)
#View(normTECounts3m_AcaPio)

nrow(normTECounts3m_AcaPio)
normTECounts3m_AcaPio_final <- filter(normTECounts3m_AcaPio, !TEclass == "Retroposon" & !TEclass == "SINE?")
nrow(normTECounts3m_AcaPio_final)


unique(normTECounts3m_AcaPio_final$TEclass)
unique(normTECounts3m_AcaPio_final$TEfam)

normTECounts3m_AcaPio_final$TEclass <- ordered(normTECounts3m_AcaPio_final$TEclass,
                                            levels=c("SINE","RC","LTR","LINE","DNA"))



plotTECLASSES_AcaPio <- ggplot(normTECounts3m_AcaPio_final, aes(x = TEclass, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,18)) + 
  labs(title = NULL, y = "Log-Normalized counts", x = "TE class") + 
  scale_fill_brewer(palette = "Oranges")+
  theme_bw()

plotTECLASSES_AcaPio2 <- plotTECLASSES_AcaPio + coord_flip()
plotTECLASSES_AcaPio2

# stats
fortestClassesDNA_AcaPio <- normTECounts3m_AcaPio_final[grepl("DNA", normTECounts3m_AcaPio_final$TEclass),]
nrow(fortestClassesDNA_AcaPio)

fortestClassesLINE_AcaPio <- normTECounts3m_AcaPio_final[grepl("LINE", normTECounts3m_AcaPio_final$TEclass),]
nrow(fortestClassesLINE_AcaPio)

fortestClassesLTR_AcaPio <- normTECounts3m_AcaPio_final[grepl("LTR", normTECounts3m_AcaPio_final$TEclass),]
nrow(fortestClassesLTR_AcaPio)

fortestClassesRC_AcaPio <- normTECounts3m_AcaPio_final[grepl("RC", normTECounts3m_AcaPio_final$TEclass),]
nrow(fortestClassesRC_AcaPio)

fortestClassesSINE_AcaPio <- normTECounts3m_AcaPio_final[grepl("SINE", normTECounts3m_AcaPio_final$TEclass),]
nrow(fortestClassesSINE_AcaPio)


testDNAAcaPio <- compare_means(value ~ tissue, data = fortestClassesDNA_AcaPio, p.adjust.method = "BH")
View(testDNAAcaPio)

testLINEAcaPio <- compare_means(value ~ tissue, data = fortestClassesLINE_AcaPio, p.adjust.method = "BH")
View(testLINEAcaPio)

testLTRAcaPio <- compare_means(value ~ tissue, data = fortestClassesLTR_AcaPio, p.adjust.method = "BH")
View(testLTRAcaPio)

testRCAcaPio <- compare_means(value ~ tissue, data = fortestClassesRC_AcaPio, p.adjust.method = "BH")
View(testRCAcaPio)

testSINEAcaPio <- compare_means(value ~ tissue, data = fortestClassesSINE_AcaPio, p.adjust.method = "BH")
View(testSINEAcaPio)



### PER FAMILY -> DNA AcaPiolliptera ---- 

normTEcounts_perfamALL_DNAm_AcaPio <- filter(normTECounts3m_AcaPio_final, TEclass == "DNA")
#View(normTEcounts_perfamALL_DNAm_AcaPio)

levels(factor(normTEcounts_perfamALL_DNAm_AcaPio$TEfam))
normTEcounts_perfamALL_DNAm_AcaPio$TEfam <- fct_rev(normTEcounts_perfamALL_DNAm_AcaPio$TEfam)
levels(factor(normTEcounts_perfamALL_DNAm_AcaPio$TEfam))

plotTEFamDNA_AcaPio <- ggplot(normTEcounts_perfamALL_DNAm_AcaPio, aes(x = TEfam, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.5) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,16)) + 
  labs(x = NULL, y = "Log-Normalized counts", title = "DNA transposons") + 
  #labs(title = "Brood size at 25°C", y = "Brood Size", x = NULL) + 
  scale_fill_brewer(palette = "Oranges")+
  theme_bw()

plotTEFamDNA_AcaPio2 <- plotTEFamDNA_AcaPio + coord_flip()
plotTEFamDNA_AcaPio2


# Retrotransposons A calliptera ----

normTEcounts_perfamALL_retro1_AcaPio <- filter(normTECounts3m_AcaPio_final, TEclass == "LINE")
#View(normTEcounts_perfamALL_retro1_AcaPio)

normTEcounts_perfamALL_retro2_AcaPio <- filter(normTECounts3m_AcaPio_final, TEclass == "LTR")
#View(normTEcounts_perfamALL_retro2_AcaPio)

normTEcounts_perfamALL_retro3_AcaPio <- filter(normTECounts3m_AcaPio_final, TEclass == "Retroposon")
#View(normTEcounts_perfamALL_retro3_AcaPio)

normTEcounts_perfamALL_retro4_AcaPio <- filter(normTECounts3m_AcaPio_final, TEclass == "SINE")
#View(normTEcounts_perfamALL_retro4_AcaPio)

normTEcounts_perfamALL_retroTEs_AcaPio <- rbind(normTEcounts_perfamALL_retro1_AcaPio,normTEcounts_perfamALL_retro2_AcaPio,
                                             normTEcounts_perfamALL_retro3_AcaPio,normTEcounts_perfamALL_retro4_AcaPio)

#View(normTEcounts_perfamALL_retroTEs_AcaPio)


levels(factor(normTEcounts_perfamALL_retroTEs_AcaPio$TEfam))
normTEcounts_perfamALL_retroTEs_AcaPio$TEfam <- fct_rev(normTEcounts_perfamALL_retroTEs_AcaPio$TEfam)
levels(factor(normTEcounts_perfamALL_retroTEs_AcaPio$TEfam))

plotTEFamRetro_AcaPio <- ggplot(normTEcounts_perfamALL_retroTEs_AcaPio, aes(x = TEfam, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.5) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,16)) + 
  labs(x = NULL, y = "Log-Normalized counts", title = "Retrotransposons") + 
  #labs(title = "Brood size at 25°C", y = "Brood Size", x = NULL) + 
  scale_fill_brewer(palette = "Oranges")+
  theme_bw()

plotTEFamRetro_AcaPio2 <- plotTEFamRetro_AcaPio + coord_flip()
plotTEFamRetro_AcaPio2





testRetroFAMsAcaPio <- compare_means(value ~ tissue, data = normTEcounts_perfamALL_retroTEs_AcaPio, group.by = "TEfam", p.adjust.method = "BH")
View(testRetroFAMsAcaPio)

testDNAFAMsAcaPio <- compare_means(value ~ tissue, data = normTEcounts_perfamALL_DNAm_AcaPio, group.by = "TEfam", p.adjust.method = "BH")
View(testDNAFAMsAcaPio)




### A. BURTONI  ----

Abu_normTEcounts <- read.table(file = "Abur_NormTEcounts.txt", sep = "\t")
Abu_normTEcounts
#View(Abu_normTEcounts)
nrow(Abu_normTEcounts)

normTECounts2_Abu <- as.data.frame(Abu_normTEcounts)
normTECounts2_Abu$TE <- rownames(normTECounts2_Abu)

rownames(normTECounts2_Abu) <- NULL


normTECounts3_Abu <- normTECounts2_Abu[,c(8,1:7)]
#View(normTECounts3_Abu)
nrow(normTECounts3_Abu)
normTECounts3m_Abu <- melt(normTECounts3_Abu)
#View(normTECounts3m_Abu)
nrow(normTECounts3m_Abu)


normTECounts3m_Abu$TEclass <-  sub(".*:", "", normTECounts3m_Abu$TE)
normTECounts3m_Abu$TEfam <- sub("^[^:]*:([^:]*).*", "\\1", normTECounts3m_Abu$TE)
#View(normTECounts3m_Abu)

normTECounts3m_Abu$tissue <-  sub("^[^_]*_([^_]*).*", "\\1", normTECounts3m_Abu$variable)
#View(normTECounts3m_Abu)

nrow(normTECounts3m_Abu)
normTECounts3m_Abu_final <- filter(normTECounts3m_Abu, !TEclass == "Retroposon" & !TEclass == "SINE?")
nrow(normTECounts3m_Abu_final)


unique(normTECounts3m_Abu_final$TEclass)

normTECounts3m_Abu_final$TEclass <- ordered(normTECounts3m_Abu_final$TEclass,
                                            levels=c("SINE","RC","LTR","LINE","DNA"))


#View(normTECounts3m_Abu_final)

plotTECLASSES_Abu <- ggplot(normTECounts3m_Abu_final, aes(x = TEclass, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,18)) + 
  labs(title = NULL, y = "Log-Normalized counts", x = "TE class") + 
  scale_fill_brewer(palette = "Blues")+
  theme_bw()

plotTECLASSES_Abu2 <- plotTECLASSES_Abu + coord_flip()
plotTECLASSES_Abu2

# stats
fortestClassesDNA_Abu <- normTECounts3m_Abu_final[grepl("DNA", normTECounts3m_Abu_final$TEclass),]
nrow(fortestClassesDNA_Abu)

fortestClassesLINE_Abu <- normTECounts3m_Abu_final[grepl("LINE", normTECounts3m_Abu_final$TEclass),]
nrow(fortestClassesLINE_Abu)

fortestClassesLTR_Abu <- normTECounts3m_Abu_final[grepl("LTR", normTECounts3m_Abu_final$TEclass),]
nrow(fortestClassesLTR_Abu)

fortestClassesRC_Abu <- normTECounts3m_Abu_final[grepl("RC", normTECounts3m_Abu_final$TEclass),]
nrow(fortestClassesRC_Abu)

fortestClassesSINE_Abu <- normTECounts3m_Abu_final[grepl("SINE", normTECounts3m_Abu_final$TEclass),]
nrow(fortestClassesSINE_Abu)


testDNAAbu <- compare_means(value ~ tissue, data = fortestClassesDNA_Abu, p.adjust.method = "BH")
View(testDNAAbu)

testLINEAbu <- compare_means(value ~ tissue, data = fortestClassesLINE_Abu, p.adjust.method = "BH")
View(testLINEAbu)

testLTRAbu <- compare_means(value ~ tissue, data = fortestClassesLTR_Abu, p.adjust.method = "BH")
View(testLTRAbu)

testRCAbu <- compare_means(value ~ tissue, data = fortestClassesRC_Abu, p.adjust.method = "BH")
View(testRCAbu)

testSINEAbu <- compare_means(value ~ tissue, data = fortestClassesSINE_Abu, p.adjust.method = "BH")
View(testSINEAbu)



### PER FAMILY -> DNA A. burtoni ---- 

normTEcounts_perfamALL_DNAm_Abu <- filter(normTECounts3m_Abu_final, TEclass == "DNA")
#View(normTEcounts_perfamALL_DNAm_Abu)

levels(factor(normTEcounts_perfamALL_DNAm_Abu$TEfam))
normTEcounts_perfamALL_DNAm_Abu$TEfam <- fct_rev(normTEcounts_perfamALL_DNAm_Abu$TEfam)
levels(factor(normTEcounts_perfamALL_DNAm_Abu$TEfam))

plotTEFamDNA_Abu <- ggplot(normTEcounts_perfamALL_DNAm_Abu, aes(x = TEfam, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.5) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,16)) + 
  labs(x = NULL, y = "Log-Normalized counts", title = "DNA transposons") + 
  #labs(title = "Brood size at 25°C", y = "Brood Size", x = NULL) + 
  scale_fill_brewer(palette = "Blues")+
  theme_bw()

plotTEFamDNA_Abu2 <- plotTEFamDNA_Abu + coord_flip()
plotTEFamDNA_Abu2


# Retrotransposons A burtoni ----

normTEcounts_perfamALL_retro1_Abu <- filter(normTECounts3m_Abu_final, TEclass == "LINE")
#View(normTEcounts_perfamALL_retro1_Abu)

normTEcounts_perfamALL_retro2_Abu <- filter(normTECounts3m_Abu_final, TEclass == "LTR")
#View(normTEcounts_perfamALL_retro2_Abu)

normTEcounts_perfamALL_retro3_Abu <- filter(normTECounts3m_Abu_final, TEclass == "Retroposon")
#View(normTEcounts_perfamALL_retro3_Abu)

normTEcounts_perfamALL_retro4_Abu <- filter(normTECounts3m_Abu_final, TEclass == "SINE")
#View(normTEcounts_perfamALL_retro4_Abu)

normTEcounts_perfamALL_retroTEs_Abu <- rbind(normTEcounts_perfamALL_retro1_Abu,normTEcounts_perfamALL_retro2_Abu,
                                             normTEcounts_perfamALL_retro3_Abu,normTEcounts_perfamALL_retro4_Abu)

#View(normTEcounts_perfamALL_retroTEs_Abu)


levels(factor(normTEcounts_perfamALL_retroTEs_Abu$TEfam))
normTEcounts_perfamALL_retroTEs_Abu$TEfam <- fct_rev(normTEcounts_perfamALL_retroTEs_Abu$TEfam)
levels(factor(normTEcounts_perfamALL_retroTEs_Abu$TEfam))

plotTEFamRetro_Abu <- ggplot(normTEcounts_perfamALL_retroTEs_Abu, aes(x = TEfam, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.5) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,16)) + 
  labs(x = NULL, y = "Log-Normalized counts", title = "Retrotransposons") + 
  #labs(title = "Brood size at 25°C", y = "Brood Size", x = NULL) + 
  scale_fill_brewer(palette = "Blues")+
  theme_bw()

plotTEFamRetro_Abu2 <- plotTEFamRetro_Abu + coord_flip()
plotTEFamRetro_Abu2


testRetroFAMsAbu <- compare_means(value ~ tissue, data = normTEcounts_perfamALL_retroTEs_Abu, group.by = "TEfam", p.adjust.method = "BH")
View(testRetroFAMsAbu)

testDNAFAMsAbu <- compare_means(value ~ tissue, data = normTEcounts_perfamALL_DNAm_Abu, group.by = "TEfam", p.adjust.method = "BH")
View(testDNAFAMsAbu)




### P. NYEREREI  ----

Pny_normTEcounts <- read.table(file = "Pnye_NormTEcounts.txt", sep = "\t")
Pny_normTEcounts
#View(Pny_normTEcounts)
nrow(Pny_normTEcounts)

normTECounts2_Pny <- as.data.frame(Pny_normTEcounts)
normTECounts2_Pny$TE <- rownames(normTECounts2_Pny)

rownames(normTECounts2_Pny) <- NULL


normTECounts3_Pny <- normTECounts2_Pny[,c(8,1:7)]
#View(normTECounts3_Pny)
nrow(normTECounts3_Pny)
normTECounts3m_Pny <- melt(normTECounts3_Pny)
#View(normTECounts3m_Pny)
nrow(normTECounts3m_Pny)


normTECounts3m_Pny$TEclass <-  sub(".*:", "", normTECounts3m_Pny$TE)
normTECounts3m_Pny$TEfam <- sub("^[^:]*:([^:]*).*", "\\1", normTECounts3m_Pny$TE)
#View(normTECounts3m_Pny)

normTECounts3m_Pny$tissue <-  sub("^[^_]*_([^_]*).*", "\\1", normTECounts3m_Pny$variable)
#View(normTECounts3m_Pny)

nrow(normTECounts3m_Pny)
normTECounts3m_Pny_final <- filter(normTECounts3m_Pny, !TEclass == "Retroposon" & !TEclass == "SINE?")
nrow(normTECounts3m_Pny_final)


unique(normTECounts3m_Pny_final$TEclass)
unique(normTECounts3m_Pny_final$variable)

normTECounts3m_Pny_final$TEclass <- ordered(normTECounts3m_Pny_final$TEclass,
                                            levels=c("SINE","RC","LTR","LINE","DNA"))



plotTECLASSES_Pny <- ggplot(normTECounts3m_Pny_final, aes(x = TEclass, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,18)) + 
  labs(title = NULL, y = "Log-Normalized counts", x = "TE class") + 
  scale_fill_brewer(palette = "Reds")+
  theme_bw()

plotTECLASSES_Pny2 <- plotTECLASSES_Pny + coord_flip()
plotTECLASSES_Pny2

# stats
fortestClassesDNA_Pny <- normTECounts3m_Pny_final[grepl("DNA", normTECounts3m_Pny_final$TEclass),]
nrow(fortestClassesDNA_Pny)

fortestClassesLINE_Pny <- normTECounts3m_Pny_final[grepl("LINE", normTECounts3m_Pny_final$TEclass),]
nrow(fortestClassesLINE_Pny)

fortestClassesLTR_Pny <- normTECounts3m_Pny_final[grepl("LTR", normTECounts3m_Pny_final$TEclass),]
nrow(fortestClassesLTR_Pny)

fortestClassesRC_Pny <- normTECounts3m_Pny_final[grepl("RC", normTECounts3m_Pny_final$TEclass),]
nrow(fortestClassesRC_Pny)

fortestClassesSINE_Pny <- normTECounts3m_Pny_final[grepl("SINE", normTECounts3m_Pny_final$TEclass),]
nrow(fortestClassesSINE_Pny)


testDNAPny <- compare_means(value ~ tissue, data = fortestClassesDNA_Pny, p.adjust.method = "BH")
View(testDNAPny)

testLINEPny <- compare_means(value ~ tissue, data = fortestClassesLINE_Pny, p.adjust.method = "BH")
View(testLINEPny)

testLTRPny <- compare_means(value ~ tissue, data = fortestClassesLTR_Pny, p.adjust.method = "BH")
View(testLTRPny)

testRCPny <- compare_means(value ~ tissue, data = fortestClassesRC_Pny, p.adjust.method = "BH")
View(testRCPny)

testSINEPny <- compare_means(value ~ tissue, data = fortestClassesSINE_Pny, p.adjust.method = "BH")
View(testSINEPny)



### PER FAMILY -> DNA P. Nyererei ---- 

normTEcounts_perfamALL_DNAm_Pny <- filter(normTECounts3m_Pny_final, TEclass == "DNA")
#View(normTEcounts_perfamALL_DNAm_Pny)

levels(factor(normTEcounts_perfamALL_DNAm_Pny$TEfam))
normTEcounts_perfamALL_DNAm_Pny$TEfam <- fct_rev(normTEcounts_perfamALL_DNAm_Pny$TEfam)
levels(factor(normTEcounts_perfamALL_DNAm_Pny$TEfam))

plotTEFamDNA_Pny <- ggplot(normTEcounts_perfamALL_DNAm_Pny, aes(x = TEfam, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.5) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,16)) + 
  labs(x = NULL, y = "Log-Normalized counts", title = "DNA transposons") + 
  #labs(title = "Brood size at 25°C", y = "Brood Size", x = NULL) + 
  scale_fill_brewer(palette = "Reds")+
  theme_bw()

plotTEFamDNA_Pny2 <- plotTEFamDNA_Pny + coord_flip()
plotTEFamDNA_Pny2



# Retrotransposons P. nyererei ----

normTEcounts_perfamALL_retro1_Pny <- filter(normTECounts3m_Pny_final, TEclass == "LINE")
#View(normTEcounts_perfamALL_retro1_Pny)

normTEcounts_perfamALL_retro2_Pny <- filter(normTECounts3m_Pny_final, TEclass == "LTR")
#View(normTEcounts_perfamALL_retro2_Pny)

normTEcounts_perfamALL_retro3_Pny <- filter(normTECounts3m_Pny_final, TEclass == "Retroposon")
#View(normTEcounts_perfamALL_retro3_Pny)

normTEcounts_perfamALL_retro4_Pny <- filter(normTECounts3m_Pny_final, TEclass == "SINE")
#View(normTEcounts_perfamALL_retro4_Pny)

normTEcounts_perfamALL_retroTEs_Pny <- rbind(normTEcounts_perfamALL_retro1_Pny,normTEcounts_perfamALL_retro2_Pny,
                                             normTEcounts_perfamALL_retro3_Pny,normTEcounts_perfamALL_retro4_Pny)

#View(normTEcounts_perfamALL_retroTEs_Pny)


levels(factor(normTEcounts_perfamALL_retroTEs_Pny$TEfam))
normTEcounts_perfamALL_retroTEs_Pny$TEfam <- fct_rev(normTEcounts_perfamALL_retroTEs_Pny$TEfam)
levels(factor(normTEcounts_perfamALL_retroTEs_Pny$TEfam))

plotTEFamRetro_Pny <- ggplot(normTEcounts_perfamALL_retroTEs_Pny, aes(x = TEfam, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.5) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,16)) + 
  labs(x = NULL, y = "Log-Normalized counts", title = "Retrotransposons") + 
  #labs(title = "Brood size at 25°C", y = "Brood Size", x = NULL) + 
  scale_fill_brewer(palette = "Reds")+
  theme_bw()

plotTEFamRetro_Pny2 <- plotTEFamRetro_Pny + coord_flip()
plotTEFamRetro_Pny2




testRetroFAMsPny <- compare_means(value ~ tissue, data = normTEcounts_perfamALL_retroTEs_Pny, group.by = "TEfam", p.adjust.method = "BH")
View(testRetroFAMsPny)

testDNAFAMsPny <- compare_means(value ~ tissue, data = normTEcounts_perfamALL_DNAm_Pny, group.by = "TEfam", p.adjust.method = "BH")
View(testDNAFAMsPny)





### O. NILOTICUS  ----

Oni_normTEcounts <- read.table(file = "Onil_NormTEcounts.txt", sep = "\t")
Oni_normTEcounts
#View(Oni_normTEcounts)
nrow(Oni_normTEcounts)

normTECounts2_Oni <- as.data.frame(Oni_normTEcounts)
normTECounts2_Oni$TE <- rownames(normTECounts2_Oni)

rownames(normTECounts2_Oni) <- NULL


normTECounts3_Oni <- normTECounts2_Oni[,c(7,1:6)]
#View(normTECounts3_Oni)
nrow(normTECounts3_Oni)
normTECounts3m_Oni <- melt(normTECounts3_Oni)
#View(normTECounts3m_Oni)
nrow(normTECounts3m_Oni)


normTECounts3m_Oni$TEclass <-  sub(".*:", "", normTECounts3m_Oni$TE)
normTECounts3m_Oni$TEfam <- sub("^[^:]*:([^:]*).*", "\\1", normTECounts3m_Oni$TE)
#View(normTECounts3m_Oni)

normTECounts3m_Oni$tissue <-  sub("^[^_]*_([^_]*).*", "\\1", normTECounts3m_Oni$variable)
#View(normTECounts3m_Oni)

nrow(normTECounts3m_Oni)
normTECounts3m_Oni_final <- filter(normTECounts3m_Oni, !TEclass == "Retroposon" & !TEclass == "SINE?")
nrow(normTECounts3m_Oni_final)


unique(normTECounts3m_Oni_final$TEclass)

normTECounts3m_Oni_final$TEclass <- ordered(normTECounts3m_Oni_final$TEclass,
                                            levels=c("SINE","RC","LTR","LINE","DNA"))



plotTECLASSES_Oni <- ggplot(normTECounts3m_Oni_final, aes(x = TEclass, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,18)) + 
  labs(title = NULL, y = "Log-Normalized counts", x = "TE class") + 
  scale_fill_brewer(palette = "Greys")+
  theme_bw()

plotTECLASSES_Oni2 <- plotTECLASSES_Oni + coord_flip()
plotTECLASSES_Oni2

# stats
fortestClassesDNA_Oni <- normTECounts3m_Oni_final[grepl("DNA", normTECounts3m_Oni_final$TEclass),]
nrow(fortestClassesDNA_Oni)

fortestClassesLINE_Oni <- normTECounts3m_Oni_final[grepl("LINE", normTECounts3m_Oni_final$TEclass),]
nrow(fortestClassesLINE_Oni)

fortestClassesLTR_Oni <- normTECounts3m_Oni_final[grepl("LTR", normTECounts3m_Oni_final$TEclass),]
nrow(fortestClassesLTR_Oni)

fortestClassesRC_Oni <- normTECounts3m_Oni_final[grepl("RC", normTECounts3m_Oni_final$TEclass),]
nrow(fortestClassesRC_Oni)

fortestClassesSINE_Oni <- normTECounts3m_Oni_final[grepl("SINE", normTECounts3m_Oni_final$TEclass),]
nrow(fortestClassesSINE_Oni)


testDNAOni <- compare_means(value ~ tissue, data = fortestClassesDNA_Oni, p.adjust.method = "BH")
View(testDNAOni)

testLINEOni <- compare_means(value ~ tissue, data = fortestClassesLINE_Oni, p.adjust.method = "BH")
View(testLINEOni)

testLTROni <- compare_means(value ~ tissue, data = fortestClassesLTR_Oni, p.adjust.method = "BH")
View(testLTROni)

testRCOni <- compare_means(value ~ tissue, data = fortestClassesRC_Oni, p.adjust.method = "BH")
View(testRCOni)

testSINEOni <- compare_means(value ~ tissue, data = fortestClassesSINE_Oni, p.adjust.method = "BH")
View(testSINEOni)



### PER FAMILY -> DNA O. niloticus ---- 

normTEcounts_perfamALL_DNAm_Oni <- filter(normTECounts3m_Oni_final, TEclass == "DNA")
#View(normTEcounts_perfamALL_DNAm_Oni)

levels(factor(normTEcounts_perfamALL_DNAm_Oni$TEfam))
normTEcounts_perfamALL_DNAm_Oni$TEfam <- fct_rev(normTEcounts_perfamALL_DNAm_Oni$TEfam)
levels(factor(normTEcounts_perfamALL_DNAm_Oni$TEfam))

plotTEFamDNA_Oni <- ggplot(normTEcounts_perfamALL_DNAm_Oni, aes(x = TEfam, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.5) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,16)) + 
  labs(x = NULL, y = "Log-Normalized counts", title = "DNA transposons") + 
  #labs(title = "Brood size at 25°C", y = "Brood Size", x = NULL) + 
  scale_fill_brewer(palette = "Greys")+
  theme_bw()

plotTEFamDNA_Oni2 <- plotTEFamDNA_Oni + coord_flip()
plotTEFamDNA_Oni2




# Retrotransposons O niloticus ----

normTEcounts_perfamALL_retro1_Oni <- filter(normTECounts3m_Oni_final, TEclass == "LINE")
#View(normTEcounts_perfamALL_retro1_Oni)

normTEcounts_perfamALL_retro2_Oni <- filter(normTECounts3m_Oni_final, TEclass == "LTR")
#View(normTEcounts_perfamALL_retro2_Oni)

normTEcounts_perfamALL_retro3_Oni <- filter(normTECounts3m_Oni_final, TEclass == "Retroposon")
#View(normTEcounts_perfamALL_retro3_Oni)

normTEcounts_perfamALL_retro4_Oni <- filter(normTECounts3m_Oni_final, TEclass == "SINE")
#View(normTEcounts_perfamALL_retro4_Oni)

normTEcounts_perfamALL_retroTEs_Oni <- rbind(normTEcounts_perfamALL_retro1_Oni,normTEcounts_perfamALL_retro2_Oni,
                                             normTEcounts_perfamALL_retro3_Oni,normTEcounts_perfamALL_retro4_Oni)

#View(normTEcounts_perfamALL_retroTEs_Oni)


levels(factor(normTEcounts_perfamALL_retroTEs_Oni$TEfam))
normTEcounts_perfamALL_retroTEs_Oni$TEfam <- fct_rev(normTEcounts_perfamALL_retroTEs_Oni$TEfam)
levels(factor(normTEcounts_perfamALL_retroTEs_Oni$TEfam))

plotTEFamRetro_Oni <- ggplot(normTEcounts_perfamALL_retroTEs_Oni, aes(x = TEfam, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.5) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,16)) + 
  labs(x = NULL, y = "Log-Normalized counts", title = "Retrotransposons") + 
  #labs(title = "Brood size at 25°C", y = "Brood Size", x = NULL) + 
  scale_fill_brewer(palette = "Greys")+
  theme_bw()

plotTEFamRetro_Oni2 <- plotTEFamRetro_Oni + coord_flip()
plotTEFamRetro_Oni2


testRetroFAMsOni <- compare_means(value ~ tissue, data = normTEcounts_perfamALL_retroTEs_Oni, group.by = "TEfam", p.adjust.method = "BH")
View(testRetroFAMsOni)

testDNAFAMsOni <- compare_means(value ~ tissue, data = normTEcounts_perfamALL_DNAm_Oni, group.by = "TEfam", p.adjust.method = "BH")
View(testDNAFAMsOni)



##########################################################################################
########  Final plots together ----  #####################################################
##########################################################################################

# All classes together

plotTECLASSES_Aca2 + plotTECLASSES_Abu2 + plotTECLASSES_Pny2 + plotTECLASSES_Oni2 +
  plot_layout(ncol = 4)

plotTECLASSES_AcaPio3 <- plotTECLASSES_AcaPio2 + theme(legend.position = "none")
plotTECLASSES_Aca3 <- plotTECLASSES_Aca2 + theme(legend.position = "none")
plotTECLASSES_Abu3 <- plotTECLASSES_Abu2 + theme(legend.position = "none")
plotTECLASSES_Pny3 <- plotTECLASSES_Pny2 + theme(legend.position = "none")
plotTECLASSES_Oni3 <- plotTECLASSES_Oni2 + theme(legend.position = "none")



plotTEFamDNA_AcaPio3 + plotTEFamDNA_Aca3 + plotTEFamDNA_Abu3 + plotTEFamDNA_Pny3 + plotTEFamDNA_Oni3  +
  plotTEFamRetro_AcaPio3 + plotTEFamRetro_Aca3 + plotTEFamRetro_Abu3 + plotTEFamRetro_Pny3 + plotTEFamRetro_Oni3 + 
  plot_layout(nrow = 2,ncol = 5)


# 13.85 x 8.55



##################################################################################################
######## In Acal Pio dataset, check expression of Active vs inactive families   ----  ############
##################################################################################################


normTECountsm_final_AcalGondas_forFilt <- normTECounts3m_AcaPio_final
normTECountsm_final_AcalGondas_forFilt$TEsubfam <- sub(":.*", "", normTECountsm_final_AcalGondas_forFilt$TE)

normTECountsm_final_AcalGonads_gypsy <- filter(normTECountsm_final_AcalGondas_forFilt, TEfam == "Gypsy")
#View(normTECountsm_final_AcalGonads_gypsy)

normTECountsm_final_AcalGonads_gypsy$NewTEsubfam <- sub("\\_\\[[A-Z]+\\]$", "", normTECountsm_final_AcalGonads_gypsy$TEsubfam)
#View(normTECountsm_final_AcalGonads_gypsy)


gypsy_active <- read.table(file = "gypsy_families_potentiallyActiveNEW.txt", header = FALSE)


gypsy_active_reads <- filter(normTECountsm_final_AcalGonads_gypsy, NewTEsubfam %in% gypsy_active$V1)


normTECountsm_final_AcalGonads_gypsy <- anti_join(normTECountsm_final_AcalGonads_gypsy,gypsy_active_reads)


gypsy_active_reads$MobStatus <- "Active"
normTECountsm_final_AcalGonads_gypsy$MobStatus <- "Inactive"
#View(normTECountsm_final_AcalGonads_gypsy)
#View(gypsy_active_reads)

normTECountsm_final_AcalGonads_gypsy_final <- bind_rows(gypsy_active_reads,normTECountsm_final_AcalGonads_gypsy)
#View(normTECountsm_final_AcalGonads_gypsy_final)

plotTE_Gypsy_Active <- ggplot(normTECountsm_final_AcalGonads_gypsy_final, aes(x = tissue, y = value, fill = MobStatus)) +
  geom_boxplot(aes(fill = MobStatus), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(-3,14)) + 
  labs(title = "Gypsy", y = "Log-Normalized counts", x = NULL) + 
  scale_fill_brewer(palette = "Accent") +
  theme_bw()

plotTE_Gypsy_Active 




Acal_gypsy_ov <- normTECountsm_final_AcalGonads_gypsy_final[grepl("ov", normTECountsm_final_AcalGonads_gypsy_final$tissue),]
nrow(Acal_gypsy_ov)

Acal_gypsy_tes <- normTECountsm_final_AcalGonads_gypsy_final[grepl("tes", normTECountsm_final_AcalGonads_gypsy_final$tissue),]
nrow(Acal_gypsy_tes)





testGypsy_mobStatus_tes <- compare_means(value ~ MobStatus, data = Acal_gypsy_tes, p.adjust.method = "BH")
View(testGypsy_mobStatus_tes)

testGypsy_mobStatus_ov <- compare_means(value ~ MobStatus, data = Acal_gypsy_ov, p.adjust.method = "BH")
View(testGypsy_mobStatus_ov)




normTECountsm_final_AcalGonads_ERV1 <- filter(normTECountsm_final_AcalGondas_forFilt, TEfam == "ERV1")
#View(normTECountsm_final_AcalGonads_ERV1)

normTECountsm_final_AcalGonads_ERV1$NewTEsubfam <- sub("\\_\\[[A-Z]+\\]$", "", normTECountsm_final_AcalGonads_ERV1$TEsubfam)
#View(normTECountsm_final_AcalGonads_ERV1)


ERV1_active <- read.table(file = "ERV_families_potentiallyActiveNEW.txt", header = FALSE)


ERV1_active_reads <- filter(normTECountsm_final_AcalGonads_ERV1, NewTEsubfam %in% ERV1_active$V1)


normTECountsm_final_AcalGonads_ERV1 <- anti_join(normTECountsm_final_AcalGonads_ERV1,ERV1_active_reads)


ERV1_active_reads$MobStatus <- "Active"
normTECountsm_final_AcalGonads_ERV1$MobStatus <- "Inactive"
#View(normTECountsm_final_AcalGonads_ERV1)
#View(ERV1_active_reads)

normTECountsm_final_AcalGonads_ERV1_final <- bind_rows(ERV1_active_reads,normTECountsm_final_AcalGonads_ERV1)


plotTE_ERV1_Active <- ggplot(normTECountsm_final_AcalGonads_ERV1_final, aes(x = tissue, y = value, fill = MobStatus)) +
  geom_boxplot(aes(fill = MobStatus), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(0,17)) + 
  labs(title = "ERV1", y = "Log-Normalized counts", x = NULL) + 
  scale_fill_brewer(palette = "Accent") +
  theme_bw()

plotTE_ERV1_Active 




Acal_ERV1_ov <- normTECountsm_final_AcalGonads_ERV1_final[grepl("ov", normTECountsm_final_AcalGonads_ERV1_final$tissue),]
nrow(Acal_ERV1_ov)

Acal_ERV1_tes <- normTECountsm_final_AcalGonads_ERV1_final[grepl("tes", normTECountsm_final_AcalGonads_ERV1_final$tissue),]
nrow(Acal_ERV1_tes)





testERV1_mobStatus_tes <- compare_means(value ~ MobStatus, data = Acal_ERV1_tes, p.adjust.method = "BH")
View(testERV1_mobStatus_tes)

testERV1_mobStatus_ov <- compare_means(value ~ MobStatus, data = Acal_ERV1_ov, p.adjust.method = "BH")
View(testERV1_mobStatus_ov)





#View(normTECountsm_final_AcalGondas_forFilt)

normTECountsm_final_AcalGonads_CMC <- filter(normTECountsm_final_AcalGondas_forFilt, TEfam == "CMC-EnSpm")
#View(normTECountsm_final_AcalGonads_CMC)

normTECountsm_final_AcalGonads_CMC$NewTEsubfam <- sub("\\_\\[[A-Z]+\\]$", "", normTECountsm_final_AcalGonads_CMC$TEsubfam)
#View(normTECountsm_final_AcalGonads_CMC)


CMC_active <- read.table(file = "CMC_families_potentiallyActiveNEW.txt", header = FALSE)


CMC_active_reads <- filter(normTECountsm_final_AcalGonads_CMC, NewTEsubfam %in% CMC_active$V1)


normTECountsm_final_AcalGonads_CMC <- anti_join(normTECountsm_final_AcalGonads_CMC,CMC_active_reads)


CMC_active_reads$MobStatus <- "Active"
normTECountsm_final_AcalGonads_CMC$MobStatus <- "Inactive"
#View(normTECountsm_final_AcalGonads_CMC)
#View(CMC_active_reads)

normTECountsm_final_AcalGonads_CMC_final <- bind_rows(CMC_active_reads,normTECountsm_final_AcalGonads_CMC)
#View(normTECountsm_final_AcalGonads_CMC_final)


plotTE_CMC_Active <- ggplot(normTECountsm_final_AcalGonads_CMC_final, aes(x = tissue, y = value, fill = MobStatus)) +
  geom_boxplot(aes(fill = MobStatus), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(0,12)) + 
  labs(title = "CMC-EnSpm", y = "Log-Normalized counts", x = NULL) + 
  scale_fill_brewer(palette = "Accent") +
  theme_bw()

plotTE_CMC_Active 




Acal_CMC_ov <- normTECountsm_final_AcalGonads_CMC_final[grepl("ov", normTECountsm_final_AcalGonads_CMC_final$tissue),]
nrow(Acal_CMC_ov)

Acal_CMC_tes <- normTECountsm_final_AcalGonads_CMC_final[grepl("tes", normTECountsm_final_AcalGonads_CMC_final$tissue),]
nrow(Acal_CMC_tes)





testCMC_mobStatus_tes <- compare_means(value ~ MobStatus, data = Acal_CMC_tes, p.adjust.method = "BH")
View(testCMC_mobStatus_tes)

testCMC_mobStatus_ov <- compare_means(value ~ MobStatus, data = Acal_CMC_ov, p.adjust.method = "BH")
View(testCMC_mobStatus_ov)





normTECountsm_final_AcalGonads_Rex <- filter(normTECountsm_final_AcalGondas_forFilt, TEfam == "Rex-Babar")
#View(normTECountsm_final_AcalGonads_Rex)

normTECountsm_final_AcalGonads_Rex$NewTEsubfam <- sub("\\_\\[[A-Z]+\\]$", "", normTECountsm_final_AcalGonads_Rex$TEsubfam)
#View(normTECountsm_final_AcalGonads_Rex)


Rex_active <- read.table(file = "Rex_families_potentiallyActiveNEW.txt", header = FALSE)


Rex_active_reads <- filter(normTECountsm_final_AcalGonads_Rex, NewTEsubfam %in% Rex_active$V1)
## NOTHING ON REFERENCE. DO NOT GO ON WITH THIS






normTECountsm_final_AcalGonads_Maverick <- filter(normTECountsm_final_AcalGondas_forFilt, TEfam == "Maverick")
#View(normTECountsm_final_AcalGonads_Maverick)

normTECountsm_final_AcalGonads_Maverick$NewTEsubfam <- sub("\\_\\[[A-Z]+\\]$", "", normTECountsm_final_AcalGonads_Maverick$TEsubfam)
#View(normTECountsm_final_AcalGonads_Maverick)


Maverick_active <- read.table(file = "Maverick_families_potentiallyActiveNEW.txt", header = FALSE)


Maverick_active_reads <- filter(normTECountsm_final_AcalGonads_Maverick, NewTEsubfam %in% Maverick_active$V1)


normTECountsm_final_AcalGonads_Maverick <- anti_join(normTECountsm_final_AcalGonads_Maverick,Maverick_active_reads)


Maverick_active_reads$MobStatus <- "Active"
normTECountsm_final_AcalGonads_Maverick$MobStatus <- "Inactive"
#View(normTECountsm_final_AcalGonads_Maverick)
#View(Maverick_active_reads)

normTECountsm_final_AcalGonads_Maverick_final <- bind_rows(Maverick_active_reads,normTECountsm_final_AcalGonads_Maverick)
#View(normTECountsm_final_AcalGonads_Maverick_final)


plotTE_Maverick_Active <- ggplot(normTECountsm_final_AcalGonads_Maverick_final, aes(x = tissue, y = value, fill = MobStatus)) +
  geom_boxplot(aes(fill = MobStatus), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(0,13)) + 
  labs(title = "Maverick", y = "Log-Normalized counts", x = NULL) + 
  scale_fill_brewer(palette = "Accent") +
  theme_bw()

plotTE_Maverick_Active 




Acal_Maverick_ov <- normTECountsm_final_AcalGonads_Maverick_final[grepl("ov", normTECountsm_final_AcalGonads_Maverick_final$tissue),]
nrow(Acal_Maverick_ov)

Acal_Maverick_tes <- normTECountsm_final_AcalGonads_Maverick_final[grepl("tes", normTECountsm_final_AcalGonads_Maverick_final$tissue),]
nrow(Acal_Maverick_tes)



testMaverick_mobStatus_tes <- compare_means(value ~ MobStatus, data = Acal_Maverick_tes, p.adjust.method = "BH")
View(testMaverick_mobStatus_tes)

testMaverick_mobStatus_ov <- compare_means(value ~ MobStatus, data = Acal_Maverick_ov, p.adjust.method = "BH")
View(testMaverick_mobStatus_ov)



### split mavericks into 3 families

normTECountsm_final_AcalGonads_Maverick_2 <- filter(normTECountsm_final_AcalGondas_forFilt, TEfam == "Maverick")
#View(normTECountsm_final_AcalGonads_Maverick_2)

normTECountsm_final_AcalGonads_Maverick_2$NewTEsubfam <- sub("\\_\\[[A-Z]+\\]$", "", normTECountsm_final_AcalGonads_Maverick_2$TEsubfam)
#View(normTECountsm_final_AcalGonads_Maverick_2)



plotTE_Maverick_Active_split <- ggplot(normTECountsm_final_AcalGonads_Maverick_2, aes(x = tissue, y = value, fill = NewTEsubfam)) +
  geom_boxplot(aes(fill = NewTEsubfam), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(0,13)) + 
  labs(title = "Maverick", y = "Log-Normalized counts", x = NULL) + 
  scale_fill_brewer(palette = "Accent") +
  theme_bw()

plotTE_Maverick_Active_split 


Acal_Maverick_ov_split <- normTECountsm_final_AcalGonads_Maverick_2[grepl("ov", normTECountsm_final_AcalGonads_Maverick_2$tissue),]
nrow(Acal_Maverick_ov_split)

Acal_Maverick_tes_split <- normTECountsm_final_AcalGonads_Maverick_2[grepl("tes", normTECountsm_final_AcalGonads_Maverick_2$tissue),]
nrow(Acal_Maverick_tes_split)



testMaverick_mobStatus_tes_split <- compare_means(value ~ NewTEsubfam, data = Acal_Maverick_tes_split, p.adjust.method = "BH")
View(testMaverick_mobStatus_tes_split)

testMaverick_mobStatus_ov_split <- compare_means(value ~ NewTEsubfam, data = Acal_Maverick_ov_split, p.adjust.method = "BH")
View(testMaverick_mobStatus_ov_split)

### decided not to use that



normTECountsm_final_AcalGonads_TcMarTc2 <- filter(normTECountsm_final_AcalGondas_forFilt, TEfam == "TcMar-Tc2")
#View(normTECountsm_final_AcalGonads_TcMarTc2)

normTECountsm_final_AcalGonads_TcMarTc2$NewTEsubfam <- sub("\\_\\[[A-Z]+\\]$", "", normTECountsm_final_AcalGonads_TcMarTc2$TEsubfam)
#View(normTECountsm_final_AcalGonads_TcMarTc2)


TcMarTc2_active <- read.table(file = "TcMar_families_potentiallyActiveNEW.txt", header = FALSE)


TcMarTc2_active_reads <- filter(normTECountsm_final_AcalGonads_TcMarTc2, NewTEsubfam %in% TcMarTc2_active$V1)


normTECountsm_final_AcalGonads_TcMarTc2 <- anti_join(normTECountsm_final_AcalGonads_TcMarTc2,TcMarTc2_active_reads)


TcMarTc2_active_reads$MobStatus <- "Active"
normTECountsm_final_AcalGonads_TcMarTc2$MobStatus <- "Inactive"
#View(normTECountsm_final_AcalGonads_TcMarTc2)
#View(TcMarTc2_active_reads)

normTECountsm_final_AcalGonads_TcMarTc2_final <- bind_rows(TcMarTc2_active_reads,normTECountsm_final_AcalGonads_TcMarTc2)
#View(normTECountsm_final_AcalGonads_TcMarTc2_final)


plotTE_TcMarTc2_Active <- ggplot(normTECountsm_final_AcalGonads_TcMarTc2_final, aes(x = tissue, y = value, fill = MobStatus)) +
  geom_boxplot(aes(fill = MobStatus), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(0,14)) + 
  labs(title = "TcMar-Tc2", y = "Log-Normalized counts", x = NULL) + 
  scale_fill_brewer(palette = "Accent") +
  theme_bw()

plotTE_TcMarTc2_Active 




Acal_TcMarTc2_ov <- normTECountsm_final_AcalGonads_TcMarTc2_final[grepl("ov", normTECountsm_final_AcalGonads_TcMarTc2_final$tissue),]
nrow(Acal_TcMarTc2_ov)

Acal_TcMarTc2_tes <- normTECountsm_final_AcalGonads_TcMarTc2_final[grepl("tes", normTECountsm_final_AcalGonads_TcMarTc2_final$tissue),]
nrow(Acal_TcMarTc2_tes)





testTcMarTc2_mobStatus_tes <- compare_means(value ~ MobStatus, data = Acal_TcMarTc2_tes, p.adjust.method = "BH")
View(testTcMarTc2_mobStatus_tes)

testTcMarTc2_mobStatus_ov <- compare_means(value ~ MobStatus, data = Acal_TcMarTc2_ov, p.adjust.method = "BH")
View(testTcMarTc2_mobStatus_ov)







normTECountsm_final_AcalGonads_PiggyBac5 <- filter(normTECountsm_final_AcalGondas_forFilt, TEfam == "PiggyBac")
#View(normTECountsm_final_AcalGonads_PiggyBac5)

normTECountsm_final_AcalGonads_PiggyBac5$NewTEsubfam <- sub("\\_\\[[A-Z]+\\]$", "", normTECountsm_final_AcalGonads_PiggyBac5$TEsubfam)
#View(normTECountsm_final_AcalGonads_PiggyBac5)


PiggyBac5_active <- read.table(file = "Piggy_families_potentiallyActiveNEW.txt", header = FALSE)


PiggyBac5_active_reads <- filter(normTECountsm_final_AcalGonads_PiggyBac5, NewTEsubfam %in% PiggyBac5_active$V1)


normTECountsm_final_AcalGonads_PiggyBac5 <- anti_join(normTECountsm_final_AcalGonads_PiggyBac5,PiggyBac5_active_reads)


PiggyBac5_active_reads$MobStatus <- "Active"
normTECountsm_final_AcalGonads_PiggyBac5$MobStatus <- "Inactive"
#View(normTECountsm_final_AcalGonads_PiggyBac5)
#View(PiggyBac5_active_reads)

normTECountsm_final_AcalGonads_PiggyBac5_final <- bind_rows(PiggyBac5_active_reads,normTECountsm_final_AcalGonads_PiggyBac5)
#View(normTECountsm_final_AcalGonads_PiggyBac5_final)


plotTE_PiggyBac5_Active <- ggplot(normTECountsm_final_AcalGonads_PiggyBac5_final, aes(x = tissue, y = value, fill = MobStatus)) +
  geom_boxplot(aes(fill = MobStatus), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(0,11)) + 
  labs(title = "PiggyBac", y = "Log-Normalized counts", x = NULL) + 
  scale_fill_brewer(palette = "Accent") +
  theme_bw()

plotTE_PiggyBac5_Active 




Acal_PiggyBac5_ov <- normTECountsm_final_AcalGonads_PiggyBac5_final[grepl("ov", normTECountsm_final_AcalGonads_PiggyBac5_final$tissue),]
nrow(Acal_PiggyBac5_ov)

Acal_PiggyBac5_tes <- normTECountsm_final_AcalGonads_PiggyBac5_final[grepl("tes", normTECountsm_final_AcalGonads_PiggyBac5_final$tissue),]
nrow(Acal_PiggyBac5_tes)





testPiggyBac5_mobStatus_tes <- compare_means(value ~ MobStatus, data = Acal_PiggyBac5_tes, p.adjust.method = "BH")
View(testPiggyBac5_mobStatus_tes)

testPiggyBac5_mobStatus_ov <- compare_means(value ~ MobStatus, data = Acal_PiggyBac5_ov, p.adjust.method = "BH")
View(testPiggyBac5_mobStatus_ov)





normTECountsm_final_AcalGonads_DIRS1 <- filter(normTECountsm_final_AcalGondas_forFilt, TEfam == "DIRS")
#View(normTECountsm_final_AcalGonads_DIRS1)

normTECountsm_final_AcalGonads_DIRS1$NewTEsubfam <- sub("\\_\\[[A-Z]+\\]$", "", normTECountsm_final_AcalGonads_DIRS1$TEsubfam)
#View(normTECountsm_final_AcalGonads_DIRS1)


DIRS1_active <- read.table(file = "DIRS_families_potentiallyActiveNEW.txt", header = FALSE)


DIRS1_active_reads <- filter(normTECountsm_final_AcalGonads_DIRS1, NewTEsubfam %in% DIRS1_active$V1)


normTECountsm_final_AcalGonads_DIRS1 <- anti_join(normTECountsm_final_AcalGonads_DIRS1,DIRS1_active_reads)


DIRS1_active_reads$MobStatus <- "Active"
normTECountsm_final_AcalGonads_DIRS1$MobStatus <- "Inactive"
#View(normTECountsm_final_AcalGonads_DIRS1)
#View(DIRS1_active_reads)

normTECountsm_final_AcalGonads_DIRS1_final <- bind_rows(DIRS1_active_reads,normTECountsm_final_AcalGonads_DIRS1)
#View(normTECountsm_final_AcalGonads_DIRS1_final)


plotTE_DIRS1_Active <- ggplot(normTECountsm_final_AcalGonads_DIRS1_final, aes(x = tissue, y = value, fill = MobStatus)) +
  geom_boxplot(aes(fill = MobStatus), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(0,10)) + 
  labs(title = "DIRS", y = "Log-Normalized counts", x = NULL) + 
  scale_fill_brewer(palette = "Accent") +
  theme_bw()

plotTE_DIRS1_Active 




Acal_DIRS1_ov <- normTECountsm_final_AcalGonads_DIRS1_final[grepl("ov", normTECountsm_final_AcalGonads_DIRS1_final$tissue),]
nrow(Acal_DIRS1_ov)

Acal_DIRS1_tes <- normTECountsm_final_AcalGonads_DIRS1_final[grepl("tes", normTECountsm_final_AcalGonads_DIRS1_final$tissue),]
nrow(Acal_DIRS1_tes)





testDIRS1_mobStatus_tes <- compare_means(value ~ MobStatus, data = Acal_DIRS1_tes, p.adjust.method = "BH")
View(testDIRS1_mobStatus_tes)

testDIRS1_mobStatus_ov <- compare_means(value ~ MobStatus, data = Acal_DIRS1_ov, p.adjust.method = "BH")
View(testDIRS1_mobStatus_ov)





normTECountsm_final_AcalGonads_hATTiP100 <- filter(normTECountsm_final_AcalGondas_forFilt, TEfam == "hAT-Tip100")
#View(normTECountsm_final_AcalGonads_hATTiP100)

normTECountsm_final_AcalGonads_hATTiP100$NewTEsubfam <- sub("\\_\\[[A-Z]+\\]$", "", normTECountsm_final_AcalGonads_hATTiP100$TEsubfam)
#View(normTECountsm_final_AcalGonads_hATTiP100)


hATTiP100_active <- read.table(file = "hATtip100_families_potentiallyActiveNEW.txt", header = FALSE)


hATTiP100_active_reads <- filter(normTECountsm_final_AcalGonads_hATTiP100, NewTEsubfam %in% hATTiP100_active$V1)


normTECountsm_final_AcalGonads_hATTiP100 <- anti_join(normTECountsm_final_AcalGonads_hATTiP100,hATTiP100_active_reads)


hATTiP100_active_reads$MobStatus <- "Active"
normTECountsm_final_AcalGonads_hATTiP100$MobStatus <- "Inactive"
#View(normTECountsm_final_AcalGonads_hATTiP100)
#View(hATTiP100_active_reads)

normTECountsm_final_AcalGonads_hATTiP100_final <- bind_rows(hATTiP100_active_reads,normTECountsm_final_AcalGonads_hATTiP100)
#View(normTECountsm_final_AcalGonads_hATTiP100_final)


plotTE_hATTiP100_Active <- ggplot(normTECountsm_final_AcalGonads_hATTiP100_final, aes(x = tissue, y = value, fill = MobStatus)) +
  geom_boxplot(aes(fill = MobStatus), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(0,11)) + 
  labs(title = "hAT-Tip100", y = "Log-Normalized counts", x = NULL) + 
  scale_fill_brewer(palette = "Accent") +
  theme_bw()

plotTE_hATTiP100_Active 




Acal_hATTiP100_ov <- normTECountsm_final_AcalGonads_hATTiP100_final[grepl("ov", normTECountsm_final_AcalGonads_hATTiP100_final$tissue),]
nrow(Acal_hATTiP100_ov)

Acal_hATTiP100_tes <- normTECountsm_final_AcalGonads_hATTiP100_final[grepl("tes", normTECountsm_final_AcalGonads_hATTiP100_final$tissue),]
nrow(Acal_hATTiP100_tes)



testhATTiP100_mobStatus_tes <- compare_means(value ~ MobStatus, data = Acal_hATTiP100_tes, p.adjust.method = "BH")
View(testhATTiP100_mobStatus_tes)

testhATTiP100_mobStatus_ov <- compare_means(value ~ MobStatus, data = Acal_hATTiP100_ov, p.adjust.method = "BH")
View(testhATTiP100_mobStatus_ov)





normTECountsm_final_AcalGonads_hATAc <- filter(normTECountsm_final_AcalGondas_forFilt, TEfam == "hAT-Ac")
#View(normTECountsm_final_AcalGonads_hATAc)

normTECountsm_final_AcalGonads_hATAc$NewTEsubfam <- sub("\\_\\[[A-Z]+\\]$", "", normTECountsm_final_AcalGonads_hATAc$TEsubfam)
#View(normTECountsm_final_AcalGonads_hATAc)


hATAc_active <- read.table(file = "hATAc_families_potentiallyActiveNEW.txt", header = FALSE)


hATAc_active_reads <- filter(normTECountsm_final_AcalGonads_hATAc, NewTEsubfam %in% hATAc_active$V1)


normTECountsm_final_AcalGonads_hATAc <- anti_join(normTECountsm_final_AcalGonads_hATAc,hATAc_active_reads)


hATAc_active_reads$MobStatus <- "Active"
normTECountsm_final_AcalGonads_hATAc$MobStatus <- "Inactive"
#View(normTECountsm_final_AcalGonads_hATAc)
#View(hATAc_active_reads)

normTECountsm_final_AcalGonads_hATAc_final <- bind_rows(hATAc_active_reads,normTECountsm_final_AcalGonads_hATAc)
#View(normTECountsm_final_AcalGonads_hATAc_final)


plotTE_hATAc_Active <- ggplot(normTECountsm_final_AcalGonads_hATAc_final, aes(x = tissue, y = value, fill = MobStatus)) +
  geom_boxplot(aes(fill = MobStatus), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(0,14)) + 
  labs(title = "hAT-Ac", y = "Log-Normalized counts", x = NULL) + 
  scale_fill_brewer(palette = "Accent") +
  theme_bw()

plotTE_hATAc_Active 




Acal_hATAc_ov <- normTECountsm_final_AcalGonads_hATAc_final[grepl("ov", normTECountsm_final_AcalGonads_hATAc_final$tissue),]
nrow(Acal_hATAc_ov)

Acal_hATAc_tes <- normTECountsm_final_AcalGonads_hATAc_final[grepl("tes", normTECountsm_final_AcalGonads_hATAc_final$tissue),]
nrow(Acal_hATAc_tes)





testhATAc_mobStatus_tes <- compare_means(value ~ MobStatus, data = Acal_hATAc_tes, p.adjust.method = "BH")
View(testhATAc_mobStatus_tes)

testhATAc_mobStatus_ov <- compare_means(value ~ MobStatus, data = Acal_hATAc_ov, p.adjust.method = "BH")
View(testhATAc_mobStatus_ov)



### plot all together


## first remove the legengs:


plotTE_TcMarTc2_Active2 <- plotTE_TcMarTc2_Active + theme(legend.position = "none")
plotTE_hATTiP100_Active2 <- plotTE_hATTiP100_Active + theme(legend.position = "none")
plotTE_PiggyBac5_Active2 <- plotTE_PiggyBac5_Active + theme(legend.position = "none")
plotTE_DIRS1_Active2 <- plotTE_DIRS1_Active + theme(legend.position = "none")
plotTE_CMC_Active2 <- plotTE_CMC_Active + theme(legend.position = "none")
plotTE_Maverick_Active2 <- plotTE_Maverick_Active + theme(legend.position = "none")
plotTE_hATAc_Active2 <- plotTE_hATAc_Active + theme(legend.position = "none")
plotTE_Gypsy_Active2 <- plotTE_Gypsy_Active + theme(legend.position = "none")
plotTE_ERV1_Active2 <- plotTE_ERV1_Active + theme(legend.position = "none")


wrap_plots(plotTE_DIRS1_Active2, plotTE_hATTiP100_Active2, plotTE_PiggyBac5_Active2,plotTE_Gypsy_Active2,plotTE_hATAc_Active2,
           plotTE_CMC_Active2,plotTE_ERV1_Active2,plotTE_Maverick_Active2,plotTE_TcMarTc2_Active2,plot_spacer(), nrow = 2, ncol = 5)

## 10.70 x 4.66




