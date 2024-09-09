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



# load files with filtered TE counts. Filtering for more than 10 counts in more than 2 samples.

# for A. burtoni

Abur_expressedTEs <- read.table("Abur_filtered_notNormTEcounts.txt", sep = "\t")
#View(Abur_expressedTEs)
Abur_expressedTEs$TE <- rownames(Abur_expressedTEs)
rownames(Abur_expressedTEs) <- NULL
#View(Abur_expressedTEs)

Abur_expressedTEs$TEfam <- sub("^[^:]*:([^:]*).*", "\\1", Abur_expressedTEs$TE)
Abur_expressedTEs$TEclass <-  sub(".*:", "", Abur_expressedTEs$TE)

#View(Abur_expressedTEs)

Abur_expressedTEs <- filter(Abur_expressedTEs, !TEclass == "Retroposon" & !TEclass == "SINE?")
#View(Abur_expressedTEs)
nrow(Abur_expressedTEs)


colourCount = length(unique(Abur_expressedTEs$TEfam))
getPalette = colorRampPalette(brewer.pal(9, "Blues"))



ABU_EXPRESSEDTES_bp <- ggplot(Abur_expressedTEs, aes(x = TEclass, fill = TEfam)) +
  geom_bar(color = "black") + 
  scale_fill_manual(values = getPalette(colourCount)) +
  labs(y = NULL, x = NULL) + 
  theme_bw()
ABU_EXPRESSEDTES_bp <- ABU_EXPRESSEDTES_bp + guides(fill=guide_legend(title = "TE family"))  +
  theme(axis.text.x = element_text(angle=-45, vjust = 0.5, hjust = 0.1))
ABU_EXPRESSEDTES_bp

abu_dnaexp <- dplyr::filter(Abur_expressedTEs, grepl("DNA", TEclass, ignore.case=TRUE))
nrow(abu_dnaexp) # 423 DNA subfamilies expressed

abu_lineexp <- dplyr::filter(Abur_expressedTEs, grepl("line", TEclass, ignore.case=TRUE))
nrow(abu_lineexp) # 213 LINE subfamilies expressed

abu_ltrexp <- dplyr::filter(Abur_expressedTEs, grepl("LTR", TEclass, ignore.case=TRUE))
nrow(abu_ltrexp) # 39 LTR subfamilies expressed

abu_helitronexp <- dplyr::filter(Abur_expressedTEs, grepl("helitron", TEfam, ignore.case=TRUE))
nrow(abu_helitronexp) # 4 helitron subfamilies expressed

abu_sineexp <- dplyr::filter(Abur_expressedTEs, grepl("SINE", TEclass, ignore.case=TRUE))
nrow(abu_sineexp) # 10 SINE subfamilies expressed

# for A. calliptera


Acal_expressedTEs <- read.table("Acal_filtered_notNormTEcounts.txt", sep = "\t")
#View(Acal_expressedTEs)
Acal_expressedTEs$TE <- rownames(Acal_expressedTEs)
rownames(Acal_expressedTEs) <- NULL
#View(Acal_expressedTEs)

Acal_expressedTEs$TEfam <- sub("^[^:]*:([^:]*).*", "\\1", Acal_expressedTEs$TE)
Acal_expressedTEs$TEclass <-  sub(".*:", "", Acal_expressedTEs$TE)

#View(Acal_expressedTEs)

Acal_expressedTEs <- filter(Acal_expressedTEs, !TEclass == "Retroposon" & !TEclass == "SINE?")
#View(Acal_expressedTEs)
nrow(Acal_expressedTEs)


colourCount = length(unique(Acal_expressedTEs$TEfam))
getPalette = colorRampPalette(brewer.pal(9, "Oranges"))



ACA_EXPRESSEDTES_bp <- ggplot(Acal_expressedTEs, aes(x = TEclass, fill = TEfam)) +
  geom_bar(color = "black") + 
  scale_fill_manual(values = getPalette(colourCount)) +
  labs(y = NULL, x = NULL) + 
  theme_bw()
ACA_EXPRESSEDTES_bp <- ACA_EXPRESSEDTES_bp + guides(fill=guide_legend(title = "TE family"))  +
  theme(axis.text.x = element_text(angle=-45, vjust = 0.5, hjust = 0.1))
ACA_EXPRESSEDTES_bp

aca_dnaexp <- dplyr::filter(Acal_expressedTEs, grepl("DNA", TEclass, ignore.case=TRUE))
nrow(aca_dnaexp) # 396 DNA subfamilies expressed

aca_lineexp <- dplyr::filter(Acal_expressedTEs, grepl("line", TEclass, ignore.case=TRUE))
nrow(aca_lineexp) # 178 LINE subfamilies expressed

aca_ltrexp <- dplyr::filter(Acal_expressedTEs, grepl("LTR", TEclass, ignore.case=TRUE))
nrow(aca_ltrexp) # 69 LTR subfamilies expressed

aca_helitronexp <- dplyr::filter(Acal_expressedTEs, grepl("helitron", TEfam, ignore.case=TRUE))
nrow(aca_helitronexp) # 1 helitron subfamilies expressed

aca_sineexp <- dplyr::filter(Acal_expressedTEs, grepl("SINE", TEclass, ignore.case=TRUE))
nrow(aca_sineexp) # 12 SINE subfamilies expressed



# for A. calliptera with curated annotation


Acal_Pio_expressedTEs <- read.table("Acal_filtered_notNormTEcounts_Pio.txt", sep = "\t")
#View(Acal_Pio_expressedTEs)
Acal_Pio_expressedTEs$TE <- rownames(Acal_Pio_expressedTEs)
rownames(Acal_Pio_expressedTEs) <- NULL
#View(Acal_Pio_expressedTEs)

Acal_Pio_expressedTEs$TEfam <- sub("^[^:]*:([^:]*).*", "\\1", Acal_Pio_expressedTEs$TE)
Acal_Pio_expressedTEs$TEclass <-  sub(".*:", "", Acal_Pio_expressedTEs$TE)

#View(Acal_Pio_expressedTEs)

Acal_Pio_expressedTEs <- filter(Acal_Pio_expressedTEs, !TEclass == "Retroposon" & !TEclass == "SINE?")
#View(Acal_Pio_expressedTEs)
nrow(Acal_Pio_expressedTEs)


colourCount = length(unique(Acal_Pio_expressedTEs$TEfam))
getPalette = colorRampPalette(brewer.pal(9, "Oranges"))



ACA_Pio_EXPRESSEDTES_bp <- ggplot(Acal_Pio_expressedTEs, aes(x = TEclass, fill = TEfam)) +
  geom_bar(color = "black") + 
  scale_fill_manual(values = getPalette(colourCount)) +
  labs(y = NULL, x = NULL) + 
  theme_bw()
ACA_Pio_EXPRESSEDTES_bp <- ACA_EXPRESSEDTES_bp + guides(fill=guide_legend(title = "TE family"))  +
  theme(axis.text.x = element_text(angle=-45, vjust = 0.5, hjust = 0.1))
ACA_Pio_EXPRESSEDTES_bp

aca_Pio_dnaexp <- dplyr::filter(Acal_Pio_expressedTEs, grepl("DNA", TEclass, ignore.case=TRUE))
nrow(aca_Pio_dnaexp) # 132 DNA subfamilies expressed

aca_Pio_lineexp <- dplyr::filter(Acal_Pio_expressedTEs, grepl("line", TEclass, ignore.case=TRUE))
nrow(aca_Pio_lineexp) # 70 LINE subfamilies expressed

aca_Pio_ltrexp <- dplyr::filter(Acal_Pio_expressedTEs, grepl("LTR", TEclass, ignore.case=TRUE))
nrow(aca_Pio_ltrexp) # 309 LTR subfamilies expressed

aca_Pio_helitronexp <- dplyr::filter(Acal_Pio_expressedTEs, grepl("helitron", TEfam, ignore.case=TRUE))
nrow(aca_Pio_helitronexp) # 1 helitron subfamilies expressed

aca_Pio_sineexp <- dplyr::filter(Acal_Pio_expressedTEs, grepl("SINE", TEclass, ignore.case=TRUE))
nrow(aca_Pio_sineexp) # 3 SINE subfamilies expressed


# for O. niloticus


Onil_expressedTEs <- read.table("Onil_filtered_notNormTEcounts.txt", sep = "\t")
#View(Onil_expressedTEs)
Onil_expressedTEs$TE <- rownames(Onil_expressedTEs)
rownames(Onil_expressedTEs) <- NULL
#View(Onil_expressedTEs)

Onil_expressedTEs$TEfam <- sub("^[^:]*:([^:]*).*", "\\1", Onil_expressedTEs$TE)
Onil_expressedTEs$TEclass <-  sub(".*:", "", Onil_expressedTEs$TE)

#View(Onil_expressedTEs)

Onil_expressedTEs <- filter(Onil_expressedTEs, !TEclass == "Retroposon" & !TEclass == "SINE?")
#View(Onil_expressedTEs)
nrow(Onil_expressedTEs)


colourCount = length(unique(Onil_expressedTEs$TEfam))
getPalette = colorRampPalette(brewer.pal(9, "Greys"))


ONI_EXPRESSEDTES_bp <- ggplot(Onil_expressedTEs, aes(x = TEclass, fill = TEfam)) +
  geom_bar(color = "black") + 
  scale_fill_manual(values = getPalette(colourCount)) +
  labs(y = NULL, x = NULL) + 
  theme_bw()
ONI_EXPRESSEDTES_bp <- ONI_EXPRESSEDTES_bp + guides(fill=guide_legend(title = "TE family"))  +
  theme(axis.text.x = element_text(angle=-45, vjust = 0.5, hjust = 0.1))
ONI_EXPRESSEDTES_bp

oni_dnaexp <- dplyr::filter(Onil_expressedTEs, grepl("DNA", TEclass, ignore.case=TRUE))
nrow(oni_dnaexp) # 450 DNA subfamilies expressed

oni_lineexp <- dplyr::filter(Onil_expressedTEs, grepl("line", TEclass, ignore.case=TRUE))
nrow(oni_lineexp) # 202 LINE subfamilies expressed

oni_ltrexp <- dplyr::filter(Onil_expressedTEs, grepl("LTR", TEclass, ignore.case=TRUE))
nrow(oni_ltrexp) # 80 LTR subfamilies expressed

oni_helitronexp <- dplyr::filter(Onil_expressedTEs, grepl("helitron", TEfam, ignore.case=TRUE))
nrow(oni_helitronexp) # 6 helitron subfamilies expressed

oni_sineexp <- dplyr::filter(Onil_expressedTEs, grepl("SINE", TEclass, ignore.case=TRUE))
nrow(oni_sineexp) # 8 SINE subfamilies expressed


# For P. nyererei

Pnye_expressedTEs <- read.table("Pnye_filtered_notNormTEcounts.txt", sep = "\t")
#View(Pnye_expressedTEs)
Pnye_expressedTEs$TE <- rownames(Pnye_expressedTEs)
rownames(Pnye_expressedTEs) <- NULL
#View(Pnye_expressedTEs)


Pnye_expressedTEs$TEfam <- sub("^[^:]*:([^:]*).*", "\\1", Pnye_expressedTEs$TE)
Pnye_expressedTEs$TEclass <-  sub(".*:", "", Pnye_expressedTEs$TE)

#View(Pnye_expressedTEs)

Pnye_expressedTEs <- filter(Pnye_expressedTEs, !TEclass == "Retroposon" & !TEclass == "SINE?")
#View(Pnye_expressedTEs)
nrow(Pnye_expressedTEs)


colourCount = length(unique(Pnye_expressedTEs$TEfam))
getPalette = colorRampPalette(brewer.pal(9, "Reds"))


PNY_EXPRESSEDTES_bp <- ggplot(Pnye_expressedTEs, aes(x = TEclass, fill = TEfam)) +
  geom_bar(color = "black") + 
  scale_fill_manual(values = getPalette(colourCount)) +
  labs(y = NULL, x = NULL) + 
  theme_bw()
PNY_EXPRESSEDTES_bp <- PNY_EXPRESSEDTES_bp + guides(fill=guide_legend(title = "TE family"))  +
  theme(axis.text.x = element_text(angle=-45, vjust = 0.5, hjust = 0.1))
PNY_EXPRESSEDTES_bp

pny_dnaexp <- dplyr::filter(Pnye_expressedTEs, grepl("DNA", TEclass, ignore.case=TRUE))
nrow(pny_dnaexp) # 379 DNA subfamilies expressed

pny_lineexp <- dplyr::filter(Pnye_expressedTEs, grepl("line", TEclass, ignore.case=TRUE))
nrow(pny_lineexp) # 202 LINE subfamilies expressed

pny_ltrexp <- dplyr::filter(Pnye_expressedTEs, grepl("LTR", TEclass, ignore.case=TRUE))
nrow(pny_ltrexp) # 28 LTR subfamilies expressed

pny_helitronexp <- dplyr::filter(Pnye_expressedTEs, grepl("helitron", TEfam, ignore.case=TRUE))
nrow(pny_helitronexp) # 1 helitron subfamilies expressed

pny_sineexp <- dplyr::filter(Pnye_expressedTEs, grepl("SINE", TEclass, ignore.case=TRUE))
nrow(pny_sineexp) # 7 SINE subfamilies expressed




#### ---- plot all together the summary stats, not by superfamily



species <- c("Acal_Pio","Acal_Pio","Acal_Pio","Acal_Pio","Acal_Pio",
             "Acal","Acal","Acal","Acal","Acal",
             "Abur","Abur","Abur","Abur","Abur",
             "Pnye","Pnye","Pnye","Pnye","Pnye",
             "Onil","Onil","Onil","Onil","Onil")

TEclasses <- c("DNA","LINE","LTR","RC","SINE",
               "DNA","LINE","LTR","RC","SINE",
               "DNA","LINE","LTR","RC","SINE",
               "DNA","LINE","LTR","RC","SINE",
               "DNA","LINE","LTR","RC","SINE")
  


Expressed_TEfams <- c(132,70,309,1,3,
                      396,178,69,1,12,
                      423,213,39,4,10,
                      379,202,28,1,7,
                      450,202,80,6,8)

## check total families:

Aca_normTEcounts
#View(Aca_normTEcounts)
nrow(Aca_normTEcounts)

AcaPio_normTEcounts
#View(AcaPio_normTEcounts)
nrow(AcaPio_normTEcounts)


Abu_normTEcounts
#View(Abu_normTEcounts)
nrow(Abu_normTEcounts)


Pny_normTEcounts
#View(Pny_normTEcounts)
nrow(Pny_normTEcounts)


length(unique(normTECounts3m_Aca_final$TE))
length(unique(normTECounts3m_AcaPio_final$TE))
length(unique(normTECounts3m_Abu_final$TE))
length(unique(normTECounts3m_Pny_final$TE))
length(unique(normTECounts3m_Oni_final$TE))

normTECounts3m_Aca_toCountFams <- normTECounts3m_Aca_final[normTECounts3m_Aca_final$variable=="Aca_tes_rep1",]

normTECounts3m_Aca_toCountFams_DNA <- normTECounts3m_Aca_toCountFams[normTECounts3m_Aca_toCountFams$TEclass == "DNA",]
nrow(normTECounts3m_Aca_toCountFams_DNA)

normTECounts3m_Aca_toCountFams_LINE <- normTECounts3m_Aca_toCountFams[normTECounts3m_Aca_toCountFams$TEclass == "LINE",]
nrow(normTECounts3m_Aca_toCountFams_LINE)

normTECounts3m_Aca_toCountFams_LTR <- normTECounts3m_Aca_toCountFams[normTECounts3m_Aca_toCountFams$TEclass == "LTR",]
nrow(normTECounts3m_Aca_toCountFams_LTR)

normTECounts3m_Aca_toCountFams_RC <- normTECounts3m_Aca_toCountFams[normTECounts3m_Aca_toCountFams$TEclass == "RC",]
nrow(normTECounts3m_Aca_toCountFams_RC)

normTECounts3m_Aca_toCountFams_SINE <- normTECounts3m_Aca_toCountFams[normTECounts3m_Aca_toCountFams$TEclass == "SINE",]
nrow(normTECounts3m_Aca_toCountFams_SINE)



normTECounts3m_Aca_Pio_toCountFams <- normTECounts3m_AcaPio_final[normTECounts3m_AcaPio_final$variable=="Aca_tes_rep1",]

normTECounts3m_Aca_Pio_toCountFams_DNA <- normTECounts3m_Aca_Pio_toCountFams[normTECounts3m_Aca_Pio_toCountFams$TEclass == "DNA",]
nrow(normTECounts3m_Aca_Pio_toCountFams_DNA)

normTECounts3m_Aca_Pio_toCountFams_LINE <- normTECounts3m_Aca_Pio_toCountFams[normTECounts3m_Aca_Pio_toCountFams$TEclass == "LINE",]
nrow(normTECounts3m_Aca_Pio_toCountFams_LINE)

normTECounts3m_Aca_Pio_toCountFams_LTR <- normTECounts3m_Aca_Pio_toCountFams[normTECounts3m_Aca_Pio_toCountFams$TEclass == "LTR",]
nrow(normTECounts3m_Aca_Pio_toCountFams_LTR)

normTECounts3m_Aca_Pio_toCountFams_RC <- normTECounts3m_Aca_Pio_toCountFams[normTECounts3m_Aca_Pio_toCountFams$TEclass == "RC",]
nrow(normTECounts3m_Aca_Pio_toCountFams_RC)

normTECounts3m_Aca_Pio_toCountFams_SINE <- normTECounts3m_Aca_Pio_toCountFams[normTECounts3m_Aca_Pio_toCountFams$TEclass == "SINE",]
nrow(normTECounts3m_Aca_Pio_toCountFams_SINE)



normTECounts3m_Abu_toCountFams <- normTECounts3m_Abu_final[normTECounts3m_Abu_final$variable=="Abu_tes_rep1",]

normTECounts3m_Abu_toCountFams_DNA <- normTECounts3m_Abu_toCountFams[normTECounts3m_Abu_toCountFams$TEclass == "DNA",]
nrow(normTECounts3m_Abu_toCountFams_DNA)

normTECounts3m_Abu_toCountFams_LINE <- normTECounts3m_Abu_toCountFams[normTECounts3m_Abu_toCountFams$TEclass == "LINE",]
nrow(normTECounts3m_Abu_toCountFams_LINE)

normTECounts3m_Abu_toCountFams_LTR <- normTECounts3m_Abu_toCountFams[normTECounts3m_Abu_toCountFams$TEclass == "LTR",]
nrow(normTECounts3m_Abu_toCountFams_LTR)

normTECounts3m_Abu_toCountFams_RC <- normTECounts3m_Abu_toCountFams[normTECounts3m_Abu_toCountFams$TEclass == "RC",]
nrow(normTECounts3m_Abu_toCountFams_RC)

normTECounts3m_Abu_toCountFams_SINE <- normTECounts3m_Abu_toCountFams[normTECounts3m_Abu_toCountFams$TEclass == "SINE",]
nrow(normTECounts3m_Abu_toCountFams_SINE)


normTECounts3m_Pny_toCountFams <- normTECounts3m_Pny_final[normTECounts3m_Pny_final$variable=="Pny_tes_rep1",]

normTECounts3m_Pny_toCountFams_DNA <- normTECounts3m_Pny_toCountFams[normTECounts3m_Pny_toCountFams$TEclass == "DNA",]
nrow(normTECounts3m_Pny_toCountFams_DNA)

normTECounts3m_Pny_toCountFams_LINE <- normTECounts3m_Pny_toCountFams[normTECounts3m_Pny_toCountFams$TEclass == "LINE",]
nrow(normTECounts3m_Pny_toCountFams_LINE)

normTECounts3m_Pny_toCountFams_LTR <- normTECounts3m_Pny_toCountFams[normTECounts3m_Pny_toCountFams$TEclass == "LTR",]
nrow(normTECounts3m_Pny_toCountFams_LTR)

normTECounts3m_Pny_toCountFams_RC <- normTECounts3m_Pny_toCountFams[normTECounts3m_Pny_toCountFams$TEclass == "RC",]
nrow(normTECounts3m_Pny_toCountFams_RC)

normTECounts3m_Pny_toCountFams_SINE <- normTECounts3m_Pny_toCountFams[normTECounts3m_Pny_toCountFams$TEclass == "SINE",]
nrow(normTECounts3m_Pny_toCountFams_SINE)



normTECounts3m_Oni_toCountFams <- normTECounts3m_Oni_final[normTECounts3m_Oni_final$variable=="Oni_tes_rep1",]

normTECounts3m_Oni_toCountFams_DNA <- normTECounts3m_Oni_toCountFams[normTECounts3m_Oni_toCountFams$TEclass == "DNA",]
nrow(normTECounts3m_Oni_toCountFams_DNA)

normTECounts3m_Oni_toCountFams_LINE <- normTECounts3m_Oni_toCountFams[normTECounts3m_Oni_toCountFams$TEclass == "LINE",]
nrow(normTECounts3m_Oni_toCountFams_LINE)

normTECounts3m_Oni_toCountFams_LTR <- normTECounts3m_Oni_toCountFams[normTECounts3m_Oni_toCountFams$TEclass == "LTR",]
nrow(normTECounts3m_Oni_toCountFams_LTR)

normTECounts3m_Oni_toCountFams_RC <- normTECounts3m_Oni_toCountFams[normTECounts3m_Oni_toCountFams$TEclass == "RC",]
nrow(normTECounts3m_Oni_toCountFams_RC)

normTECounts3m_Oni_toCountFams_SINE <- normTECounts3m_Oni_toCountFams[normTECounts3m_Oni_toCountFams$TEclass == "SINE",]
nrow(normTECounts3m_Oni_toCountFams_SINE)

Total_TEfams <- c(133,70,350,1,3,
                  460,212,75,1,12,
                  457,228,39,4,10,
                  434,235,29,1,10,
                  503,226,83,7,9)

  
AllExpressed_families <- data.frame(species,TEclasses,Expressed_TEfams,Total_TEfams)

AllExpressed_families$NotExpressed_TEfams <- AllExpressed_families$Total_TEfams - AllExpressed_families$Expressed_TEfams
AllExpressed_families


# first plot total number of expressed families

AllExpressed_families$species <- ordered(AllExpressed_families$species, levels = c("Acal_Pio","Acal","Abur",
                                                    "Pnye","Onil"))

NumberExpFamilies <- ggplot(AllExpressed_families, aes(x = species, y = Expressed_TEfams, fill = TEclasses)) +
  geom_bar(position="dodge", stat="identity", color = "black") +
  labs(x = NULL, y = "Number of expressed TE families") + 
  scale_fill_manual(values = c("#B77DB3","#7EA6A9","#B4D395","#A7D1F0","#E4ADB5")) +
  theme_bw()
NumberExpFamilies


## then plot proportion
FractionExpFamilies <- ggplot(AllExpressed_families, aes(x = species, y = Expressed_TEfams/Total_TEfams, fill = TEclasses)) +
  geom_bar(position="dodge", stat="identity", color = "black") +
  labs(x = NULL, y = "Fraction of expressed TE families") + 
  scale_fill_manual(values = c("#B77DB3","#7EA6A9","#B4D395","#A7D1F0","#E4ADB5")) +
  theme_bw()
FractionExpFamilies

NumberExpFamilies + FractionExpFamilies



### Overall TE family number plot


species2 <- c("Acal_Pio","Acal_Pio",
              "Acal","Acal",
              "Abur","Abur",
              "Pnye","Pnye",
              "Onil","Onil")

ExpressedStatus <- c("Yes","No",
                     "Yes","No",
                     "Yes","No",
                     "Yes","No",
                     "Yes","No")

# check numbers
AllExpressed_families  %>% group_by(species) %>% summarise(expressed = sum(Expressed_TEfams), notexpressed = sum(NotExpressed_TEfams))


values <- c(515,42,
            656,104,
            689,49,
            617,92,
            746,82)

ExpressedFamilies_General <- data.frame(species2,ExpressedStatus,values)


ExpressedFamilies_General$species2 <- ordered(ExpressedFamilies_General$species2, levels = c("Acal_Pio","Acal","Abur",
                                                                                             "Pnye","Onil"))

AllFams <- ggplot(ExpressedFamilies_General, aes(x = species2, y = values, fill = ExpressedStatus)) +
  geom_bar(position="fill", stat="identity", color = "black") +
  labs(x = NULL, y = "Fraction of expressed TE families") + 
  scale_fill_manual(values = c("#514F4F","#02A704")) +
  theme_bw()
AllFams


AllFams_Stack <- ggplot(ExpressedFamilies_General, aes(x = species2, y = values, fill = ExpressedStatus)) +
  geom_bar(position="stack", stat="identity", color = "black") +
  labs(x = NULL, y = "Fraction of expressed TE families") + 
  scale_fill_manual(values = c("#514F4F","#02A704")) +
  theme_bw()
AllFams_Stack

(AllFams_Stack + AllFams + NumberExpFamilies + FractionExpFamilies) + plot_layout(ncol = 4)


## 15.66 x 2.32
