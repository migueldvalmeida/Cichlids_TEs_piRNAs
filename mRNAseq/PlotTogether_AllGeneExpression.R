library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(ggpubr)



################################################################
########  Load  dataframes ----  ############################### 
################################################################

Aca_allgenes_melt <- read.table(file = "NormCounts_Calliptera_AllGenes_melted.txt",
                                sep = "\t", header = TRUE)
Aca_allgenes_melt




Abu_allgenes_melt <- read.table(file = "NormCounts_Aburtoni_AllGenes_melted.txt",
                                sep = "\t", header = TRUE)
Abu_allgenes_melt




Pny_allgenes_melt <- read.table(file = "NormCounts_Pnyererei_AllGenes_melted.txt",
                                sep = "\t", header = TRUE)
Pny_allgenes_melt



Oni_allgenes_melt <- read.table(file = "NormCounts_Oniloticus_AllGenes_melted.txt",
                                sep = "\t", header = TRUE)
Oni_allgenes_melt




Allspecies_allgenes <- rbind(Aca_allgenes_melt,Abu_allgenes_melt,
                             Pny_allgenes_melt,Oni_allgenes_melt)

#View(Allspecies_allgenes)


Allspecies_allgenes$species <- ordered(Allspecies_allgenes$species,
                                       levels = c("Aca","Abur",
                                                  "Pnye","Onil"))



Allspecies_allgenes$tissue_terms_changed <- paste(Allspecies_allgenes$species, Allspecies_allgenes$tissue_terms, sep = "_")
#View(Allspecies_allgenes)

Allspecies_allgenes$tissue_terms_changed <- ordered(Allspecies_allgenes$tissue_terms_changed,
                                            levels = c("Aca_tes","Aca_ov",
                                                       "Abur_tes","Abur_ov",
                                                       "Pnye_tes","Pnye_ov",
                                                       "Onil_tes","Onil_ov"))

################################################################
########  Make plots ----  #####################################
################################################################


testes_vs_ov_plot <- ggplot(Allspecies_allgenes, aes(x = species, y = value)) +
  geom_boxplot(aes(fill = tissue_terms_changed), alpha = 1, outlier.shape = NA,
               width = 0.5) +
  ylim(c(-5,24)) + 
  labs(x = NULL, y = "Normalized counts", title = "All protein-coding genes") + 
  scale_fill_manual(values = c("#F7AD6D", "#FEE6CF",
                               "#9FCAE1", "#DFEBF8",
                               "#F29074", "#FCE0D3",
                               "#BDBDBC", "#F0F0F0"))+
  theme_bw()

testes_vs_ov_plot 




################################################################
########  Stats ----  ##########################################
################################################################


Test_fortest_AC <- compare_means(value ~ tissue_terms, data = Allspecies_allgenes[Allspecies_allgenes$species == "Aca",], p.adjust.method = "BH")
View(Test_fortest_AC)
length(unique(Allspecies_allgenes[Allspecies_allgenes$species == "Aca",]$ensembl_gene_id))

Test_fortest_AB <- compare_means(value ~ tissue_terms, data = Allspecies_allgenes[Allspecies_allgenes$species == "Abur",], p.adjust.method = "BH")
View(Test_fortest_AB)
length(unique(Allspecies_allgenes[Allspecies_allgenes$species == "Abur",]$ensembl_gene_id))

Test_fortest_PN <- compare_means(value ~ tissue_terms, data = Allspecies_allgenes[Allspecies_allgenes$species == "Pnye",], p.adjust.method = "BH")
View(Test_fortest_PN)
length(unique(Allspecies_allgenes[Allspecies_allgenes$species == "Pnye",]$ensembl_gene_id))


Test_fortest_ON <- compare_means(value ~ tissue_terms, data = Allspecies_allgenes[Allspecies_allgenes$species == "Onil",], p.adjust.method = "BH")
View(Test_fortest_ON)
length(unique(Allspecies_allgenes[Allspecies_allgenes$species == "Onil",]$ensembl_gene_id))


