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
library(RColorBrewer)
library(ashr)
library(ggrepel)
library(ggpubr)
library(patchwork)


#######################################
#Load Counts  ---- 
#######################################


# Choose the files 
TEcounts = read.table("./PathToYourFeatureCounts_Results", fill = TRUE, header = TRUE)

View(TEcounts)

# manually add names to columns

colnames(TEcounts)[1] <- c("TEfam") 
  
TEcounts <- TEcounts[,c(1,7:21)]

# change colnames according to need.
#colnames(TEcounts)[2:16] <- str_replace(colnames(TEcounts)[2:16],
#                                       '.+20_mapped_STAR_filtered24_35Reads..(.+)', '\\1')
#colnames(TEcounts)

#colnames(TEcounts)[2:16] <- str_replace(colnames(TEcounts)[2:16],
 #                                      '(.*)_T_T_.*', '\\1')
  
colnames(TEcounts)



###############################################
# make sample info  and TE counts matrix ---- #
###############################################

sampleInfoTE = as.data.frame(colnames(TEcounts[2:16]))
#sampleInfoTE$group = sub("_.*", "", sampleInfoTE$`colnames(TEcounts[2:16])`)
sampleInfoTE$group = sub("_n.*", "", sampleInfoTE$`colnames(TEcounts[2:16])`)
sampleInfoTE$group = sub(".*_", "", sampleInfoTE$group)

colnames(sampleInfoTE) = c("sample", "group")
row.names(sampleInfoTE) = sampleInfoTE$sample

TEcounts_matrix <- as.data.frame(TEcounts)
rownames(TEcounts_matrix) <- TEcounts_matrix$TE
head(TEcounts_matrix)

TEcounts_matrix <- TEcounts_matrix[,2:16]
head(TEcounts_matrix)


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


write.table(normTECounts,
            "NormTEcounts.txt",sep="\t",
            col.names=T, quote=F, row.names=T)


#######################################
# PCA  ---- 
#######################################

plotPCAEx(rlogObjTE,PCx = 1,PCy = 2,cond = "group",ntop = 500, T)
plotPCAEx(rlogObjTE,PCx = 2,PCy = 3,cond = "group",ntop = 500, T)


#######################################
# TE expression per tissue  ---- 
#######################################


normTECounts2 <- as.data.frame(normTECounts)
#View(normTECounts2)
normTECounts2$TE <- rownames(normTECounts2)
rownames(normTECounts2) <- NULL

normTECounts2 <- left_join(x = normTECounts2, y = TEFullNametoLeftJoin, by = c("TE" = "TEsubfam"))
normTECounts2 <- normTECounts2[,c(2:15,17)]
#View(normTECounts2)
#View(TEFullNametoLeftJoin)

nrow(normTECounts2)
normTECountsm <- melt(normTECounts2)
#View(normTECountsm)
nrow(normTECountsm) 


normTECountsm$TEclass <-  sub(".*:", "", normTECountsm$TE_list)
normTECountsm$TEfam <- sub("^[^:]*:([^:]*).*", "\\1", normTECountsm$TE_list)
#View(normTECountsm)

normTECountsm$tissue <-  sub("_n.*", "", normTECountsm$variable)
normTECountsm$tissue <-  sub(".*_", "", normTECountsm$tissue)

#View(normTECountsm)

nrow(normTECountsm)
normTECountsm_final <- filter(normTECountsm, !TEclass == "Retroposon" & !TEclass == "SINE?")
nrow(normTECountsm_final)


unique(normTECountsm_final$TEclass)

normTECountsm_final$TEclass <- ordered(normTECountsm_final$TEclass,
                                        levels=c("DNA","LINE","LTR",
                                                 "RC","SINE"))


write.table(normTECountsm_final, file = "Malawi_sRNA_TEnormcounts_ALL_melted_PIO.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

normTECountsm_final_AcalGonads <- filter(normTECountsm_final, !tissue == "mus")
#View(normTECountsm_final_AcalGonads)
normTECountsm_final_AcalGonads <- normTECountsm_final_AcalGonads[!grepl("Tma",normTECountsm_final_AcalGonads$variable),]
normTECountsm_final_AcalGonads <- normTECountsm_final_AcalGonads[!grepl("Mze",normTECountsm_final_AcalGonads$variable),]
#View(normTECountsm_final_AcalGonads)

#View(normTECountsm_final)

plotTECLASSES <- ggplot(normTECountsm_final_AcalGonads, aes(x = TEclass, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(0,19)) + 
  labs(title = "Expression per TEclass", y = "Log-Normalized counts", x = "TE class") + 
  scale_fill_brewer(palette = "Oranges") +
   theme_bw()

plotTECLASSES 



# 6.69 x 5.72

# stats
fortestClassesDNA <- normTECountsm_final_AcalGonads[grepl("DNA", normTECountsm_final_AcalGonads$TEclass),]
nrow(fortestClassesDNA)

fortestClassesLINE <- normTECountsm_final_AcalGonads[grepl("LINE", normTECountsm_final_AcalGonads$TEclass),]
nrow(fortestClassesLINE)

fortestClassesLTR <- normTECountsm_final_AcalGonads[grepl("LTR", normTECountsm_final_AcalGonads$TEclass),]
nrow(fortestClassesLTR)

fortestClassesRC <- normTECountsm_final_AcalGonads[grepl("RC", normTECountsm_final_AcalGonads$TEclass),]
nrow(fortestClassesRC)

fortestClassesSINE <- normTECountsm_final_AcalGonads[grepl("SINE", normTECountsm_final_AcalGonads$TEclass),]
nrow(fortestClassesSINE)


testDNA <- compare_means(value ~ tissue, data = fortestClassesDNA, p.adjust.method = "BH")
View(testDNA)

testLINE <- compare_means(value ~ tissue, data = fortestClassesLINE, p.adjust.method = "BH")
View(testLINE)

testLTR <- compare_means(value ~ tissue, data = fortestClassesLTR, p.adjust.method = "BH")
View(testLTR)

testRC <- compare_means(value ~ tissue, data = fortestClassesRC, p.adjust.method = "BH")
View(testRC)

testSINE <- compare_means(value ~ tissue, data = fortestClassesSINE, p.adjust.method = "BH")
View(testSINE)




### PER FAMILY ->DNA

normTEcounts_perfamALL_DNAm <- filter(normTECountsm_final_AcalGonads, TEclass == "DNA")
#View(normTEcounts_perfamALL_DNAm)


plotTEFamDNA <- ggplot(normTEcounts_perfamALL_DNAm, aes(x = tissue, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.5) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(0,19)) + 
  labs(x = NULL, y = "Log-Normalized counts", title = "DNA transposons") + 
  #labs(title = "Brood size at 25°C", y = "Brood Size", x = NULL) + 
  scale_fill_brewer(palette = "Oranges") +
  theme_bw()

plotTEFamDNA + facet_wrap(~TEfam) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# 6.69 x 5.72

write.table(normTEcounts_perfamALL_DNAm, file = "DNATEs_NormCounts.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)



# Retrotransposons

normTEcounts_perfamALL_retro1 <- filter(normTECountsm_final_AcalGonads, TEclass == "LINE")
#View(normTEcounts_perfamALL_retro1)

normTEcounts_perfamALL_retro2 <- filter(normTECountsm_final_AcalGonads, TEclass == "LTR")
#View(normTEcounts_perfamALL_retro2)

normTEcounts_perfamALL_retro3 <- filter(normTECountsm_final_AcalGonads, TEclass == "Retroposon")
#View(normTEcounts_perfamALL_retro3)

normTEcounts_perfamALL_retro4 <- filter(normTECountsm_final_AcalGonads, TEclass == "SINE")
#View(normTEcounts_perfamALL_retro4)

normTEcounts_perfamALL_retroTEs <- rbind(normTEcounts_perfamALL_retro1,normTEcounts_perfamALL_retro2,
                                         normTEcounts_perfamALL_retro3,normTEcounts_perfamALL_retro4)

#View(normTEcounts_perfamALL_retroTEs)


plotTEFamRetro <- ggplot(normTEcounts_perfamALL_retroTEs, aes(x = tissue, y = value)) +
  geom_boxplot(aes(fill = tissue), alpha = 1, outlier.shape = NA,
               width = 0.5) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(0,16)) + 
  labs(x = NULL, y = "Log-Normalized counts", title = "Retrotransposons") + 
  #labs(title = "Brood size at 25°C", y = "Brood Size", x = NULL) + 
  scale_fill_brewer(palette = "Oranges") +
  theme_bw()

plotTEFamRetro + facet_wrap(~TEfam) + theme(axis.text.x = element_text(angle = 315, vjust = .4, hjust=0))


# 6.69 x 5.72


write.table(normTEcounts_perfamALL_retroTEs, file = "retroTEs_NormCounts.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)



#####################################################################
# TE expression active vs inactive ----  ############################
#####################################################################

normTECountsm_final_AcalGondas_forFilt <- normTECountsm_final_AcalGonads
normTECountsm_final_AcalGondas_forFilt$TEsubfam <- sub(":.*", "", normTECountsm_final_AcalGondas_forFilt$TE_list)



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


plotTE_Gypsy_Active <- ggplot(normTECountsm_final_AcalGonads_gypsy_final, aes(x = tissue, y = value, fill = MobStatus)) +
  geom_boxplot(aes(fill = MobStatus), alpha = 1, outlier.shape = NA,
               width = 0.75) +
  #geom_point(position = position_jitter(w = 0.2, h = 0), aes(fill= variable),
  #shape = 21) +
  ylim(c(0,15)) + 
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
  ylim(c(0,15)) + 
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
  ylim(c(0,15)) + 
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
### NO REX-BABAR-1 in reference... do not move on with this one.





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
  ylim(c(0,17)) + 
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
  ylim(c(0,15)) + 
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
  ylim(c(0,15)) + 
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
  ylim(c(0,17)) + 
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
  ylim(c(0,15)) + 
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
  ylim(c(0,15)) + 
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
plotTE_Gypsy_Active2 <- plotTE_Gypsy_Active + theme(legend.position = "none")
plotTE_ERV1_Active2 <- plotTE_ERV1_Active + theme(legend.position = "none")
plotTE_hATAc_Active2 <- plotTE_hATAc_Active + theme(legend.position = "none")

### alltogether:

wrap_plots(plotTE_DIRS1_Active2, plotTE_hATTiP100_Active2, plotTE_PiggyBac5_Active2,plotTE_Gypsy_Active2,plotTE_hATAc_Active2,
           plotTE_CMC_Active2,plotTE_ERV1_Active2,plotTE_Maverick_Active2,plotTE_TcMarTc2_Active2,plot_spacer(), nrow = 2)

## 10.70 x 4.66


