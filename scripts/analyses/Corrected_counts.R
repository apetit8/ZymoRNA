#Script to make corrected counts
source("scripts/functions/RNAseq_analyses.R")
library(corrplot)
library(RColorBrewer)
################################################################################
#ONLY 01
################################################################################
min <-  10 #min = reads mean min to be considered as outliers
#Threshold
FCmin <- 2
pj <- 0.05
#Data
quantData <- list.files(path = "All_Data/COUNT/", full.names=TRUE, recursive = TRUE, pattern = "t_data.ctab")
quantData <- quantData[order(as.character(quantData), method = c("radix"))]
names(quantData) <- str_split(quantData, "/", simplify = TRUE)[,4] #select column corresponding to sample name
quantData <- quantData[c(1:43)]
#Design info
samp.info <- read.csv("All_Data/samples_design.csv")
samp.info <- samp.info[order(as.character(samp.info$sample_id), method = c("radix")),]
samp.info <- subset(samp.info, ancester==01)
samp.info[,2:ncol(samp.info)] <- lapply(samp.info[,2:ncol(samp.info)] , factor)
#"Group" or coloumn of ref for ComBat_seq
group <- 8
#List of ribosomes from annotation file
ribosomes <- str_split(subset.pattern(read.csv(file="All_Data/Annotation_r.csv", sep='\t'), "ribosom", 4)[,3], as.character("\\."), n=2, simplify = TRUE)[,1]
name <- "01_ComBatYearBatch"
################################################################################
dds_all <- count_read_all(quantData=quantData, samp.info=samp.info, min=min, ribosomes.array=ribosomes, filename = name, batcheffect=FALSE)
counts_raw_corrected <- ComBat_for_batch(counts(estimateSizeFactors(dds_all), normalized=FALSE), samp.info=samp.info, BiorepTemp=TRUE, group=group,
                                         T23=c(3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36,39,40,42,43),
                                         T17=c(1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38,41))
dds_all <- DESeq(DESeqDataSetFromMatrix(counts_raw_corrected, colData = samp.info, design = ~ 0 + Pop))
write.csv(counts(estimateSizeFactors(dds_all), normalized=TRUE), "scripts/data/counts_01_nomalized_corrected.csv")
write.csv(counts_raw_corrected, "scripts/data/counts_01_raw_corrected.csv")

################################################################################
#ONLY 44
################################################################################
#Data
quantData <- list.files(path = "All_Data/COUNT/", full.names=TRUE, recursive = TRUE, pattern = "t_data.ctab")
quantData <- quantData[order(as.character(quantData), method = c("radix"))]
names(quantData) <- str_split(quantData, "/", simplify = TRUE)[,4] #select column corresponding to sample name
quantData <- quantData[c(44:79)]
#Design info
samp.info <- read.csv("All_Data/samples_design.csv")
samp.info <- samp.info[order(as.character(samp.info$sample_id), method = c("radix")),]
samp.info <- subset(samp.info, ancester==44)
samp.info[,2:ncol(samp.info)] <- lapply(samp.info[,2:ncol(samp.info)] , factor)
#List of ribosomes from annotation file
name <- "44_ComBatYearBatch"
################################################################################
dds_all <- count_read_all(quantData=quantData, samp.info=samp.info, min=min, ribosomes.array=ribosomes, filename = name)
counts_raw_corrected <- ComBat_for_batch(counts(estimateSizeFactors(dds_all), normalized=FALSE), samp.info=samp.info, BiorepTemp=TRUE, group=group,
                                         T23=c(3,4,7,8,11,12,15,18,19,21,22,25,26,29,30,31,32,35,36),
                                         T17=c(1,2,5,6,9,10,13,14,16,17,20,23,24,27,28,33,34))
dds_all <-DESeq(DESeqDataSetFromMatrix(counts_raw_corrected, colData = samp.info, design = ~ 0 + Pop))
write.csv(counts(estimateSizeFactors(dds_all), normalized=TRUE), "scripts/data/counts_44_nomalized_corrected.csv")
write.csv(counts_raw_corrected, "scripts/data/counts_44_raw_corrected.csv")
