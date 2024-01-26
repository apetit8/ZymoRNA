source("scripts/functions/RNAseq_analyses.R")
################################################################################
#ONLY 01
################################################################################
#List of ribosomes from annotation file
ribosomes <- str_split(subset.pattern(read.csv(file="All_Data/Annotation_r.csv", sep='\t'), "ribosom", 4)[,3], as.character("\\."), n=2, simplify = TRUE)[,1]
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

name <- "01"
################################################################################
#Make DESeq2 object with design: ~ 0 + Pop
dds_all <- count_read_all(quantData=quantData, samp.info=samp.info, min=min, ribosomes.array=ribosomes, filename = name, batcheffect=FALSE)
#Warning: "some variables in design formula are characters, converting to factors" has no consequences
resultsNames(dds_all)

#to test new design:
m1 <- model.matrix(~ 1 + Year + Year:Biorepnested + Year:Evoltemp, colData(dds_all))
deg <- DESeq(dds_all, full=m1, betaPrior = F)
resultsNames(deg)
dds_all <- deg

#Manual contrast (have to be column names of the design matrix)
results(dds_all, contrast = list(c("Year2016.EvoltempA01_23", "Year2023.EvoltempA01_23"), c("Year2016.EvoltempA01_17","Year2023.EvoltempA01_17")))



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
name <- "44"
################################################################################
#Make DESeq2 object with design: ~ 0 + Pop
dds_all <- count_read_all(quantData=quantData, samp.info=samp.info, min=min, ribosomes.array=ribosomes, filename = name, batcheffect=FALSE)
#Warning: "some variables in design formula are characters, converting to factors" has no consequences
resultsNames(dds_all)

m1 <- model.matrix(~ 0 + Year + Year:Biorepnested + Evoltemp, colData(dds_all))
deg <- DESeq(dds_all, full=m1, betaPrior = F)
resultsNames(deg)
dds_all <- deg

#Manual contrast (have to be column names of the design matrix)
resanc44 <- results(dds_all, contrast = list(c("Evoltemp44A_23"))) # one group = versus the reference (i think)

