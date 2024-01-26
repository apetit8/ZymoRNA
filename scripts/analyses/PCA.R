source("scripts/functions/RNAseq_analyses.R")
################################################################################
#List of ribosomes from annotation file
ribosomes <- str_split(subset.pattern(read.csv(file="All_Data/Annotation_r.csv", sep='\t'), "ribosom", 4)[,3], as.character("\\."), n=2, simplify = TRUE)[,1]
#outsiders filtering parameters
min <-  10 #min = reads mean min to be considered as outliers
max <-  4*10^6 #max = reads mean max to be considered as outliers
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
#List of ribosomes from annotation file
name <- "01_wtt_outliers"
################################################################################
dds_all <- count_read_all(quantData=quantData, samp.info=samp.info, min=min, max=max, ribosomes.array=ribosomes, filename = name)
#Warning: "some variables in design formula are characters, converting to factors" has no consequences
resultsNames(dds_all)

vsdata <- vst(dds_all, blind=FALSE)

pdf(paste0("figures/PCA/PCA_",name,".pdf"))
plotPCA.custom(vsdata, intgroup="Evol", PCx=1, PCy=2) 
plotPCA.custom(vsdata, intgroup="Evol", PCx=3, PCy=4)  
plotPCA.custom(vsdata, intgroup="Year", PCx=1, PCy=2) 
plotPCA.custom(vsdata, intgroup="Year", PCx=3, PCy=4) 
plotPCA.custom(vsdata, intgroup="Batch", PCx=1, PCy=2) 
plotPCA.custom(vsdata, intgroup="Batch", PCx=3, PCy=4) 
plotPCA.custom(vsdata, intgroup="temp", PCx=1, PCy=2) 
plotPCA.custom(vsdata, intgroup="temp", PCx=3, PCy=4) 
plotPCA.custom(vsdata, intgroup="regime1", PCx=1, PCy=2) 
plotPCA.custom(vsdata, intgroup="regime1", PCx=3, PCy=4) 
plotPCA.custom(vsdata, intgroup="Pop", PCx=1, PCy=2) 
plotPCA.custom(vsdata, intgroup="Pop", PCx=3, PCy=4) 
dev.off()
