source("scripts/functions/RNAseq_analyses.R")
library(ggplot2)
################################################################################
plast.evol <- function(dds, contrast_groups, dfanc, pj = 0.05, sample="Sample", FCmin=2){
  res <- results(dds, contrast = contrast_groups)
  dfres <- as.data.frame(res)
  dfres$plastevol <- rep("X", nrow(dfres))
  dfres$regime <- rep(sample, nrow(dfres))
  dfres$gene_id <- rownames(dfres)
  for (gene in 1:nrow(dfres)) {
    dfres[gene, 7] <- ifelse( (abs(dfres[gene, 2]) >= FCmin && (dfres[gene, 6] <= pj & !is.na(dfres[gene, 6])) && (dfanc[gene, 6] > pj | is.na(dfanc[gene, 6]))), "1New_PL",
                              ifelse( ( (dfres[gene, 6] > pj | is.na(dfres[gene, 6])) && (dfanc[gene, 6] > pj | is.na(dfanc[gene, 6])) ), "nonplastic",  
                                      ifelse( ((dfres[gene, 6] > pj | is.na(dfres[gene, 6])) && dfanc[gene, 6] <= pj  && !is.na(dfanc[gene, 6]) && abs(dfanc[gene, 2]) >= FCmin), "4NP_new",
                                              ifelse( ((dfres[gene, 6] > pj | is.na(dfres[gene, 6])) && abs(dfanc[gene, 2]) <= FCmin), "nonplastic",
                                                      ifelse( (abs(dfres[gene, 2]) >= FCmin && (dfres[gene, 6] <= pj && !is.na(dfres[gene, 6])) && round(abs(dfres[gene, 2]),2) > round(abs(dfanc[gene, 2]),2) && abs(dfanc[gene, 2]) >= FCmin && !is.na(dfanc[gene, 6])), "3Gain",
                                                              ifelse( (abs(dfres[gene, 2]) >= FCmin && (dfres[gene, 6] <= pj && !is.na(dfres[gene, 6])) && round(abs(dfres[gene, 2]),2) == round(abs(dfanc[gene, 2]),2) && abs(dfanc[gene, 2]) >= FCmin && !is.na(dfanc[gene, 6])), "3Same",
                                                                      ifelse( (abs(dfres[gene, 2]) >= FCmin && (dfres[gene, 6] <= pj && !is.na(dfres[gene, 6])) && round(abs(dfres[gene, 2]),2) < round(abs(dfanc[gene, 2]),2) && abs(dfanc[gene, 2]) >= FCmin && !is.na(dfanc[gene, 6])), "2Loss",
                                                                              ifelse( (abs(dfres[gene, 2]) <= FCmin && (dfres[gene, 6] <= pj && !is.na(dfres[gene, 6])) && round(abs(dfres[gene, 2]),2) < round(abs(dfanc[gene, 2]),2) && abs(dfanc[gene, 2]) >= FCmin && !is.na(dfanc[gene, 6]) && dfanc[gene, 6] <= pj), "2Loss_FC","nonplastic"
                                                                              ))))))))
  }
  return(dfres)
}
#Script to look at genes differentially expressed for ancestors and lineages
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
dds_all <- count_read_all(quantData=quantData, samp.info=samp.info, min=min, ribosomes.array=ribosomes, filename = name, batcheffect=TRUE)
resultsNames(dds_all)

# m1 <- model.matrix(~ 1 + Year + Year:Biorepnested + Year:Evoltemp, colData(dds_all))
# deg <- DESeq(dds_all, full=m1, betaPrior = F)
# resultsNames(deg)
# dds_all <- deg

resanc01 <- results(dds_all, contrast = list(c("Year2016.EvoltempA01_23", "Year2023.EvoltempA01_23"), c("Year2016.EvoltempA01_17","Year2023.EvoltempA01_17")))
dfanc01 <- as.data.frame(resanc01)
dfanc_01 <- subset(dfanc01, padj <= pj & abs(log2FoldChange) >= FCmin & !is.null(padj))
dfanc_01$plastevol <- rep("Ancplast", nrow(dfanc_01))
dfanc_01$regime <- rep("Anc01", nrow(dfanc_01))
dfanc_01$gene_id <- rownames(dfanc_01)


dfresF_01 <- plast.evol(dds=dds_all, contrast_groups=list(c("Year2016.Evoltemp01F_23", "Year2023.Evoltemp01F_23")), dfanc=dfanc01, pj = pj, sample="F_a", FCmin=FCmin)
##
dfresS17_01 <- plast.evol(dds=dds_all, contrast_groups=list(c("Year2016.Evoltemp01S17_23", "Year2023.Evoltemp01S17_23"), c("Year2016.Evoltemp01S17_17","Year2023.Evoltemp01S17_17")), dfanc=dfanc01, pj = pj, sample="S17_a", FCmin=FCmin)
##
dfresS23_01 <- plast.evol(dds=dds_all, contrast_groups=list(c("Year2016.Evoltemp01S23_23", "Year2023.Evoltemp01S23_23"), c("Year2016.Evoltemp01S23_17","Year2023.Evoltemp01S23_17")), dfanc=dfanc01, pj = pj, sample="S23_a", FCmin=FCmin)
####
gg1 <- ggplot(subset(rbind(dfanc_01, dfresF_01,dfresS17_01, dfresS23_01), plastevol!="nonplastic"),
              aes(x=regime, y=after_stat(count), fill=plastevol)) + geom_bar()+ 
  labs(x = "Regime", y = "Gene number", color = "vs ancst", title = "01 : plasticity evolution")+
  theme_bw()+ scale_fill_manual(values=c(  "#56B4E9","#E69F00","yellow2","yellowgreen","darkblue","grey","maroon"))
gg1

pdf(paste0("figures/new_data_plast_evolBatch01_FC",FCmin,".pdf"))
gg1
dev.off()


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
name <- "44"
################################################################################
dds_all <- count_read_all(quantData=quantData, samp.info=samp.info, min=min, ribosomes.array=ribosomes, filename = name, batcheffect=FALSE)
resultsNames(dds_all)

m1 <- model.matrix(~ 0 + Year + Year:Biorepnested + Evoltemp, colData(dds_all))
deg <- DESeq(dds_all, full=m1, betaPrior = F)
resultsNames(deg)
dds_all <- deg
##
resanc44 <- results(dds_all, contrast = list(c("Evoltemp44A_23"))) # "Evoltemp44A_17" is the reference (not in the column names of m1)
dfanc44 <- as.data.frame(resanc44)
sum(is.na(dfanc44[,5]))
dfanc_44 <- subset(dfanc44, padj <= pj & abs(log2FoldChange) >= FCmin  & !is.na(padj))
dfanc_44$plastevol <- rep("Ancplast", nrow(dfanc_44))
dfanc_44$regime <- rep("Anc44", nrow(dfanc_44))
dfanc_44$gene_id <- rownames(dfanc_44)
##

dfresF_44 <- plast.evol(dds=dds_all, contrast_groups=list(c("Evoltemp44F_23"), c("Evoltemp44F_17")), dfanc=dfanc44, pj = pj, sample="F_a", FCmin=FCmin)
##
dfresS17_44 <- plast.evol(dds=dds_all, contrast_groups=list(c("Evoltemp44S17_23"), c("Evoltemp44S17_17")), dfanc=dfanc44, pj = pj, sample="S17_a", FCmin=FCmin)
##
dfresS23_44 <- plast.evol(dds=dds_all, contrast_groups=list(c("Evoltemp44S23_23"), c("Evoltemp44S23_17")), dfanc=dfanc44, pj = pj, sample="S23_a", FCmin=FCmin)
####
gg44 <- ggplot(subset(rbind(dfanc_44, dfresF_44,dfresS17_44, dfresS23_44), plastevol!="nonplastic"),
               aes(x=regime, y=after_stat(count), fill=plastevol)) + geom_bar()+ 
  labs(x = "Regime", y = "Gene number", color = "vs ancst", title = "44 : plasticity evolution")+
  theme_bw()+ scale_fill_manual(values=c(  "#56B4E9","#E69F00","yellow2","yellowgreen","darkblue","grey","maroon"))
gg44

pdf(paste0("figures/new_data_plast_evolBatch_FC",FCmin,".pdf"))
gg1
gg44
dev.off()




