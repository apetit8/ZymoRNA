source("scripts/functions/RNAseq_analyses.R")
library(ggplot2)
################################################################################
plast.evol <- function(dds, contrast_groups, dfanc, pj = 0.05, sample="Sample", FCmin=2){
  res <- results(dds, contrast = contrast_groups, listValues=c(1/length(contrast_groups[[1]]), -1/length(contrast_groups[[2]])))
  dfres <- as.data.frame(res)
  dfres$plastevol <- rep("X", nrow(dfres))
  dfres$regime <- rep(sample, nrow(dfres))
  dfres$gene_id <- rownames(dfres)
  for (gene in 1:nrow(dfres)) {
    dfres[gene, 7] <- ifelse( (abs(dfres[gene, 2]) >= FCmin && (dfres[gene, 6] <= pj & !is.na(dfres[gene, 6])) && (dfanc[gene, 6] > pj | is.na(dfanc[gene, 6]))), "1New_PL",
                              ifelse( ( (dfres[gene, 6] > pj | is.na(dfres[gene, 6])) && (dfanc[gene, 6] > pj | is.na(dfanc[gene, 6])) ), "nonplastic",  
                                      ifelse( ((dfres[gene, 6] > pj | is.na(dfres[gene, 6])) && dfanc[gene, 6] <= pj  && abs(dfanc[gene, 2]) >= FCmin), "4NP_new",
                                              ifelse( ((dfres[gene, 6] > pj | is.na(dfres[gene, 6])) && dfanc[gene, 6] <= pj  && abs(dfanc[gene, 2]) <= FCmin), "nonplastic",
                                                      ifelse( (abs(dfres[gene, 2]) >= FCmin && (dfres[gene, 6] <= pj & !is.na(dfres[gene, 6])) && round(abs(dfres[gene, 2]),2) > round(abs(dfanc[gene, 2]),2) && abs(dfanc[gene, 2]) >= FCmin), "3Gain",
                                                              ifelse( (abs(dfres[gene, 2]) >= FCmin && (dfres[gene, 6] <= pj & !is.na(dfres[gene, 6])) && round(abs(dfres[gene, 2]),2) == round(abs(dfanc[gene, 2]),2) && abs(dfanc[gene, 2]) >= FCmin), "3Same",
                                                                      ifelse( (abs(dfres[gene, 2]) >= FCmin && (dfres[gene, 6] <= pj & !is.na(dfres[gene, 6])) && round(abs(dfres[gene, 2]),2) < round(abs(dfanc[gene, 2]),2) && abs(dfanc[gene, 2]) >= FCmin), "2Loss",
                                                                              ifelse( (abs(dfres[gene, 2]) <= FCmin && (dfres[gene, 6] <= pj & !is.na(dfres[gene, 6])) && round(abs(dfres[gene, 2]),2) < round(abs(dfanc[gene, 2]),2) && abs(dfanc[gene, 2]) >= FCmin && dfanc[gene, 6] <= pj), "2Loss_FC","nonplastic"
                                                                              ))))))))
  }
  return(dfres)
}
#Script to look at genes differentially expressed for ancestors and lineages
#Threshold
FCmin <- 2
pj <- 0.05
################################################################################
#ONLY 01
################################################################################
samp.info <- read.csv("All_Data/samples_design.csv")
samp.info <- samp.info[order(as.character(samp.info$sample_id), method = c("radix")),]
samp.info <- subset(samp.info, ancester==01)
samp.info[,2:ncol(samp.info)] <- lapply(samp.info[,2:ncol(samp.info)] , factor)
#Counts
raw_counts <- read.csv("scripts/data/counts_01_raw_corrected.csv")
rownames(raw_counts) <- raw_counts$X
raw_counts <- raw_counts[,2:ncol(raw_counts)]
dds01 <- DESeq(DESeqDataSetFromMatrix(raw_counts, colData = samp.info, design = ~ 0 + Pop))
################################################################################
resultsNames(dds01)

contrast_groups <- list(c("PopA01_a_23","Pop01A_b_23_n"), c("PopA01_a_17","Pop01A_b_17_n"))
resanc01 <- results(dds_all, contrast = contrast_groups, listValues=c(1/length(contrast_groups[[1]]), -1/length(contrast_groups[[2]])))
dfanc01 <- as.data.frame(resanc01)
dfanc_01 <- subset(dfanc01, padj <= pj & abs(log2FoldChange) >= FCmin)
dfanc_01$plastevol <- rep("Ancplast", nrow(dfanc_01))
dfanc_01$regime <- rep("Anc01", nrow(dfanc_01))
dfanc_01$gene_id <- rownames(dfanc_01)
##
dfresFa_01 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop01F_a_23_n"), c("Pop01F_a_17_n")), dfanc=dfanc01, pj = pj, sample="F_a", FCmin=FCmin)
dfresFb_01 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop01F_b_23"), c("Pop01F_b_17")), dfanc=dfanc01, pj = pj, sample="F_b", FCmin=FCmin)
dfresFc_01 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop01F_c_23"), c("Pop01F_c_17")), dfanc=dfanc01, pj = pj, sample="F_c", FCmin=FCmin)
##
dfresS17a_01 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop01S17_a_23_n"), c("Pop01S17_a_17_n")), dfanc=dfanc01, pj = pj, sample="S17_a", FCmin=FCmin)
dfresS17b_01 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop01S17_b_23"), c("Pop01S17_b_17")), dfanc=dfanc01, pj = pj, sample="S17_b", FCmin=FCmin)
dfresS17c_01 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop01S17_c_23_n"), c("Pop01S17_c_17_n")), dfanc=dfanc01, pj = pj, sample="S17_c", FCmin=FCmin)
##
dfresS23a_01 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop01S23_a_23"), c("Pop01S23_a_17")), dfanc=dfanc01, pj = pj, sample="S23_a", FCmin=FCmin)
dfresS23b_01 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop01S23_b_23_n"), c("Pop01S23_b_17_n")), dfanc=dfanc01, pj = pj, sample="S23_b", FCmin=FCmin)
dfresS23c_01 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop01S23_c_23_n"), c("Pop01S23_c_17_n")), dfanc=dfanc01, pj = pj, sample="S23_c", FCmin=FCmin)
####
gg1 <- ggplot(subset(rbind(dfanc_01, dfresFa_01,dfresFb_01,dfresFc_01, dfresS17a_01,dfresS17b_01,dfresS17c_01, dfresS23a_01,dfresS23b_01,dfresS23c_01), plastevol!="nonplastic"),
              aes(x=regime, y=after_stat(count), fill=plastevol)) + geom_bar()+ 
  labs(x = "Regime", y = "Gene number", color = "vs ancst", title = "01 : plasticity evolution")+
  theme_bw()+ scale_fill_manual(values=c(  "#56B4E9","#E69F00","yellow2","yellowgreen","darkblue","grey","maroon"))
gg1



################################################################################
#ONLY 44
################################################################################
samp.info <- read.csv("All_Data/samples_design.csv")
samp.info <- samp.info[order(as.character(samp.info$sample_id), method = c("radix")),]
samp.info <- subset(samp.info, ancester==44)
samp.info[,2:ncol(samp.info)] <- lapply(samp.info[,2:ncol(samp.info)], factor)
#Counts
raw_counts <- read.csv("scripts/data/counts_44_raw_corrected.csv")
rownames(raw_counts) <- raw_counts$X
raw_counts <- raw_counts[,2:ncol(raw_counts)]
dds44 <- DESeq(DESeqDataSetFromMatrix(raw_counts, colData = samp.info, design = ~ 0 + Pop))
################################################################################
dds44 <- count_read_all(quantData=quantData, samp.info=samp.info, min=min, ribosomes.array=ribosomes, filename = name, batcheffect=FALSE)
resultsNames(dds44)

contrast_groups <- list(c("Pop44A_a_23", "Pop44A_b_23_n"), c("Pop44A_a_17"))
resanc44 <- results(dds_all, contrast = contrast_groups, listValues=c(1/length(contrast_groups[[1]]), -1/length(contrast_groups[[2]])))
dfanc44 <- as.data.frame(resanc44)
dfanc_44 <- subset(dfanc44, padj <= pj & abs(log2FoldChange) >= FCmin)
dfanc_44$plastevol <- rep("Ancplast", nrow(dfanc_44))
dfanc_44$regime <- rep("Anc44", nrow(dfanc_44))
dfanc_44$gene_id <- rownames(dfanc_44)
##
dfresFa_44 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop44F_a_23"), c("Pop44F_a_17")), dfanc=dfanc44, pj = pj, sample="F_a", FCmin=FCmin)
dfresFb_44 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop44F_b_23_n"), c("Pop44F_b_17_n")), dfanc=dfanc44, pj = pj, sample="F_b", FCmin=FCmin)
dfresFc_44 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop44F_c_23"), c("Pop44F_c_17")), dfanc=dfanc44, pj = pj, sample="F_c", FCmin=FCmin)
##
dfresS17a_44 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop44S17_a_23"), c("Pop44S17_a_17")), dfanc=dfanc44, pj = pj, sample="S17_a", FCmin=FCmin)
dfresS17b_44 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop44S17_b_23_n"), c("Pop44S17_b_17_n")), dfanc=dfanc44, pj = pj, sample="S17_b", FCmin=FCmin)
##
dfresS23a_44 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop44S23_a_23_n"), c("Pop44S23_a_17_n")), dfanc=dfanc44, pj = pj, sample="S23_a", FCmin=FCmin)
dfresS23b_44 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop44S23_b_23_n"), c("Pop44S23_b_17_n")), dfanc=dfanc44, pj = pj, sample="S23_b", FCmin=FCmin)
dfresS23c_44 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop44S23_c_23"), c("Pop44S23_c_17")), dfanc=dfanc44, pj = pj, sample="S23_c", FCmin=FCmin)
####
gg44 <- ggplot(subset(rbind(dfanc_44, dfresFa_44,dfresFb_44,dfresFc_44, dfresS17a_44,dfresS17b_44, dfresS23a_44,dfresS23b_44,dfresS23c_44), plastevol!="nonplastic"),
               aes(x=regime, y=after_stat(count), fill=plastevol)) + geom_bar()+ 
  labs(x = "Regime", y = "Gene number", color = "vs ancst", title = "44 : plasticity evolution")+
  theme_bw()+ scale_fill_manual(values=c(  "#56B4E9","#E69F00","yellow2","yellowgreen","darkblue","grey","maroon"))
gg44
#
dfresF_44 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop44F_a_23", "Pop44F_b_23_n", "Pop44F_c_23"),c("Pop44F_a_17", "Pop44F_b_17_n","Pop44F_c_17")), dfanc=dfanc44, pj = pj, sample="F_a", FCmin=FCmin)
dfresS17_44 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop44S17_a_23", "Pop44S17_b_23_n"), c("Pop44S17_a_17","Pop44S17_b_17_n")), dfanc=dfanc44, pj = pj, sample="S17_a", FCmin=FCmin)
dfresS23_44 <- plast.evol(dds=dds_all, contrast_groups=list(c("Pop44S23_a_23_n", "Pop44S23_b_23_n","Pop44S23_c_23"), c("Pop44S23_a_17_n","Pop44S23_b_17_n","Pop44S23_c_17")), dfanc=dfanc44, pj = pj, sample="S23_a", FCmin=FCmin)
ggplot(subset(rbind(dfanc_44, dfresF_44,dfresS17_44, dfresS23_44), plastevol!="nonplastic"),
       aes(x=regime, y=after_stat(count), fill=plastevol)) + geom_bar()+ 
  labs(x = "Regime", y = "Gene number", color = "vs ancst", title = "44 : plasticity evolution")+
  theme_bw()+ scale_fill_manual(values=c(  "#56B4E9","#E69F00","yellow2","yellowgreen","darkblue","grey","maroon"))

pdf(paste0("figures/new_data_plast_",name,"_evol_FC",FCmin,".pdf"))
gg1
gg44
dev.off()
