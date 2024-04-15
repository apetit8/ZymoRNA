library(readr)
library(DESeq2)
library(tximport)
library(tximportData)
library(stringr)
library(ggplot2)
################################################################################

count_read_all <- function(quantData, samp.info, min=10, ribosomes.array=ribosomes, filename="all", readLength=75, batcheffect=FALSE){
  #Function to run a DESeq2 analyses on raw read counts.
  #readLength =-s argument of stringtie, default is 75
  tmp <- read_tsv(quantData[1])
  tx2gene <- tmp[, c("t_name", "gene_name")]
  txi1 <- tximport(files=quantData, type = "stringtie", tx2gene = tx2gene, readLength=readLength) 
  tx.count1 <- as.data.frame(txi1[["counts"]])
  
  
  if(identical(samp.info$sample_id, names(quantData))==TRUE){
    samp.info[,2:ncol(samp.info)] <- lapply(samp.info[,2:ncol(samp.info)] , factor)
    #Tximport data with : DeSEQ2
    if(batcheffect) reads.raw.count <- DESeqDataSetFromTximport(txi = txi1, colData=samp.info,
                                                                   design =  ~ 1 + Year + Year:Batch2 + Year:Evol) #batch effect accounted for ~ Batch + Evol
    else reads.raw.count <- DESeqDataSetFromTximport(txi = txi1, colData=samp.info,
                                                design = ~ 0 + Pop) #Contrast null
    #Remove genes with a read mean inferior at x
    reads.raw.count <- reads.raw.count[rowMeans(counts(reads.raw.count))>=min,]
    # #Remove outliers
    #reads.raw.count <- reads.raw.count[rowSums(counts(reads.raw.count))<=max,]
    #Remove ribosomes
    reads.raw.count <- reads.raw.count[!names(reads.raw.count)%in%ribosomes,] 
    
    write.csv(counts(reads.raw.count), paste0("scripts/data/",filename,"_raw_read_counts.csv"))
    
    # reads.raw.count2 <- estimateSizeFactors(reads.raw.count)
    # sizeFactors(reads.raw.count2)
    reads.normalized.counts <- counts(estimateSizeFactors(reads.raw.count), normalized=TRUE)
    #Filter with edgeR to improve normalization
    m_design <- model.matrix(~ 1 + Year, samp.info)
    print(length(which(edgeR::filterByExpr(reads.normalized.counts, design=m_design)==FALSE)))
    reads.normalized.counts <- reads.normalized.counts[edgeR::filterByExpr(reads.normalized.counts, design=m_design),]
    write.csv(reads.normalized.counts, file=paste0("scripts/data/",filename,"_normalized_counts.csv"))

    dds <- DESeq(reads.raw.count)

    return(dds)
  }else(print("Sample_info and sample files are not in the same order"))
}

#Function to subset by string pattern
subset.pattern <- function(data, patterns, var, keep.found = TRUE, useBytes = TRUE){
  data$y <- grepl(pattern = paste0(patterns, collapse="|"), x = data[, var],
                  useBytes = useBytes)
  subdata <- subset(data, y == keep.found)
  subdata$y <- NULL
  subdata
}


plotPCA.custom <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE, PCx=1, PCy=2) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, PCx], PC2 = pca$x[, PCy], group = group, intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[PCx:PCy]
    return(d)
  }
  ggplot(data=d, aes(x=PC1, y=PC2, color=group)) +
    geom_point(size=3) + 
    xlab(paste0("PC",PCx,": ", round(percentVar[PCx] * 100), "% variance")) +
    ylab(paste0("PC",PCy,": ", round(percentVar[PCy] * 100), "% variance")) +
    coord_fixed()+theme_bw()
}



ComBat_for_batch <- function(counts_raw, BiorepTemp=FALSE, samp.info, group=8, T23=c(3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36,39,40,42,43), 
                             T17=c(1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38,41)){
  #group=samp.info$Evol initially but can be changed
  samp.info[,2:ncol(samp.info)] <- lapply(samp.info[,2:ncol(samp.info)] , factor)
  gene_exp_adj <- sva::ComBat_seq(counts = as.matrix(counts_raw), batch = samp.info$Year, 
                                  group = samp.info[,group]) #CORRECTS REALLY WELL THE YEAR EFFECT ! YAHOU
  dds_all <- DESeqDataSetFromMatrix(gene_exp_adj, colData = samp.info, design = ~ 0 + Pop)
  counts_raw <- counts(estimateSizeFactors(dds_all), normalized=FALSE)
  
  if(BiorepTemp==TRUE){
    #Correct batch effect within 23°C
    counts_raw23 <- counts_raw[,T23]
    colnames(counts_raw23)
    samp.info23 <- subset(samp.info, temp==23)
    samp.info23$Biorep <- factor(samp.info23$Biorep)
    samp.info23[,group] <- factor(samp.info23[,group])
    gene_exp_adj23 <- sva::ComBat_seq(counts = counts_raw23, batch = samp.info23$Biorep, 
                                      group = samp.info23[,group])
    #Correct batch effect within 23°C
    counts_raw17 <- counts_raw[,T17]
    colnames(counts_raw17)
    samp.info17 <- subset(samp.info, temp==17)
    samp.info17$Biorep <- factor(samp.info17$Biorep)
    samp.info17[,group] <- factor(samp.info17[,group])
    gene_exp_adj17 <- sva::ComBat_seq(counts = counts_raw17, batch = samp.info17$Biorep, 
                                      group = samp.info17[,group])
    #
    counts_raw_corrected <- cbind(as.data.frame(gene_exp_adj23), as.data.frame(gene_exp_adj17))
    counts_raw <- counts_raw_corrected[order(colnames(counts_raw_corrected), method = c("radix"))]
  }
  
  return(counts_raw)
}



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