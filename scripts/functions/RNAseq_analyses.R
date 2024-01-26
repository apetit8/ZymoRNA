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
  #gicler ribosomes ici
  
  
  if(identical(samp.info$sample_id, names(quantData))==TRUE){
    
    #Tximport data with : DeSEQ2
    if(batcheffect) reads.raw.count <- DESeqDataSetFromTximport(txi = txi1, colData=samp.info,
                                                                   design =  ~ 1 + Year + Year:Biorepnested + Year:Evoltemp) #batch effect accounted for ~ Batch + Evol
    else reads.raw.count <- DESeqDataSetFromTximport(txi = txi1, colData=samp.info,
                                                design = ~ 0 + Pop) #Contrast null
    #Remove genes with a read mean inferior at x
    reads.raw.count <- reads.raw.count[rowMeans(counts(reads.raw.count))>=min,]
    #Remove ribosomes
    reads.raw.count <- reads.raw.count[!names(reads.raw.count)%in%ribosomes,] 
    
    write.csv(counts(reads.raw.count), paste0("scripts/data/",filename,"_raw_read_counts.csv"))
    
    # reads.raw.count2 <- estimateSizeFactors(reads.raw.count)
    # sizeFactors(reads.raw.count2)
    reads.normalized.counts <- counts(estimateSizeFactors(reads.raw.count), normalized=TRUE)
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

