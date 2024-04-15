#DF genes
# Evol between ancester and evolved lineages
source("scripts/functions/RNAseq_analyses.R")
#WGCNA setup####################################################################
################################################################################
Tpadj <- 0.2  #Threshold padj 
################################################################################
#MGGP01
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
#Ancesters vs Evolved
res <- as.data.frame(results(dds01, contrast = list(c("PopA01_a_23","Pop01A_b_23_n","PopA01_a_17","Pop01A_b_17_n"),
                                                    c("Pop01F_a_23_n","Pop01F_a_17_n","Pop01F_b_23","Pop01F_b_17","Pop01F_c_23","Pop01F_c_17",
                                                      "Pop01S17_a_23_n","Pop01S17_a_17_n", "Pop01S17_b_23","Pop01S17_b_17", "Pop01S17_c_23_n","Pop01S17_c_17_n",
                                                      "Pop01S23_a_23","Pop01S23_a_17","Pop01S23_b_23_n","Pop01S23_b_17_n","Pop01S23_c_23_n","Pop01S23_c_17_n")), listValues=c(1/4, -1/18)))

plot(res$log2FoldChange, -log10(res$padj), col = ifelse(res$padj < Tpadj & abs(res$log2FoldChange) >2,'blue','black'), pch = 19, main="MGGP01 Anc vs evolved" )
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("red","black"),
       legend=c( "Significant P-value","Not significant") )
text(5, 60, paste0(nrow(subset(res, padj < Tpadj & abs(log2FoldChange) >2 ))," DEG"))

DEG_01_anc <- rownames(subset(res, padj < Tpadj & abs(log2FoldChange) >2 ))

#Temperature 17°C vs 23°C
res1 <- as.data.frame(results(dds01, contrast = list(c("Pop01F_a_17_n","PopA01_a_17","Pop01A_b_17_n","Pop01F_b_17",
                                                       "Pop01F_c_17","Pop01S17_a_17_n","Pop01S17_b_17","Pop01S17_c_17_n",
                                                       "Pop01S23_a_17","Pop01S23_b_17_n","Pop01S23_c_17_n"),
                                                     c("PopA01_a_23","Pop01A_b_23_n", "Pop01F_a_23_n","Pop01F_b_23","Pop01F_c_23",
                                                       "Pop01S17_a_23_n", "Pop01S17_b_23", "Pop01S17_c_23_n",
                                                       "Pop01S23_a_23","Pop01S23_b_23_n","Pop01S23_c_23_n")), listValues=c(1/11, -1/11)))

plot(res1$log2FoldChange, -log10(res1$padj), col = ifelse(res1$padj < Tpadj & abs(res1$log2FoldChange) > 2,'red3','black'), pch = 19, main="MGGP01 17 vs 23" )
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("red3","black"),
       legend=c( "Significant P-value","Not significant") )
text(-10, 60, paste0(nrow(subset(res1, padj < Tpadj & abs(log2FoldChange) >2 ))," DEG"))

#Very few repeatable DEG.
DEG_01_temp <- rownames(subset(res1, padj < Tpadj & abs(log2FoldChange) >2 ))


################################################################################
#MGGP44
################################################################################
samp.info <- read.csv("All_Data/samples_design.csv")
samp.info <- samp.info[order(as.character(samp.info$sample_id), method = c("radix")),]
samp.info <- subset(samp.info, ancester==44)
samp.info[,2:ncol(samp.info)] <- lapply(samp.info[,2:ncol(samp.info)] , factor)
#Counts
raw_counts <- read.csv("scripts/data/counts_44_raw_corrected.csv")
rownames(raw_counts) <- raw_counts$X
raw_counts <- raw_counts[,2:ncol(raw_counts)]
dds44 <- DESeq(DESeqDataSetFromMatrix(raw_counts, colData = samp.info, design = ~ 0 + Pop))
################################################################################
#Ancesters vs Evolved
res3 <- as.data.frame(results(dds44, contrast = list(c("Pop44A_a_17","Pop44A_a_23","Pop44A_b_23_n"),
                                                     c("Pop44F_a_17","Pop44F_a_23","Pop44F_b_17_n","Pop44F_b_23_n","Pop44F_c_17",    
                                                       "Pop44F_c_23","Pop44S17_a_17","Pop44S17_a_23","Pop44S17_b_17_n","Pop44S17_b_23_n",
                                                       "Pop44S23_a_17_n","Pop44S23_a_23_n","Pop44S23_b_17_n",
                                                       "Pop44S23_b_23_n","Pop44S23_c_17","Pop44S23_c_23")), listValues=c(1/3, -1/16)))

plot(res3$log2FoldChange, -log10(res3$padj), col = ifelse(res3$padj < Tpadj & abs(res3$log2FoldChange) >2,'blue','black'), pch = 19, main="MGGP44 Anc vs evolved" )
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("red","black"),
       legend=c( "Significant P-value","Not significant") )
text(-5, 30, paste0(nrow(subset(res3, padj < Tpadj & abs(log2FoldChange) >2 ))," DEG"))

DEG_44_anc <- rownames(subset(res3, padj < Tpadj & abs(log2FoldChange) >2 ))

#Temperature 17°C vs 23°C
res4 <- as.data.frame(results(dds44, contrast = list(c("Pop44A_a_17","Pop44F_a_17","Pop44F_b_17_n","Pop44F_c_17","Pop44S17_a_17","Pop44S17_b_17_n","Pop44S23_a_17_n",
                                                       "Pop44S23_b_17_n","Pop44S23_c_17"),
                                                     c("Pop44A_b_23_n","Pop44A_b_23_n","Pop44F_a_23","Pop44F_b_23_n","Pop44F_c_23","Pop44S17_a_23","Pop44S17_b_23_n","Pop44S23_a_23_n",
                                                       "Pop44S23_b_23_n","Pop44S23_c_23")), listValues=c(1/9, -1/9)))

plot(res4$log2FoldChange, -log10(res4$padj), col = ifelse(res4$padj < Tpadj & abs(res4$log2FoldChange) > 2,'red3','black'), pch = 19, main="MGGP44 17 vs 23" )
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("red3","black"),
       legend=c( "Significant P-value","Not significant") )
text(-8, 50, paste0(nrow(subset(res4, padj < Tpadj & abs(log2FoldChange) >2 ))," DEG"))
#
#Very few repeatable DEG.
DEG_44_temp <- rownames(subset(res4, padj < Tpadj & abs(log2FoldChange) >2 ))


##
pdf("figures/Volcano_plot.pdf", width=10, height=4.5)
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
plot(res$log2FoldChange, -log10(res$padj), col = ifelse(res$padj < Tpadj & abs(res$log2FoldChange) >2,'green3','black'), pch = 19, main="MGGP01 Anc vs evolved", cex=0.6 )
legend("topright", box.lty=0,  bg="transparent", fill=c("green3","black"),
       legend=c( "Significant P-value","Not significant") )
text(0, 60, paste0("Padj <",Tpadj,": ", nrow(subset(res, padj < Tpadj & abs(log2FoldChange) >2 ))," DEG"))
#
plot(res1$log2FoldChange, -log10(res1$padj), col = ifelse(res1$padj < Tpadj & abs(res1$log2FoldChange) > 2,'green3','black'), pch = 19, main="MGGP01 17 vs 23", cex=0.6 )
text(-7, 60, paste0("Padj <",Tpadj,": ",nrow(subset(res1, padj < Tpadj & abs(log2FoldChange) >2 ))," DEG"))
#
plot(res3$log2FoldChange, -log10(res3$padj), col = ifelse(res3$padj < Tpadj & abs(res3$log2FoldChange) >2,'green3','black'), pch = 19, main="MGGP44 Anc vs evolved", cex=0.6 )
legend("topright", box.lty=0,  bg="transparent", fill=c("green3","black"),
       legend=c( "Significant P-value","Not significant") )
text(-5, 25, paste0("Padj <",Tpadj,": ",nrow(subset(res3, padj < Tpadj & abs(log2FoldChange) >2 ))," DEG"))
#
plot(res4$log2FoldChange, -log10(res4$padj), col = ifelse(res4$padj < Tpadj & abs(res4$log2FoldChange) > 2,'green3','black'), pch = 19, main="MGGP44 17 vs 23", cex=0.6 )
text(-5, 50, paste0("Padj <",Tpadj,": ",nrow(subset(res4, padj < Tpadj & abs(log2FoldChange) >2 ))," DEG"))
#
#barplot(as.matrix(Df.bar4), col = c("yellowgreen"),las=2, main="DEG to Temp (for any pop)", ylim = c(0, max(Df.bar3[1,])))
dev.off()


length(c(DEG_44_anc, DEG_01_anc))- length(unique(c(DEG_44_anc, DEG_01_anc)))
length(c(DEG_44_temp, DEG_01_temp))- length(unique(c(DEG_44_temp, DEG_01_temp)))

DEG_44_anc[ DEG_44_anc %in% DEG_01_anc]
DEG_44_temp[ DEG_44_temp %in% DEG_01_temp]

################################################################################
#GO terms
annot <- read.csv("All_Data/Annotations-metadata.csv", sep="\t")
subset(annot, transcript %in% DEG_01_temp)$Seq..Description
subset(annot, transcript %in% DEG_44_temp)$Seq..Description

subset(annot, transcript %in% DEG_44_anc[ DEG_44_anc %in% DEG_01_anc])$Seq..Description
subset(annot, transcript %in% DEG_44_temp[ DEG_44_temp %in% DEG_01_temp])$Seq..Description



