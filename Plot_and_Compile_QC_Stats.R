# Damon Polioudakis
# 2016-02-26
# Process QC statistics from RNAseq
# TODO
#   gcbias section is incomplete
################################################################################

rm(list=ls())
sessionInfo()
library(ggplot2)
library(reshape2)

# Load data and assign variables

# alignment_stats.txt
picAlign <- read.table("../data/QC/merged/alignment_stats.txt"
                       , row.names = 1, skip = 1)
# remove unused column names from alignment_stats.txt
cNames <- unlist(strsplit(readLines(con = "../data/QC/merged/alignment_stats.txt"
                             , n = 1), split = " "))
colnames(picAlign) <- cNames[ -c(1, (length(cNames)-2):length(cNames))]

# rnaseq_stats.txt
picSeq <- read.table("../data/QC/merged/rnaseq_stats.txt", skip = 1, row.names = 1)
# remove unused column names from alignment_stats.txt
# (Picard outputs column names for statistics even if they are not calculated)
cNames <- unlist(strsplit(readLines(con = "../data/QC/merged/rnaseq_stats.txt"
                                    , n = 1), split = " "))
colnames(picSeq) <- cNames[!cNames %in% c("RIBOSOMAL_BASES"
                  , "PCT_RIBOSOMAL_BASES", "SAMPLE", "LIBRARY", "READ_GROUP")]

# rnaseq_stats_Transcript_Coverage.txt
# Add NA to blanks - I think these are from samples with almost no coverage
txCovDF <- read.table("../data/QC/merged/rnaseq_stats_Transcript_Coverage.txt"
                      , fill = TRUE)

# gcbias_summary.txt
gcBiasDF <- read.table("../data/QC/merged/gcbias_summary.txt", header = TRUE)

# duplication_stats.txt
# Add NA to blanks
dupDF <- read.table("../data/QC/merged/duplication_stats.txt"
                    , fill = TRUE, header = TRUE)
# Shift columns so they align with headers
dupDF[ ,c(3:ncol(dupDF))] <- dupDF[ ,c(4:ncol(dupDF))]

# Reads aligned for to human only (excluded those that also map to mouse)
nMapHumanDF <- read.table("../metadata/Reads_Aligned_Only_Human.txt"
                          , header = TRUE, sep = "\t")
################################################################################

# RNAseq Percentage of Bases Aligned by Location

pctHq <- picAlign[,"PF_HQ_ALIGNED_BASES"]/picAlign[,"PF_ALIGNED_BASES"]

pdf("../analysis/graphs/Plot_QC_Percent_Location.pdf")
boxplot(cbind(pctHq, picSeq[,c(14, 10:13)]), ylab = "Percent of Bases"
        , main = "Plot_and_Compile_QC_Stats.R
RNA-seq QC metrics across samples"
        , names = c("High\nQuality"
                    , "mRNA\nBases"
                    , "Protein\nCoding"
                    , "Untranslated\nRegion"
                    , "Intronic\nRegion"
                    , "Intergenic\nRegion")
        , ylim = c(0, 1))
dev.off()
################################################################################

# Boxplot number of fragments mapping to different types of location

ggDF <- picSeq[, c("PF_BASES", "CODING_BASES", "INTRONIC_BASES"
           , "INTERGENIC_BASES")]/100
colnames(ggDF) <- c("Total", "Coding", "Intronic", "Intergenic")
ggDF <- melt(ggDF)

ggplot(ggDF, aes(x = variable, y = value, col = variable)) +
  geom_boxplot(outlier.shape = NA) +
  # Adjust limits after outlier removal
  coord_cartesian(ylim = range(boxplot(ggDF[ggDF$variable == "Total", ]$value, plot = FALSE)$stats) * c(.9, 1.1)) +
  theme_bw(base_size = 18) +
  guides(col = FALSE) +
  ylab("Number of Fragments Aligned") +
  xlab("Location") +
  ggtitle("Plot_and_Compile_QC_Stats.R
Number of Fragments Aligned by Location")
ggsave("../analysis/graphs/Plot_QC_Aligned_Location.pdf")
################################################################################

# Transcript Coverage

keep <- match(unique(txCovDF[ ,1]), txCovDF[ ,1])
txCovDF <- txCovDF[keep, ]
rownames(txCovDF) <- txCovDF[ ,1]
txCovDF.cols <- as.character(txCovDF[1, seq(4, 204, by = 2)])
txCovDF <- txCovDF[ ,seq(5, 205, by = 2)]
colnames(txCovDF) <- txCovDF.cols
txCovqtlDF <- apply(txCovDF, 2, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)

pdf("../analysis/graphs/Plot_QC_Transcript_Coverage.pdf")
plot(x = as.numeric(colnames(txCovqtlDF)), y = txCovqtlDF[2,]
     , xlab = "Percentile of gene body (5' -> 3')"
     , ylab = "Coverage Relative to Whole Transcript"
     , pch = 19, ylim = c(0, 7)
     , main = "Plot_and_Compile_QC_Stats.R
Relative transcript coverage - median with 95% CIs across samples")
lines(x = as.numeric(colnames(txCovqtlDF)), y = txCovqtlDF[2,], col = "black")
lines(x = as.numeric(colnames(txCovqtlDF)), y = txCovqtlDF[1,], col = "grey")
lines(x = as.numeric(colnames(txCovqtlDF)), y = txCovqtlDF[3,], col = "grey")
dev.off()
################################################################################

# GC Bias

gcBias <- matrix(NA, nrow = nrow(gcBiasDF), ncol = 101)
colnames(gcBias) <- as.character(gcBiasDF[1, seq(2, 602, by = 6)])
rownames(gcBias) <- gcBiasDF[ ,1]
gcBias <- gcBiasDF[ ,seq(6, 607, by = 6)]
gc.coverage.quant <- apply(gcBias,2,quantile,c(0.025,0.5,0.975))
colnames(gc.coverage.quant) <- as.character(gcBiasDF[1,seq(2,602,by=6)])

gcBiasqual <- matrix(NA,nrow=nrow(gcBiasDF),ncol=101)
colnames(gcBiasqual) <- gcBiasDF[1,seq(2,602,by=6)]
rownames(gcBiasqual) <- gcBiasDF[,1]
gcBiasqual <- gcBiasDF[,seq(5,607,by=6)]
gc.qual.quant <- apply(gcBiasqual,2,quantile,c(0.025,0.5,0.975))

gcBiasbin <- matrix(NA,nrow=nrow(gcBiasDF),ncol=101)
colnames(gcBiasbin) <- gcBiasDF[1,seq(2,602,by=6)]
rownames(gcBiasbin) <- gcBiasDF[,1]
gcBiasbin <- gcBiasDF[,seq(3,607,by=6)]
gc.bins.quant <- apply(gcBiasbin,2,quantile,c(0.5))
gc.bins.quant <- gc.bins.quant/sum(gc.bins.quant)

par(mfrow=c(3,1))
plot(x=as.numeric(colnames(gc.coverage.quant)),y=gc.coverage.quant[2,],xlab="% GC content",ylab="Normalized coverage",pch=19,ylim=c(0,15),main="Coverage by GC percentage - median with 95% CIs across samples")
lines(x=as.numeric(colnames(gc.coverage.quant)),y=gc.coverage.quant[2,],col="black")
lines(x=as.numeric(colnames(gc.coverage.quant)),y=gc.coverage.quant[1,],col="grey")
lines(x=as.numeric(colnames(gc.coverage.quant)),y=gc.coverage.quant[3,],col="grey")

plot(x=as.numeric(colnames(gc.coverage.quant)),y=gc.qual.quant[2,],xlab="%GC content",ylab="Read quality score",pch=19,ylim=c(0,36),main="Read quality by GC percent - median with 95% CIs across samples")
lines(x=as.numeric(colnames(gc.coverage.quant)),y=gc.qual.quant[2,],col="black")
lines(x=as.numeric(colnames(gc.coverage.quant)),y=gc.qual.quant[1,],col="grey")
lines(x=as.numeric(colnames(gc.coverage.quant)),y=gc.qual.quant[3,],col="grey")

plot(density(x=as.numeric(colnames(gc.coverage.quant)),y=gc.bins.quant),xlab="%GC content",ylab="Proportion of 100bp bins",pch=19,ylim=c(0,max(gc.bins.quant)+0.01),main="Proportion of bins corresponding to each GC percentile")
################################################################################

### Duplication Metrics

## Percent Duplicate for all capture sites

# For Single End RNAseq
ggDF <- data.frame(PERCENT_DUPLICATION = dupDF$PERCENT_DUPLICATION
                   , SAMPLE = rownames(dupDF))
# Order by percent duplication
ggDF <- ggDF[order(ggDF$PERCENT_DUPLICATION), ]
ggDF$SAMPLE <- factor(ggDF$SAMPLE, levels = ggDF$SAMPLE)
mPctDup <- mean(ggDF$PERCENT_DUPLICATION)

ggplot(ggDF, aes(x = SAMPLE, y = PERCENT_DUPLICATION)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Percentage Duplication") +
  xlab("Samples") +
  ggtitle(paste0("Plot_and_Compile_QC_Stats.R"
                , "\nPercentage of Mapped Sequence Marked as Duplicate"
                , "\nMean Percent Duplicate: ", mPctDup))
ggsave("../analysis/graphs/Plot_QC_Percent_Duplication.pdf")


## Percent Duplicate for capture sites with > 50,000 mapped reads that align
## to only human genome

# Collect duplication statistics
ggDF <- data.frame(PERCENT_DUPLICATION = dupDF$PERCENT_DUPLICATION
                   , SAMPLE = dupDF$SAMPLE)
# Add number of reads aligned to human only
ggDF <- merge(ggDF, nMapHumanDF[ ,c("SampleID", "Uniquely_Mapped")]
      , by.x = "SAMPLE", by.y = "SampleID")
# Filter for samples with > 50,000 mapped reads
ggDF <- ggDF[ggDF$Uniquely_Mapped > 5*10^4, ]

# Order by percent duplication
ggDF <- ggDF[order(ggDF$PERCENT_DUPLICATION), ]
ggDF$SAMPLE <- factor(ggDF$SAMPLE, levels = ggDF$SAMPLE)
mPctDup <- mean(ggDF$PERCENT_DUPLICATION)

ggplot(ggDF, aes(x = SAMPLE, y = PERCENT_DUPLICATION)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Percentage Duplication") +
  xlab("Samples") +
  ggtitle(paste0("Plot_and_Compile_QC_Stats.R"
                 , "\nPercentage of Mapped Sequence Marked as Duplicate"
                 , "\nFiltered for capture sites with > 50,000 reads mapped to only human"
                 , "\nMean Percent Duplicate: ", mPctDup))
ggsave("../analysis/graphs/Plot_QC_Percent_Duplication_Ft50000mapHuman.pdf")
################################################################################

# Compile QC metrics

## Compile most important QC metrics - choose the metrics that you think are valuable to the analysis and to sequencing statistics
## Metrics concering read depth (ie. Total Reads) and GC bias (ie. 5'/3' bias) should always be included

PQCdat <- cbind(picAlign[ ,c(2, 6, 9, 21)], picSeq, GC_DROPOUT = gcBiasDF[ ,c(6)]
                , dupDF[ ,c(7, 8)])
write.csv(PQCdat,"../metadata/PicardToolsQC.csv")
dev.off()
