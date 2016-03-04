
# Damon Polioudakis
# 2016-02-29
# Plot % Unmapped Reads vs Mapped Reads and ERCC % vs Mapped Reads

################################################################################

rm(list=ls())
sessionInfo()

library(ggplot2)

# Picard Sequencing Statistics
picStatsDF <- read.csv("../metadata/PicardToolsQC.csv")

# Metadata
metDatDF <- read.table("../metadata/Compiled_Metadata_20160229.txt"
                       , header = TRUE, sep = "\t")

exDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon.csv", row.names = 1)

colnames(picStatsDF)

metDatDF <- metDatDF[order(as.factor(metDatDF$HTseq_IDs)), ]
metDatDF$HTseq_IDs

exDatDF <- exDatDF[ ,order(colnames(exDatDF))]
colnames(exDatDF)

picStatsDF <- picStatsDF[order(picStatsDF$X), ]
picStatsDF$X

picStatsDF <- picStatsDF[picStatsDF$X %in% metDatDF$HTseq_IDs, ]

picStatsDF[! picStatsDF$X %in% metDatDF$HTseq_IDs, ]

# pctUM <- picStatsDF$PF_READS_ALIGNED/picStatsDF$TOTAL_READS
nMapGenes <- apply(head(exDatDF, -97), 2, sum)
nMapERCC <- apply(head(tail(exDatDF, 97), -5), 2, sum)

pctERCC <- nMapERCC/picStatsDF$PF_READS_ALIGNED
round(pctERCC, 2)

plot(picStatsDF$PF_READS_ALIGNED/(10^6), pctERCC)

plot(picStatsDF$PF_READS_ALIGNED/(10^6), (1 -metDatDF$Pct_Uniquely_Mapped))
