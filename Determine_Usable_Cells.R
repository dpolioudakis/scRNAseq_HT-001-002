
# Damon Polioudakis
# 2016-02-29
# Plot % Unmapped Reads vs Mapped Reads and ERCC % vs Mapped Reads

################################################################################

rm(list=ls())
sessionInfo()

library(ggplot2)

outGraphs <- "../analysis/graphs/Determine_Usable_Cells_"
graphCodeTitle <- "Determine_Usable_Cells.R"

# Picard Sequencing Statistics
picStatsDF <- read.csv("../metadata/PicardToolsQC.csv")

# Metadata
metDatDF <- read.table("../metadata/Compiled_Metadata_20160229.txt"
                       , header = TRUE, sep = "\t")

exDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon.csv", row.names = 1)
################################################################################

## Format

metDatDF <- metDatDF[order(as.factor(metDatDF$HTseq_IDs)), ]
metDatDF$HTseq_IDs

exDatDF <- exDatDF[ ,order(colnames(exDatDF))]

picStatsDF <- picStatsDF[order(picStatsDF$X), ]

picStatsDF <- picStatsDF[picStatsDF$X %in% metDatDF$HTseq_IDs, ]

picStatsDF[! picStatsDF$X %in% metDatDF$HTseq_IDs, ]


## Sum ERCCs for each sample
nMapERCC <- apply(head(tail(exDatDF, 97), -5), 2, sum)


## Plots

# Histogram: Percentage of uniquely mapped reads - Mouse
pdf(paste0(outGraphs, "Bargraph_ERCC_Counts.pdf"))
barplot(nMapERCC
     , main = paste(graphCodeTitle
                    , "Raw ERCC counts for each capture site"
                    , sep = "\n")
     , xlab = "Capture Sites"
     , ylab = "Counts")
dev.off()

pdf(paste0(outGraphs, "Histogram_ERCC_Counts.pdf"))
hist(nMapERCC, breaks = 0:130
     , main = paste(graphCodeTitle
                    , "Histogram: ERCC Counts Per Sample"
                    , sep = "\n")
     , xlab = "ERCC Counts Per Sample")
dev.off()

