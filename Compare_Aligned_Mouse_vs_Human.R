
# Damon Polioudakis
# 2016-02-18
# Add some RNA STAR statistics and Picard Statistics to sample metadata
# Gene Length and GC bias were calculated with Calc_Gene_Length_and_GC.R
################################################################################

rm(list=ls())

library(xlsx)
library(gridExtra)

# Load data and assign variables

graphCodeTitle <- "Compare_Aligned_Mouse_vs_Human.R"
outGraphs <- "../analysis/graphs/Compare_Aligned_Mouse_vs_Human_"

# Picard Sequencing Statistics
picStatsDF <- read.csv("../metadata/PicardToolsQC.csv")

# RNA STAR Stats Human Lane 1 
stStatsH1dF <- read.table("../metadata/RNAstar_Stats_SxaQSEQsXap096L1.txt"
                         , sep = "\t", header = TRUE)
# RNA STAR Stats Human Lane 2
stStatsH2dF <- read.table("../metadata/RNAstar_Stats_SxaQSEQsXap096L2.txt"
                         , sep = "\t", header = TRUE)
# RNA STAR Stats Mouse Lane 1 
stStatsM1dF <- read.table("../metadata/RNAstar_Stats_Mouse_SxaQSEQsXap096L1.txt"
                         , sep = "\t", header = TRUE)
# RNA STAR Stats Mouse Lane 2
stStatsM2dF <- read.table("../metadata/RNAstar_Stats_Mouse_SxaQSEQsXap096L2.txt"
                         , sep = "\t", header = TRUE)
# Reads Aligned to Both Mouse and Human
mapBothDF <- read.table("../metadata/Reads_Aligned_Human_And_Mouse_Stats.txt"
                        , header = TRUE)
# Visual QC
visualQCdF <- read.xlsx("../metadata/ExpID3-C1-HT-CaptureData-2016-02-02.xlsx", 1)

# HTseqID to SampleID
idMapDF <- read.table("../metadata/Compiled_Metadata_20160229.txt", header = TRUE)
################################################################################

### Format

## Add HTseq_SampleID to visual QC data
visualQCdF <- merge(visualQCdF, idMapDF, by.x = c("Row", "Column")
                    , by.y = c("Fluidigm_Row", "Fluidigm_Column"))

## Format Visual QC
visualQCdF$Dead.or.Alive...[is.na(visualQCdF$Dead.or.Alive...)] <- "0"
visualQCdF$Dead.or.Alive... <- gsub(".*,.*", "Multiple", visualQCdF$Dead.or.Alive...)
visualQCdF$Dead.or.Alive... <- gsub(".*0.*", "Empty", visualQCdF$Dead.or.Alive...)
visualQCdF$Dead.or.Alive... <- gsub(".*a.*", "Alive", visualQCdF$Dead.or.Alive...)
visualQCdF$Dead.or.Alive... <- gsub(".*x.*", "Dead", visualQCdF$Dead.or.Alive...)
visualQCdF$Dead.or.Alive... <- gsub("\\?", "Stained for Both", visualQCdF$Dead.or.Alive...)

## Remove 1 missing cell in Lane 2 from Lane 1
stStatsH1dF <- stStatsH1dF[stStatsH1dF$SampleID %in% stStatsH2dF$SampleID, ]
stStatsM1dF <- stStatsM1dF[stStatsM1dF$SampleID %in% stStatsH2dF$SampleID, ]
stStatsM2dF <- stStatsM2dF[stStatsM2dF$SampleID %in% stStatsH2dF$SampleID, ]
mapBothDF <- mapBothDF[mapBothDF$SampleID %in% stStatsH2dF$SampleID, ]
visualQCdF <- visualQCdF[visualQCdF$HTseq_IDs %in% stStatsH2dF$SampleID, ]
################################################################################

### Combine statistics from both lanes for mouse and for human

##
# Calculate and add total reads, percent uniquely mapped reads, percent mapped
# to Multipleple loci, and percent unmapped reads to metadata table

# Function to combine statistics from 2 lanes
Compile_STAR_Lanes <- function(stStats1dF, stStats2dF) {
  
  metDatDF <- data.frame(SampleID = stStats1dF$SampleID)

  # Total Reads
  totReads <- (stStats1dF$Number.of.input.reads
                  + stStats2dF$Number.of.input.reads)
  metDatDF$Total_Reads <- totReads
  
  # Uniquely Mapped
  uniqMapReads <- (stStats1dF$Uniquely.mapped.reads.number + stStats2dF$Uniquely.mapped.reads.number)
  metDatDF$Uniquely_Mapped <- uniqMapReads
  
  # Uniquely Mapped Percent - Convert % to numeric and Multipleply by number of reads in lane
  umReads <- (((as.numeric(sub("%", "", stStats1dF$Uniquely.mapped.reads..))/100) * stStats1dF$Number.of.input.reads)
           + (((as.numeric(sub("%", "", stStats2dF$Uniquely.mapped.reads..))/100) * stStats2dF$Number.of.input.reads)))
  # Percent Uniquely Mapped
  pctUMRds <- umReads / totReads
  metDatDF$Pct_Uniquely_Mapped <- pctUMRds
  
  # Total mapped to Multipleple loci
  totMpMtLoci <- (stStats1dF$Number.of.reads.mapped.to.multiple.loci
                + stStats2dF$Number.of.reads.mapped.to.multiple.loci)
  # Percent mapped to Multipleple loci
  pctMpMtLoci <- totMpMtLoci / totReads
  metDatDF$Pct_Mapped_Multipleple_Loci <- pctMpMtLoci
  
  # Percent unmapped reads
  Calc_Unmapped <- function (data) {
    # Convert % to numeric and Multipleply by number of reads in lane
    df <- apply(data[ ,c("X..of.reads.unmapped..too.many.mismatches"
                                        , "X..of.reads.unmapped..too.short"
                                        , "X..of.reads.unmapped..other")]
                , 2, function(column) {as.numeric(sub("%", "", column)) / 100})
    # Multipleple by number of reads in lane to convert to read number
    df <- apply(df, 2, function(column) {column * data$Number.of.input.reads})
    apply(df, 1, sum)
    }
  # Lane 1
  unMap1 <- Calc_Unmapped(stStats1dF)
  # Lane 2
  unMap2 <- Calc_Unmapped(stStats2dF)
  # Percent unmapped reads
  pctUnMap <- (unMap1 + unMap2) / totReads
  metDatDF$Pct_Unmapped <- pctUnMap
  
  # Function Output
  metDatDF
}

# Combine statistics from both lanes for human
hmStatsDF <- Compile_STAR_Lanes(stStatsH1dF, stStatsH2dF)
# And for mouse
msStatsDF <- Compile_STAR_Lanes(stStatsM1dF, stStatsM2dF)
# Add visual QC
hmStatsDF <- merge(hmStatsDF, visualQCdF[ ,c("HTseq_IDs", "Dead.or.Alive...")]
      , by.x = "SampleID", by.y = "HTseq_IDs", all.x = TRUE)
msStatsDF <- merge(msStatsDF, visualQCdF[ ,c("HTseq_IDs", "Dead.or.Alive...")]
                   , by.x = "SampleID", by.y = "HTseq_IDs", all.x = TRUE)
# Change NAs (capture sites not visually QCed) to Unkown
hmStatsDF$Dead.or.Alive...[is.na(hmStatsDF$Dead.or.Alive...)] <- "Not Visually QCed"
msStatsDF$Dead.or.Alive...[is.na(msStatsDF$Dead.or.Alive...)] <- "Not Visually QCed"
# Change to Factor and set "Not Visually QCed" as first factor order so that plots as
# black later
hmStatsDF$Dead.or.Alive... <- relevel(as.factor(hmStatsDF$Dead.or.Alive...), "Not Visually QCed")
msStatsDF$Dead.or.Alive... <- relevel(as.factor(msStatsDF$Dead.or.Alive...), "Not Visually QCed")

## Subtract number of reads that map to both human and mouse

# Combine lanes 1 and 2 for map to both statistics
sampleIDs <- mapBothDF$SampleID[1:800]
mapBothDF <- mapBothDF[1:800, ] + mapBothDF[801:1600, ]
mapBothDF$SampleID <- sampleIDs

# Substract from both Total Reads and Uniquely Mapped
Subtract_Map_Both <- function (starStats, mapBothStats) {
  outDF <- data.frame(SampleID = starStats$SampleID)
  outDF$Total_Reads <- starStats$Total_Reads - mapBothStats$Align_To_Both
  outDF$Uniquely_Mapped <- starStats$Uniquely_Mapped - mapBothStats$Align_To_Both
  outDF$Pct_Uniquely_Mapped <- outDF$Uniquely_Mapped / outDF$Total_Reads
  outDF$Dead.or.Alive... <- starStats$Dead.or.Alive...
  outDF
}
hmRmbStatsDF <- Subtract_Map_Both(hmStatsDF, mapBothDF)
msRmbStatsDF <- Subtract_Map_Both(msStatsDF, mapBothDF)

## Save adjusted reads mapping as table
write.table(hmRmbStatsDF, "../metadata/Reads_Aligned_Only_Human.txt"
            , quote = FALSE, sep = "\t", row.names = FALSE)
write.table(msRmbStatsDF, "../metadata/Reads_Aligned_Only_Mouse.txt"
            , quote = FALSE, sep = "\t", row.names = FALSE)
################################################################################

### Plots

## Histograms

# Histogram: Percentage of uniquely mapped reads - Mouse
pdf(paste0(outGraphs, "Histogram_Pct_Mapped_Mouse.pdf"))
hist(msRmbStatsDF$Pct_Uniquely_Mapped
  , main = paste(graphCodeTitle
                , "Histogram: Percentage of uniquely mapped reads - Mouse"
                , "Removed reads mapped to human"
                , sep = "\n")
  , xlab = "Percent Uniquely Mapped"
  , xlim = c(0, 1))
dev.off()

# Histogram: Percentage of uniquely mapped reads - Human
pdf(paste0(outGraphs, "Histogram_Pct_Mapped_Human.pdf"))
hist(hmRmbStatsDF$Pct_Uniquely_Mapped
     , main = paste(graphCodeTitle
                    , "Histogram: Percentage of uniquely mapped reads - Human"
                    , "Removed reads mapped to mouse"
                    , sep = "\n")
     , xlab = "Percent Uniquely Mapped"
     , xlim = c(0, 1))
dev.off()

# Histogram: Percentage of uniquely mapped reads Human + Mouse
pdf(paste0(outGraphs, "Histogram_Pct_Mapped_Human_Plus_Mouse.pdf"))
hist(hmRmbStatsDF$Pct_Uniquely_Mapped + msRmbStatsDF$Pct_Uniquely_Mapped
     , main = paste(graphCodeTitle
                    , "Histogram: Percentage of uniquely mapped reads Human + Mouse"
                    , "Removed reads that map to both human and mouse"
                    , sep = "\n")
     , xlab = "Percent Uniquely Mapped"
     , xlim = c(0, 1.2))
dev.off()

# Histogram: Ratio of Uniquely Mapped Reads - Human Compared to Mouse
pdf(paste0(outGraphs, "Histogram_Ratio_Mapped_Human_Vs_Mouse.pdf"))
hist(log(hmRmbStatsDF$Uniquely_Mapped / msRmbStatsDF$Uniquely_Mapped), 2
     , breaks = seq(-6, 6, by = 0.1)
     , main = paste(graphCodeTitle
                    , "Histogram: Ratio of Uniquely Mapped Reads - Human Compared to Mouse"
                    , "Removed reads that map to both human and mouse"
                    , sep = "\n")
     , xlab = "log2(Human # reads mapped / Mouse # reads mapped)")
dev.off()

# Histogram: Ratio of Percent Uniquely Mapped Reads - Human Compared to Mouse
pdf(paste0(outGraphs, "Histogram_Ratio_Percent_Mapped_Human_Vs_Mouse.pdf"))
hist(log(hmRmbStatsDF$Pct_Uniquely_Mapped / msRmbStatsDF$Pct_Uniquely_Mapped), 2
     , breaks = seq(-6, 6, by = 0.1)
     , main = paste(graphCodeTitle
                    , "Histogram: Ratio of Percent Uniquely Mapped Reads - Human Compared to Mouse"
                    , "Removed reads that map to both human and mouse"
                    , sep = "\n")
     , xlab = "log2(Human percent reads mapped / Mouse percent reads mapped)")
dev.off()

# Histogram: Total Reads Per Sample
pdf(paste0(outGraphs, "Total_Reads_Per_Sample.pdf"))
hist(stStatsH1dF$Number.of.input.reads + stStatsH2dF$Number.of.input.reads
     , breaks = 100
     , main = paste(graphCodeTitle
                    , "\nHistogram: Total Reads Per Sample"
                    , "\nMean:", mean(stStatsH1dF$Number.of.input.reads + stStatsH2dF$Number.of.input.reads)
                    , "\nMedian:", median(stStatsH1dF$Number.of.input.reads + stStatsH2dF$Number.of.input.reads)
                    , sep = " ")
     , xlab = "Total Reads")
dev.off()

## Number of reads uniquely mapping to human vs mouse for each capture site
# Inludes reads mapping to both mouse and human
pdf(paste0(outGraphs, "Reads_Mapping_to_Human_and_or_Mouse.pdf"))
plot(data.frame(hmStatsDF$Uniquely_Mapped, msStatsDF$Uniquely_Mapped)
     , col = hmStatsDF$Dead.or.Alive...
     , cex = 0.75
     , xlim = c(0, 5*10^5)
     , ylim = c(0, 5*10^5)
     , xlab = "Reads Uniquely Mapped to Human"
     , ylab = "Reads Uniquely Mapped to Mouse"
     , main = paste0(graphCodeTitle
                     , "\nUniquely Mapped Reads to Human or Mouse"
                     , "\nIncludes Reads Mapping to Both Human and Mouse"))
points(data.frame(hmStatsDF$Uniquely_Mapped[hmStatsDF$Dead.or.Alive... == "Empty"]
                  , msStatsDF$Uniquely_Mapped[hmStatsDF$Dead.or.Alive... == "Empty"])
       , col = "blue", pch = 3)
legend('topright', legend = paste(levels(hmStatsDF$Dead.or.Alive...), "-"
                                  , table(hmStatsDF$Dead.or.Alive...))
       , pch = 1, col = seq_along(levels(hmStatsDF$Dead.or.Alive...)))
dev.off()

# Reads mapping to both mouse and human removed
pdf(paste0(outGraphs, "Reads_Mapping_to_Human_vs_Mouse.pdf"))
plot(data.frame(hmRmbStatsDF$Uniquely_Mapped, msRmbStatsDF$Uniquely_Mapped)
     , col = hmRmbStatsDF$Dead.or.Alive...
     , cex = 0.75
     , xlim = c(0, 5*10^5)
     , ylim = c(0, 5*10^5)
     , xlab = "Reads Uniquely Mapped to Human"
     , ylab = "Reads Uniquely Mapped to Mouse"
     , main = paste0(graphCodeTitle
                     , "\nUniquely Mapped Reads to Human or Mouse"
                     , "\nReads Mapping to Both Human and Mouse Removed"))
pal = palette()[-1]
points(data.frame(hmRmbStatsDF$Uniquely_Mapped[! hmRmbStatsDF$Dead.or.Alive... == "Not Visually QCed"]
                  , msRmbStatsDF$Uniquely_Mapped[! hmRmbStatsDF$Dead.or.Alive... == "Not Visually QCed"])
       , pch = 3, col = hmRmbStatsDF$Dead.or.Alive...[! hmRmbStatsDF$Dead.or.Alive... == "Not Visually QCed"], pal = palette()[-1])
abline(v = 5*10^4)
abline(h = 5*10^4)
legend('topright', legend = paste(levels(hmRmbStatsDF$Dead.or.Alive...), "-"
                                  , table(hmRmbStatsDF$Dead.or.Alive...))
       , pch = 1, col = seq_along(levels(hmRmbStatsDF$Dead.or.Alive...)))
dev.off()

# Reads mapping to both mouse and human removed - highlight Visually QCed Empty
pdf(paste0(outGraphs, "Reads_Mapping_to_Human_vs_Mouse_VisualQCempty.pdf"))
plot(data.frame(hmRmbStatsDF$Uniquely_Mapped, msRmbStatsDF$Uniquely_Mapped)
     , col = hmRmbStatsDF$Dead.or.Alive...
     , cex = 0.75
     , xlim = c(0, 5*10^5)
     , ylim = c(0, 5*10^5)
     , xlab = "Reads Uniquely Mapped to Human"
     , ylab = "Reads Uniquely Mapped to Mouse"
     , main = paste0(graphCodeTitle
                     , "\nUniquely Mapped Reads to Human or Mouse"
                     , "\nReads Mapping to Both Human and Mouse Removed"))
points(data.frame(hmRmbStatsDF$Uniquely_Mapped[hmRmbStatsDF$Dead.or.Alive... == "Empty"]
                  , msRmbStatsDF$Uniquely_Mapped[hmRmbStatsDF$Dead.or.Alive... == "Empty"])
       , col = "blue", pch = 3)
legend('topright', legend = paste(levels(hmRmbStatsDF$Dead.or.Alive...), "-"
                                  , table(hmRmbStatsDF$Dead.or.Alive...))
       , pch = 1, col = seq_along(levels(hmRmbStatsDF$Dead.or.Alive...)))
dev.off()

## Table of number of doublets, single human cells, single mouse cells
doubletDF <- data.frame(
  Human = nrow(hmRmbStatsDF[hmRmbStatsDF$Uniquely_Mapped > 5*10^4 & msRmbStatsDF$Uniquely_Mapped < 5*10^4, ])
  , Mouse = nrow(hmRmbStatsDF[hmRmbStatsDF$Uniquely_Mapped < 5*10^4 & msRmbStatsDF$Uniquely_Mapped > 5*10^4, ])
  , Doublet = nrow(hmRmbStatsDF[hmRmbStatsDF$Uniquely_Mapped > 5*10^4 & msRmbStatsDF$Uniquely_Mapped > 5*10^4, ])
  , Low_Read_Depth = nrow(hmRmbStatsDF[hmRmbStatsDF$Uniquely_Mapped < 5*10^4 & msRmbStatsDF$Uniquely_Mapped < 5*10^4, ])
)
write.table(doubletDF, "../analysis/tables/Compare_Aligned_Mouse_vs_Human_Doublets.txt"
            , quote = FALSE, sep = "\t", row.names = FALSE)

# Percent mapping to mouse and to human vs capture sites ordered by percent
# mapping to human
pdf(paste0(outGraphs, "Percent_Reads_Mapping_Ordered_By_Human.pdf"))
par(cex = 0.75)
plotDF <- hmRmbStatsDF[order(hmRmbStatsDF$Pct_Uniquely_Mapped), ]
plot(sort(plotDF$Pct_Uniquely_Mapped)
     , pch = 1
     , col = plotDF$Dead.or.Alive...
     , ylab = "Percent Mapping"
     , xlab = "Capture Sites Sorted by Percent Mapping to Human"
     , main = paste0(graphCodeTitle
                     , "\nPercent Mapping to Mouse and to Human"))
points(data.frame(msRmbStatsDF$Pct_Uniquely_Mapped)[order(hmRmbStatsDF$Pct_Uniquely_Mapped), ]
       , pch = 2
       , col = plotDF$Dead.or.Alive...)
legend('topright', legend = c(levels(plotDF$Dead.or.Alive...), "Human", "Mouse")
       , pch = c(rep(1, 7), 2), col = c(seq_along(levels(plotDF$Dead.or.Alive...)), 1, 1))
dev.off()

# Percent mapping to mouse and to human vs capture sites ordered by percent
# mapping to mouse
pdf(paste0(outGraphs, "Percent_Reads_Mapping_Ordered_By_Mouse.pdf"))
par(cex = 0.75)
plotDF <- msRmbStatsDF[order(msRmbStatsDF$Pct_Uniquely_Mapped), ]
plot(sort(plotDF$Pct_Uniquely_Mapped)
     , pch = 2
     , col = plotDF$Dead.or.Alive...
     , ylab = "Percent Mapping"
     , xlab = "Capture Sites Sorted by Percent Mapping to Mouse"
     , main = paste0(graphCodeTitle
                     , "\nPercent Mapping to Mouse and to Human"))
points(data.frame(hmRmbStatsDF$Pct_Uniquely_Mapped)[order(msRmbStatsDF$Pct_Uniquely_Mapped), ]
       , pch = 1
       , col = plotDF$Dead.or.Alive...)
legend('topright', legend = c(levels(hmRmbStatsDF$Dead.or.Alive...), "Human", "Mouse")
       , pch = c(rep(1, 7), 2), col = c(seq_along(levels(hmRmbStatsDF$Dead.or.Alive...)), 1, 1))
dev.off()

# Percent uniquely mapping to human vs mouse
pdf(paste0(outGraphs, "Percent_Uniquely_Mapping_To_Human_Vs_Mouse.pdf"))
plot(hmRmbStatsDF$Pct_Uniquely_Mapped, msRmbStatsDF$Pct_Uniquely_Mapped
     , col = msRmbStatsDF$Dead.or.Alive...
     , ylab = "Percent Uniquely Mapped - Human"
     , xlab = "Percent Uniquely Mapped - Mouse"
     , main = paste0(graphCodeTitle
                     , "\nPercent Uniquely Mapping to Human Vs Mouse"))
legend('topright', legend = paste(levels(hmRmbStatsDF$Dead.or.Alive...), "-"
                                  , table(hmRmbStatsDF$Dead.or.Alive...))
       , pch = 1, col = seq_along(levels(hmRmbStatsDF$Dead.or.Alive...)))
dev.off()
################################################################################


library(ggplot2)

# Input transcript expression table (counts from HTseq)
exDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon.csv", row.names = 1)

# Remove ERCCs

exDatDF <- head(exDatDF, -97)


## Functions for MDS analysis

# Function to output data frame of MDS PCs values and PCs
calcMDS <- function (exprDF) {
  # dist calculates the distances between the rows of a data matrix
  # Transpose data frame so samples are rows and genes are columns
  mds = cmdscale(dist(t(exprDF)), eig = T)
  pc1 = mds$eig[1]^2 / sum(mds$eig^2)
  pc2 = mds$eig[2]^2 / sum(mds$eig^2)
  mdsAndTreatmentLDF <- data.frame(mds$points
                                   , pc1 = pc1, pc2 = pc2)
  mdsAndTreatmentLDF
}

## MDS All Genes

# Calculate MDS
mdsDF <- calcMDS(exDatDF)
# Add column with sample type info
mdsDF$VisualQC <- head(hmRmbStatsDF$Dead.or.Alive..., -1)

hmRmbStatsDF <- hmRmbStatsDF[hmRmbStatsDF$SampleID %in% row.names(mdsDF), ]
msRmbStatsDF <- msRmbStatsDF[msRmbStatsDF$SampleID %in% row.names(mdsDF), ]

mdsDF$Reads_Align_Filter <- "Not Enough Reads"
mdsDF$Reads_Align_Filter[hmRmbStatsDF$Uniquely_Mapped < 5*10^4
                         & msRmbStatsDF$Uniquely_Mapped > 5*10^4] <- "Mouse"
mdsDF$Reads_Align_Filter[hmRmbStatsDF$Uniquely_Mapped > 5*10^4
                         & msRmbStatsDF$Uniquely_Mapped < 5*10^4] <- "Human"
mdsDF$Reads_Align_Filter[hmRmbStatsDF$Uniquely_Mapped > 5*10^4
                         & msRmbStatsDF$Uniquely_Mapped > 5*10^4] <- "Mouse and Human"

ftmdsDF <- mdsDF[hmRmbStatsDF$Uniquely_Mapped > 5*10^4
  & msRmbStatsDF$Uniquely_Mapped < 5*10^4, ]

ftmdsDF <- mdsDF[msRmbStatsDF$Uniquely_Mapped < 5*10^4, ]

ggplot(mdsDF, aes(x = X1, y = X2, color = factor(VisualQC))) +
  geom_point(size = 1)

ggplot(ftmdsDF, aes(x = X1, y = X2, color = factor(VisualQC))) +
  geom_point(size = 1) +
  geom_point(size = 1) +
  xlab(paste("PC1 (", signif(100*mdsDF$pc1, 3), "%)", sep = "")) +
  ylab(paste("PC2 (", signif(100*mdsDF$pc2, 3), "%)", sep = "")) +
  labs(title = paste(graphCodeTitle, "MDS Plot: HTSeq log2(Counts + 1)"
                     , sep = "\n")) +
  theme_grey(base_size = 16)

ggplot(mdsDF, aes(x = X1, y = X2, color = factor(Reads_Align_Filter))) +
  geom_point(size = 1) +
  xlab(paste("PC1 (", signif(100*mdsDF$pc1, 3), "%)", sep = "")) +
  ylab(paste("PC2 (", signif(100*mdsDF$pc2, 3), "%)", sep = "")) +
  labs(title = paste(graphCodeTitle, "MDS Plot: HTSeq log2(Counts + 1)"
                     , sep = "\n")) +
  theme_grey(base_size = 16)


+
  geom_text(aes(label = row.names(mdsDF)), vjust = -1, size = 1.5) +
  scale_color_discrete(name = "Sample Type") +
  xlab(paste("PC1 (", signif(100*mdsDF$pc1, 3), "%)", sep = "")) +
  ylab(paste("PC2 (", signif(100*mdsDF$pc2, 3), "%)", sep = "")) +
  labs(title = paste(graphCodeTitle, "MDS Plot: HTSeq log2(Counts + 1)"
                     , sep = "\n")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(aspect.ratio = 4/5)












plot(hmRmbStatsDF$Pct_Uniquely_Mapped[hmRmbStatsDF$Total_Reads > 100000], msRmbStatsDF$Pct_Uniquely_Mapped[msRmbStatsDF$Total_Reads > 100000])


mean(stStatsH1dF$Number.of.input.reads + stStatsH2dF$Number.of.input.reads)
median(stStatsH1dF$Number.of.input.reads + stStatsH2dF$Number.of.input.reads)

mean(stStatsM1dF$Number.of.input.reads + stStatsM2dF$Number.of.input.reads)
median(stStatsM1dF$Number.of.input.reads + stStatsM2dF$Number.of.input.reads)


df <- merge(hmRmbStatsDF, idMapDF, by.x = "SampleID", by.y = "HTseq_IDs")
plot(sort(df$Pct_Uniquely_Mapped.x), col = df$Fluidigm_Row)
