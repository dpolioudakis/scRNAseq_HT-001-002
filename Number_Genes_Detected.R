# Damon Polioudakis
# 2016-02-29
# Graph of Number of Genes expressed
################################################################################

rm(list=ls())
sessionInfo()

library(ggplot2)
library(xlsx)
library(reshape2)

# Input transcript expression table (counts from HTseq) - Human
exDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon.csv", row.names = 1)

# Input transcript expression table (counts from HTseq) - Mouse
exDatMmDF <- read.csv("../data/htseq/mouse/merged/Exprs_HTSCexon.csv", row.names = 1)

# Metadata
metDatDF <- read.table("../metadata/Compiled_Metadata_20160317.txt"
                       , header = TRUE, sep = "\t")

metDatMmDF <- read.table("../metadata/Compiled_Metadata_Mouse_20160317.txt"
                       , header = TRUE, sep = "\t")

# Reads aligned to human only (excluded those that also map to mouse)
nMapHumanDF <- read.table("../metadata/Reads_Aligned_Only_Human.txt"
                          , header = TRUE, sep = "\t")
nMapMouseDF <- read.table("../metadata/Reads_Aligned_Only_Mouse.txt"
                          , header = TRUE, sep = "\t")

# Visual QC
visualQCdF <- read.xlsx("../metadata/ExpID3-C1-HT-CaptureData-2016-02-02.xlsx", 1)

# HTseqID to SampleID
idMapDF <- read.table("../metadata/Compiled_Metadata_20160229.txt", header = TRUE)

# Picard Tools QC statistics
picardDF <- read.csv("../metadata/PicardToolsQC_Mouse.csv")

dir.create("../analysis/tables/", recursive = TRUE)

# ggplot2 Theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

### Format

## Move ERCCs to separate data frame

erccDF <- tail(exDatDF, 97)
exDatDF <- head(exDatDF, -97)

erccMmDF <- tail(exDatMmDF, 97)
exDatMmDF <- head(exDatMmDF, -97)


## Filter for missing samples

nMapMouseDF <- nMapMouseDF[(nMapMouseDF$SampleID %in% colnames(exDatDF)), ]
nMapHumanDF <- nMapHumanDF[(nMapHumanDF$SampleID %in% colnames(exDatDF)), ]

picardDF <- picardDF[picardDF$X %in% colnames(exDatDF), ]

exDatMmDF <- exDatMmDF[ ,colnames(exDatMmDF) %in% metDatMmDF$HTseq_IDs]


## Add exonic reads from Picard to nMapMouseDF
# This is probably from combined mouse and human reference mapping
nMapMouseDF <- merge(nMapMouseDF, picardDF[ ,c("X", "CODING_BASES", "UTR_BASES")]
      , by.x = "SampleID", by.y = "X")
# Sum coding bases and utr bases, then divide by R2 read length (75bp)
nMapMouseDF$Exonic_Reads <- (nMapMouseDF$CODING_BASES + nMapMouseDF$UTR_BASES) / 75


## Format Visual QC

# Add HTseq_SampleID to visual QC data
visualQCdF <- merge(visualQCdF, metDatDF, by.x = c("Row", "Column")
                    , by.y = c("Fluidigm_Row", "Fluidigm_Column"))

# Format
visualQCdF$Dead.or.Alive...[is.na(visualQCdF$Dead.or.Alive...)] <- "0"
visualQCdF$Dead.or.Alive... <- gsub(".*,.*", "Multiple", visualQCdF$Dead.or.Alive...)
visualQCdF$Dead.or.Alive... <- gsub(".*0.*", "Empty", visualQCdF$Dead.or.Alive...)
visualQCdF$Dead.or.Alive... <- gsub(".*a.*", "Alive", visualQCdF$Dead.or.Alive...)
visualQCdF$Dead.or.Alive... <- gsub(".*x.*", "Dead", visualQCdF$Dead.or.Alive...)
visualQCdF$Dead.or.Alive... <- gsub("\\?", "Stained for Both", visualQCdF$Dead.or.Alive...)
################################################################################

## Calculate Number of Genes per sample with CPM > 1

# CPM (normalized by number of reads mapped to exons)
cpmDF <- t(t(exDatDF) / (metDatDF$Uniquely_Mapped / 10^6))

cpmMmDF <- t(t(exDatMmDF) / (metDatMmDF$Uniquely_Mapped / 10^6))

# Number of Genes expressed > 1 CPM
gExprd <- apply(cpmDF, 2, function(transcripts) length(subset(transcripts
                                                            , transcripts > 1)))
mean(gExprd)

gExprdMm <- apply(cpmMmDF, 2, function(transcripts) length(subset(transcripts
                                                            , transcripts > 1)))
mean(gExprdMm)
################################################################################

### Plot Number of Genes expressed per capture site after filtering for
### number of reads mapped

## Human

# Dataframe for ggplot2 of Number of Genes expressed > 1 CPM after different
# read mapped filters
df <- rbind(data.frame(Number = gExprd[nMapHumanDF$Uniquely_Mapped > 2*10^5
                                       & nMapMouseDF$Uniquely_Mapped < 10^5]
                       , Filter = "Hs_200000_Mm_50000")
  , data.frame(Number = gExprd[nMapHumanDF$Uniquely_Mapped > 10^5
                               & nMapMouseDF$Uniquely_Mapped < 10^5]
               , Filter = "Hs_100000_Mm_50000")
)
# Calculate means to add as text to boxplot
means <- aggregate(Number ~ Filter, df, mean)
means$Number <- round(means$Number, 2)

# Calculate number of samples to add as text to boxplot
nSamples <- data.frame(table(df$Filter))
colnames(nSamples) <- c("Filter", "Number_of_Samples")

# Boxplot of Number of Genes detected
ggplot(df, aes(y = Number, x = Filter)) +
  geom_boxplot() +
  geom_text(data = means, aes(label = paste("Mean:", Number), y = Number)) +
  geom_text(data = nSamples, aes(label = paste("n =", Number_of_Samples)
                                 , y = means$Number - 200)) +
  scale_x_discrete(labels = c("> 2*10^5 Hs\n< 10^5 Mm"
                              , "> 10^5 Hs\n< 10^5 Mm")) +
  xlab("Filter By Number of Reads Mapping") +
  ylab("Number of Genes > 1 CPM") +
  ggtitle("Number_Genes_Detected.R
            Number of Genes Detected
           After Filtering for Different Numbers of Reads Mapping to Human Only
           and Mouse Only")
ggsave("../analysis/graphs/Number_Genes_Detected_Boxplot.pdf")


## Plot Number of Genes expressed per capture site versus read depth

df <- data.frame(gExprd)
df <- merge(df, metDatDF[ ,c("HTseq_IDs", "Total_Reads")]
            , by.x = "row.names", by.y = "HTseq_IDs")

# Add visual QC
df <- merge(df, visualQCdF[ ,c("HTseq_IDs", "Dead.or.Alive...")]
                   , by.x = "Row.names", by.y = "HTseq_IDs", all.x = TRUE)

# Change NAs (capture sites not visually QCed) to Unkown
df$Dead.or.Alive...[is.na(df$Dead.or.Alive...)] <- "Not Visually QCed"
# Change to Factor and set "Not Visually QCed" as first factor order so that plots as
# black later
df$Dead.or.Alive... <- relevel(as.factor(df$Dead.or.Alive...), "Not Visually QCed")


df <- df[nMapHumanDF$Uniquely_Mapped > 10^5
                                 & nMapMouseDF$Uniquely_Mapped < 10^5, ]
ggplot(df, aes(x = Total_Reads, y = gExprd, color = Dead.or.Alive...)) +
  geom_point() +
  scale_color_discrete(name = "Visual QC") +
  xlab("Total Reads") +
  ylab("Number of Genes > 1 CPM") +
  ggtitle(paste("Number_Genes_Detected.R
            Number of Genes Detected Versus Read Depth
           After Filtering for > 10^5 Reads Mapping Only Human
          and < 10^5 Reads Mapping Only Mouse
          n remaining =", nrow(df)))
ggsave("../analysis/graphs/Number_Genes_Detected_Vs_Read_Depth.pdf")


## Mouse

# Dataframe for ggplot2 of Number of Genes expressed > 1 CPM after different
# read mapped filters
df <- rbind(data.frame(Number = gExprdMm[nMapHumanDF$Uniquely_Mapped < 10^5
                                        & nMapMouseDF$Uniquely_Mapped > 5*10^5]
                         , Filter = "Mm_500000_Hs_100000")
            , data.frame(Number = gExprdMm[nMapHumanDF$Uniquely_Mapped < 10^5
                                       & nMapMouseDF$Uniquely_Mapped > 2*10^5]
                       , Filter = "Mm_200000_Hs_100000")
            , data.frame(Number = gExprdMm[nMapHumanDF$Uniquely_Mapped < 10^5
                                         & nMapMouseDF$Uniquely_Mapped > 10^5]
                         , Filter = "Mm_100000_Hs_100000")
            , data.frame(Number = gExprdMm[nMapHumanDF$Uniquely_Mapped < 10^5]
                         , Filter = "Hs_100000")
)
# Calculate means to add as text to boxplot
means <- aggregate(Number ~ Filter, df, mean)
means$Number <- round(means$Number, 2)

# Calculate number of samples to add as text to boxplot
nSamples <- data.frame(table(df$Filter))
colnames(nSamples) <- c("Filter", "Number_of_Samples")

# Percent of reads total reads associated with cell barcodes remaining after
# filtering
pctTotalReadsDF <- data.frame(Mm_500000_Hs_100000 = sum(nMapMouseDF$Total_Reads[
    nMapHumanDF$Uniquely_Mapped < 10^5 & nMapMouseDF$Uniquely_Mapped > 5*10^5]) / 
    sum(nMapMouseDF$Total_Reads[nMapHumanDF$Uniquely_Mapped < 10^5]
    ) * 100
  , Mm_200000_Hs_100000 = sum(nMapMouseDF$Total_Reads[
    nMapHumanDF$Uniquely_Mapped < 10^5 & nMapMouseDF$Uniquely_Mapped > 2*10^5]) / 
    sum(nMapMouseDF$Total_Reads[nMapHumanDF$Uniquely_Mapped < 10^5]
    ) * 100
  , Mm_100000_Hs_100000 = sum(nMapMouseDF$Total_Reads[
    nMapHumanDF$Uniquely_Mapped < 10^5 & nMapMouseDF$Uniquely_Mapped > 10^5]) / 
    sum(nMapMouseDF$Total_Reads[nMapHumanDF$Uniquely_Mapped < 10^5]
    ) * 100
  , Hs_100000 = sum(nMapMouseDF$Total_Reads[
    nMapHumanDF$Uniquely_Mapped < 10^5]) / 
    sum(nMapMouseDF$Total_Reads[nMapHumanDF$Uniquely_Mapped < 10^5]
    ) * 100)
pctTotalReadsDF <- melt(pctTotalReadsDF)
pctTotalReadsDF$value <- round(pctTotalReadsDF$value, 2)
colnames(pctTotalReadsDF) <- c("Filter", "Pct_Total_Reads")

# Percent of exonic reads associated with cell barcodes remaining after
# filtering
pctExonicDF <- data.frame(Mm_500000_Hs_100000 = sum(nMapMouseDF$Exonic_Reads[
  nMapHumanDF$Uniquely_Mapped < 10^5 & nMapMouseDF$Uniquely_Mapped > 5*10^5]) / 
    sum(nMapMouseDF$Exonic_Reads[nMapHumanDF$Uniquely_Mapped < 10^5]
    ) * 100
  , Mm_200000_Hs_100000 = sum(nMapMouseDF$Exonic_Reads[
    nMapHumanDF$Uniquely_Mapped < 10^5 & nMapMouseDF$Uniquely_Mapped > 2*10^5]) / 
    sum(nMapMouseDF$Exonic_Reads[nMapHumanDF$Uniquely_Mapped < 10^5]
    ) * 100
  , Mm_100000_Hs_100000 = sum(nMapMouseDF$Exonic_Reads[
    nMapHumanDF$Uniquely_Mapped < 10^5 & nMapMouseDF$Uniquely_Mapped > 10^5]) / 
    sum(nMapMouseDF$Exonic_Reads[nMapHumanDF$Uniquely_Mapped < 10^5]
    ) * 100
  , Hs_100000 = sum(nMapMouseDF$Exonic_Reads[
    nMapHumanDF$Uniquely_Mapped < 10^5]) / 
    sum(nMapMouseDF$Exonic_Reads[nMapHumanDF$Uniquely_Mapped < 10^5]
    ) * 100)
pctExonicDF <- melt(pctExonicDF)
pctExonicDF$value <- round(pctExonicDF$value, 2)
colnames(pctExonicDF) <- c("Filter", "Pct_Exonic_Reads")

# Write out percent of reads remaining after different filters to table
outTableDF <- cbind(pctTotalReadsDF
                    , Pct_Exonic_Reads = pctExonicDF$Pct_Exonic_Reads
                    , Number_Of_Cells_After_Filter = nSamples$Number_of_Samples
                    , Mean_Number_Genes_Expressed = means$Number)
write.table(outTableDF, "../analysis/tables/Number_Genes_Detected_Vs_Percent_Of_Reads_Mouse.txt"
            , sep = "\t", quote = FALSE)

# Boxplot of Number of Genes detected
ggplot(df, aes(y = Number, x = Filter)) +
  geom_boxplot() +
  geom_text(data = means, aes(label = paste("Mean:", Number), y = Number)) +
  geom_text(data = nSamples, aes(label = paste("n =", Number_of_Samples)
                                 , y = means$Number - 200)) +
  geom_text(data = pctTotalReadsDF, aes(label = paste("Percent of total reads =", Pct_Total_Reads)
                                        , y = means$Number - 400)) +
  geom_text(data = pctExonicDF, aes(label = paste("Percent of exonic reads =", Pct_Exonic_Reads)
                                        , y = means$Number - 600)) +
  scale_x_discrete(labels = c("> 5*10^5 Mm\n< 10^5 Hs"
                              , "> 2*10^5 Mm\n< 10^5 Hs"
                              , "> 10^5 Mm\n< 10^5 Hs"
                              , "< 10^5 Hs")) +
  xlab("Filter By Number of Reads Mapping") +
  ylab("Number of Genes > 1 CPM") +
  ggtitle(paste0("Number_Genes_Detected.R"
                 , "\nNumber of genes detected for mouse"
                 , "\nFiltered for different numbers of reads uniquely mapping"
                 , "\nto human only and mouse only"
                 , "\n"))
ggsave("../analysis/graphs/Number_Genes_Detected_Mouse_Boxplot.pdf")


## Plot Number of Genes expressed per capture site versus read depth

df <- data.frame(gExprdMm)
df <- merge(df, metDatDF[ ,c("HTseq_IDs", "Total_Reads")]
            , by.x = "row.names", by.y = "HTseq_IDs")

# Add visual QC
df <- merge(df, visualQCdF[ ,c("HTseq_IDs", "Dead.or.Alive...")]
            , by.x = "Row.names", by.y = "HTseq_IDs", all.x = TRUE)

# Change NAs (capture sites not visually QCed) to Unkown
df$Dead.or.Alive...[is.na(df$Dead.or.Alive...)] <- "Not Visually QCed"
# Change to Factor and set "Not Visually QCed" as first factor order so that plots as
# black later
df$Dead.or.Alive... <- relevel(as.factor(df$Dead.or.Alive...), "Not Visually QCed")


df <- df[nMapHumanDF$Uniquely_Mapped < 10^5
         & nMapMouseDF$Uniquely_Mapped > 10^5, ]
ggplot(df, aes(x = Total_Reads, y = gExprdMm, color = Dead.or.Alive...)) +
  geom_point() +
  scale_color_discrete(name = "Visual QC") +
  xlab("Total Reads") +
  ylab("Number of Genes > 1 CPM") +
  ggtitle(paste("Number_Genes_Detected.R
          Number of Genes Detected Versus Read Depth - Mouse
          After Filtering for < 10^5 Reads Mapping Only Human
          and > 10^5 Reads Mapping Only Mouse
          n remaining =", nrow(df)))
ggsave("../analysis/graphs/Number_Genes_Detected_Vs_Read_Depth_Mouse.pdf")
