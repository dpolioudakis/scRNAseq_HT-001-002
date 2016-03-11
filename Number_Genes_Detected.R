# Damon Polioudakis
# 2016-02-29
# Graph of number of transcripts expressed
rm(list=ls())
sessionInfo()
library(ggplot2)
library(xlsx)

# Input transcript expression table (counts from HTseq) - Human
exDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon.csv", row.names = 1)

# Input transcript expression table (counts from HTseq) - Mouse
exDatMmDF <- read.csv("../data/htseq/mouse/merged/Exprs_HTSCexon.csv", row.names = 1)

# Metadata
metDatDF <- read.table("../metadata/Compiled_Metadata_20160229.txt"
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
################################################################################

# Format

# Move ERCCs to separate data frame

erccDF <- tail(exDatDF, 97)
exDatDF <- head(exDatDF, -97)

erccMmDF <- tail(exDatMmDF, 97)
exDatMmDF <- head(exDatMmDF, -97)


# Filter for missing samples

nMapMouseDF <- nMapMouseDF[(nMapMouseDF$SampleID %in% colnames(exDatDF)), ]
nMapHumanDF <- nMapHumanDF[(nMapHumanDF$SampleID %in% colnames(exDatDF)), ]

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

## Calculate number of transcripts per sample with TPM > 1

# TPM (normalized by number of reads mapped to exons)
nMapToExon <- apply(exDatDF, 2, sum)
nMapToExon <- nMapToExon / 10^6
tpmDF <- t(t(exDatDF) / nMapToExon)

nMapToExon <- apply(exDatMmDF, 2, sum)
nMapToExon <- nMapToExon / 10^6
tpmMmDF <- t(t(exDatMmDF) / nMapToExon)

# Number of transcripts expressed > 1 TPM
gExprd <- apply(tpmDF, 2, function(transcripts) length(subset(transcripts
                                                            , transcripts > 1)))
mean(gExprd)

gExprdMm <- apply(tpmMmDF, 2, function(transcripts) length(subset(transcripts
                                                            , transcripts > 1)))
mean(gExprdMm)
################################################################################

### Plot number of transcripts expressed per capture site after filtering for
### number of reads mapped

## Human

# Dataframe for ggplot2 of number of transcripts expressed > 1 TPM after different
# read mapped filters
df <- rbind(data.frame(Number = gExprd[nMapHumanDF$Uniquely_Mapped > 2*10^5
                                       & nMapMouseDF$Uniquely_Mapped < 5*10^4]
                       , Filter = "Hs_200000_Mm_50000")
  , data.frame(Number = gExprd[nMapHumanDF$Uniquely_Mapped > 10^5
                               & nMapMouseDF$Uniquely_Mapped < 5*10^4]
               , Filter = "Hs_100000_Mm_50000")
  , data.frame(Number = gExprd[nMapHumanDF$Uniquely_Mapped > 5*10^4
                               & nMapMouseDF$Uniquely_Mapped < 5*10^4]
               , Filter = "Hs_50000_Mm_50000")
)
# Calculate means to add as text to boxplot
means <- aggregate(Number ~ Filter, df, mean)
means$Number <- round(means$Number, 2)

# Calculate number of samples to add as text to boxplot
nSamples <- data.frame(table(df$Filter))
colnames(nSamples) <- c("Filter", "Number_of_Samples")

# Boxplot of number of transcripts detected
ggplot(df, aes(y = Number, x = Filter)) +
  geom_boxplot() +
  geom_text(data = means, aes(label = paste("Mean:", Number), y = Number)) +
  geom_text(data = nSamples, aes(label = paste("n =", Number_of_Samples)
                                 , y = means$Number - 100)) +
  scale_x_discrete(labels = c("> 2*10^5 Hs\n> 5*10^4 Mm"
                              , "> 10^5 Hs\n> 5*10^4 Mm"
                              , "> 5*10^4 Hs\n> 5*10^4 Mm")) +
  xlab("Filter By Number of Reads Mapping") +
  ylab("Number of Transcripts > 1 TPM") +
  ggtitle("Number_Transcripts_Detected.R
            Number of Transcripts Detected
           After Filtering for Different Numbers of Reads Mapping to Human Only
           and Mouse Only") +
  theme_grey(base_size = 12)
ggsave("../analysis/graphs/Number_Transcripts_Detected_Boxplot.pdf")


## Plot number of transcripts expressed per capture site versus read depth

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


df <- df[nMapHumanDF$Uniquely_Mapped > 5*10^4
                                 & nMapMouseDF$Uniquely_Mapped < 5*10^4, ]
ggplot(df, aes(x = Total_Reads, y = gExprd, color = Dead.or.Alive...)) +
  geom_point() +
  xlab("Total Reads") +
  ylab("Number of transcripts > 1 TPM") +
  ggtitle("Number_Transcripts_Detected.R
            Number of Transcripts Detected Versus Read Depth
           After Filtering for > 5*10^4 Reads Mapping Only Human
          and < 5*10^4 Reads Mapping Only Mouse
          n remaining =", nrow(df))) +
  theme_grey(base_size = 12)
ggsave("../analysis/graphs/Number_Transcripts_Detected_Vs_Read_Depth.pdf")


## Mouse

# Dataframe for ggplot2 of number of transcripts expressed > 1 TPM after different
# read mapped filters
df <- rbind(data.frame(Number = gExprdMm[nMapHumanDF$Uniquely_Mapped < 5*10^4
                                        & nMapMouseDF$Uniquely_Mapped > 5*10^5]
                         , Filter = "Mm_300000_Hs_50000")
            , data.frame(Number = gExprdMm[nMapHumanDF$Uniquely_Mapped < 5*10^4
                                       & nMapMouseDF$Uniquely_Mapped > 2*10^5]
                       , Filter = "Mm_200000_Hs_50000")
            , data.frame(Number = gExprdMm[nMapHumanDF$Uniquely_Mapped < 5*10^4
                                         & nMapMouseDF$Uniquely_Mapped > 10^5]
                         , Filter = "Mm_100000_Hs_50000")
            , data.frame(Number = gExprdMm[nMapHumanDF$Uniquely_Mapped < 5*10^4
                                         & nMapMouseDF$Uniquely_Mapped > 5*10^4]
                         , Filter = "Mm_50000_Hs_50000")
)
# Calculate means to add as text to boxplot
means <- aggregate(Number ~ Filter, df, mean)
means$Number <- round(means$Number, 2)

# Calculate number of samples to add as text to boxplot
nSamples <- data.frame(table(df$Filter))
colnames(nSamples) <- c("Filter", "Number_of_Samples")

# Boxplot of number of transcripts detected
ggplot(df, aes(y = Number, x = Filter)) +
  geom_boxplot() +
  geom_text(data = means, aes(label = paste("Mean:", Number), y = Number)) +
  geom_text(data = nSamples, aes(label = paste("n =", Number_of_Samples)
                                 , y = means$Number - 100)) +
  scale_x_discrete(labels = c("> 5*10^5 Mm\n> 5*10^4 Hs"
                              , "> 2*10^5 Mm\n> 5*10^4 Hs"
                              , "> 10^5 Mm\n> 5*10^4 Hs"
                              , "> 5*10^4 Mm\n> 5*10^4 Hs")) +
  xlab("Filter By Number of Reads Mapping") +
  ylab("Number of Transcripts > 1 TPM") +
  ggtitle("Number_Transcripts_Detected.R
          Number of Transcripts Detected for Mouse
          After Filtering for Different Numbers of Reads Mapping to Human Only
          and Mouse Only") +
  theme_grey(base_size = 12)
ggsave("../analysis/graphs/Number_Transcripts_Detected_Mouse_Boxplot.pdf")


## Plot number of transcripts expressed per capture site versus read depth

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


df <- df[nMapHumanDF$Uniquely_Mapped < 5*10^4
         & nMapMouseDF$Uniquely_Mapped > 5*10^4, ]
ggplot(df, aes(x = Total_Reads, y = gExprd, color = Dead.or.Alive...)) +
  geom_point() +
  xlab("Total Reads") +
  ylab("Number of transcripts > 1 TPM") +
  ggtitle(paste("Number_Transcripts_Detected.R
          Number of Transcripts Detected Versus Read Depth - Mouse
          After Filtering for < 5*10^4 Reads Mapping Only Human
          and > 5*10^4 Reads Mapping Only Mouse
          n remaining =", nrow(df))) +
  theme_grey(base_size = 12)
ggsave("../analysis/graphs/Number_Transcripts_Detected_Vs_Read_Depth_Mouse.pdf")
