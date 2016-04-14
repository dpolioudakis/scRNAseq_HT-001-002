# Damon Polioudakis
# 2016-04-13
# Make list of human only capture sites
################################################################################

rm(list=ls())
sessionInfo()

library(ggplot2)
library(xlsx)

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

exDatMmDF <- exDatMmDF[ ,colnames(exDatMmDF) %in% metDatMmDF$HTseq_IDs]


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

### Determine human only capture sites

## Select human only capture sites
names <- nMapHumanDF$SampleID[nMapHumanDF$Uniquely_Mapped > 10^5
                                       & nMapMouseDF$Uniquely_Mapped < 10^5]
length(names)
# 23

## Filter using visual QC data
# 1 capture site is multiple, removed
visualQCdF[match(names, visualQCdF$HTseq_IDs), ]
names <- names[! names == "HT_ROW20_N719"]
length(names)
# 22

## Write out as list
write.table(names
      , file = "../analysis/tables/Human_Only_Capture_Sites_10^5Hs_10^5Mm.txt"
      , quote = FALSE, row.names = FALSE, col.names = FALSE)
