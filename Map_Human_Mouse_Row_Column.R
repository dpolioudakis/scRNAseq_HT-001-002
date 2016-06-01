# Damon Polioudakis
# 2016-05-31
# Heatmap of mouse or human capture sites ordered by HT row and and column
################################################################################

rm(list=ls())
sessionInfo()

library(ggplot2)
library(xlsx)
library(reshape2)

### Load data and assign variables

## Load data

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

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 18)))
theme_update(plot.title = element_text(size = 14))
theme_update(axis.title.y = element_text(angle = 0))

## Variables
graphCodeTitle <- "Map_Human_Mouse_Row_Column.R"
outGraph <- "../analysis/graphs/Map_Human_Mouse_Row_Column_"
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

### Determine capture site species and plot according to Fluidigm chip row and
### column

## Add species information based on reads mapped
nMapHumanDF$Species <- "NA"
nMapHumanDF$Species[nMapHumanDF$Uniquely_Mapped > 10^5
                              & nMapMouseDF$Uniquely_Mapped < 10^5] <- "Human"
nMapHumanDF$Species[nMapHumanDF$Uniquely_Mapped < 10^5
                              & nMapMouseDF$Uniquely_Mapped > 10^5] <- "Mouse"
nMapHumanDF$Species[nMapHumanDF$Uniquely_Mapped < 10^5
                              & nMapMouseDF$Uniquely_Mapped < 10^5] <- "Unknown"
nMapHumanDF$Species[nMapHumanDF$Uniquely_Mapped > 10^5
                              & nMapMouseDF$Uniquely_Mapped > 10^5] <- "Both"


## Heatmap 

ggDF <- merge(metDatDF[ ,2:4], nMapHumanDF[ ,c(1,6)]
            , by.x = "HTseq_IDs", by.y = "SampleID")

ggplot(ggDF, aes(x = Fluidigm_Column, y = Fluidigm_Row, fill = Species)) +
  geom_tile() +
  ylab("Fluidigm Chip Row") +
  xlab("Fluidigm Chip Column") +
  ggtitle(paste0(graphCodeTitle
                 , "\nCapture Sites With Reads Aligned to Mouse or Human"
                 , "\nFluidigm HT, Human VZ and CP and Mouse Progenitors"
                 , "\n"))
ggsave(paste0(outGraph, "Heatmap.pdf"), height = 12)