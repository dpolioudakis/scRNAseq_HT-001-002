# Damon Polioudakis
# 2016-04-06
# Select genes expressed in human cells above VZ
################################################################################

rm(list=ls())
sessionInfo()

library(ggplot2)
library(xlsx)

# Input transcript expression table (counts from HTseq) - Human
exDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon.csv", row.names = 1)

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

# Directory to output genes lists
outDir <- "../analysis/tables/expressed_genes_10tpm"
dir.create(outDir)
################################################################################

### Format

## Move ERCCs to separate data frame

erccDF <- tail(exDatDF, 97)
exDatDF <- head(exDatDF, -97)


## Filter for missing samples

nMapMouseDF <- nMapMouseDF[(nMapMouseDF$SampleID %in% colnames(exDatDF)), ]
nMapHumanDF <- nMapHumanDF[(nMapHumanDF$SampleID %in% colnames(exDatDF)), ]


## Order metadata by expression data table columns

metDatDF <- metDatDF[match(colnames(exDatDF), metDatDF$HTseq_IDs), ]

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

### TPM

# TPM (normalized by number of reads mapped)
tpmDF <- t(t(exDatDF) / (metDatDF$Uniquely_Mapped / 10^6))
################################################################################

### Make lists of expressed genes in human cells from VZ

# Select cells with >50,000 reads mapped to human only,
# < 50,000 reads mapped to mouse only, and
# from VZ brain region
hcTpmDF <- tpmDF[, nMapHumanDF$Uniquely_Mapped > 5*10^4
      & nMapMouseDF$Uniquely_Mapped < 5*10^4 & metDatDF$Fluidigm_Column > 10]

# Filter for genes > 10 TPM
genesLL <- apply(hcTpmDF, 2, function(x) x[x > 10])

## Write out gene lists for each cell
for (i in 1:length(genesLL)) {
  write.table(names(genesLL[[i]])
              , file = paste0(outDir, "/", names(genesLL)[i])
              , quote = FALSE, row.names = FALSE, col.names = FALSE)
}
# Check number of genes above TPM cutoff for each cell
lapply(genesLL, function(genes) length(names(genes)))
