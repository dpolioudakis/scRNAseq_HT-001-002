
# Damon Polioudakis
# 2016-02-29
# Add some RNA STAR statistics to sample metadata
################################################################################

rm(list=ls())

library(xlsx)

# Load data and assign variables

# RNA STAR Stats for Run 1 
stStats1dF <- read.table("../metadata/RNAstar_Stats_Mouse_SxaQSEQsXap096L1.txt"
                         , sep = "\t", header = TRUE)
# RNA STAR Stats for Run 2
stStats2dF <- read.table("../metadata/RNAstar_Stats_Mouse_SxaQSEQsXap096L2.txt"
                         , sep = "\t", header = TRUE)

# Table of CellIDs and Fluidigm Columns and Rows
metDatDF <- read.csv("../metadata/CellID_FluidigmColRow.csv", header = TRUE)
################################################################################

### Format data and add to metadata

## Remove 1 missing cell in Lane 2 from Lane 1

stStats1dF <- stStats1dF[stStats1dF$SampleID %in% stStats2dF$SampleID, ]
stStats1dF <- stStats1dF[stStats1dF$SampleID %in% metDatDF$HTseq_IDs, ]
stStats2dF <- stStats2dF[stStats2dF$SampleID %in% metDatDF$HTseq_IDs, ]

## Remove "ROW" from Fluidigm_Row (ie ROW02 -> 02)

metDatDF$Fluidigm_Row <- as.numeric(gsub("ROW", "", metDatDF$Fluidigm_Row))

##
# Calculate and add total reads, percent uniquely mapped reads, percent mapped
# to multiple loci, and percent unmapped reads to metadata table

# Total Reads
totReads <- (stStats1dF$Number.of.input.reads
             + stStats2dF$Number.of.input.reads)
metDatDF$Total_Reads <- totReads

# Uniquely Mapped - Convert % to numeric and multiply by number of reads in lane
umReads <- (stStats1dF$Uniquely.mapped.reads.number + stStats2dF$Uniquely.mapped.reads.number)
metDatDF$Uniquely_Mapped <- umReads
# Percent Uniquely Mapped
pctUMRds <- umReads / totReads
metDatDF$Pct_Uniquely_Mapped <- pctUMRds

# Total mapped to multiple loci
totMpMtLoci <- (stStats1dF$Number.of.reads.mapped.to.multiple.loci
                + stStats2dF$Number.of.reads.mapped.to.multiple.loci)
# Percent mapped to multiple loci
pctMpMtLoci <- totMpMtLoci / totReads
metDatDF$Pct_Mapped_Multiple_Loci <- pctMpMtLoci

# Percent unmapped reads
Calc_Unmapped <- function (data) {
  # Convert % to numeric and multiply by number of reads in lane
  df <- apply(data[ ,c("X..of.reads.unmapped..too.many.mismatches"
                       , "X..of.reads.unmapped..too.short"
                       , "X..of.reads.unmapped..other")]
              , 2, function(column) {as.numeric(sub("%", "", column)) / 100})
  # Multiple by number of reads in lane to convert to read number
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

##
# Write out to tab separated table
head(metDatDF)
write.table(metDatDF, file = paste0("../metadata/Compiled_Metadata_Mouse_"
                                    , format(Sys.time(), "%Y%m%d"), ".txt")
            , quote = FALSE, row.names = FALSE, sep = "\t")
