
# Damon Polioudakis
# 2016-02-29
# Convert Row name and Nextera Indexes to unique cells IDs

################################################################################

rm(list=ls())
sessionInfo()

library(xlsx)

# Table of Nextera Indexes matched to Fluidigm Column Numbers
idxColDF <- read.xlsx("../metadata/Nextera_Indexes.xlsx", 1)

# Load expression data - raw counts output by HTseq
exExDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon.csv", row.names = 1)
exGeDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCgene.csv", row.names = 1)
################################################################################

# Sample names from HTseq expression data
cNames <- colnames(exExDatDF)

# Loop through Sample names, extract Nextera Index, and match Index to
# Fluidigm Column in Nextera_Indexes.xlsx
# Output dataframe of Fluidigm Column Numbers, Row Numbers, and Nexter Indexes
metDatDF <- data.frame()
for (cName in cNames) {
  idxName <- gsub("HT_ROW.._", "", cName)
  fluCol <- idxColDF[idxColDF$Index_Name == idxName, ]$Fluidigm_Column
  fluRow <- gsub("HT_*", "", cName)
  fluRow <- gsub("_N...", "", fluRow)
  metDatDF <- rbind(metDatDF
                    , data.frame(Fluidigm_Column = fluCol
                               , Fluidigm_Row = fluRow
                               , Nextera_Index = idxName))
}

# Add HTseq Sample IDs
metDatDF$HTseq_IDs <- cNames

# Order by Column Number
metDatDF <- metDatDF[order(metDatDF$Fluidigm_Column), ]

# Make and add CellID
metDatDF$CellID <- seq(1:length(cNames))

# Order by Row Number to prepare to change column names of Expression Data to CellID
metDatDF <- metDatDF[order(metDatDF$Fluidigm_Row), ]

# Change Column Names of Expression Data to CellID and save
colnames(exExDatDF) <- paste0("CellID_", metDatDF$CellID)
exExDatDF <- exExDatDF[ ,order(colnames(exExDatDF))]
write.csv(exExDatDF, "../data/htseq/merged/Exprs_HTSCexon_CellIDs.csv", quote = FALSE)

colnames(exGeDatDF) <- paste0("CellID_", metDatDF$CellID)
exGeDatDF <- exGeDatDF[ ,order(colnames(exGeDatDF))]
write.csv(exGeDatDF, "../data/htseq/merged/Exprs_HTSCgene_CellIDs.csv", quote = FALSE)

# Save CellID, Fluidigm Column and Row table
write.csv(metDatDF, "../metadata/CellID_FluidigmColRow.csv", quote = FALSE, row.names = FALSE)
