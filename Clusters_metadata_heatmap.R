# Damon Polioudakis
# 2016-01-05
# Heatmap of correlation and clusters by gene expression 2015-12 scRNAseq
# run SxaQSEQsXbp060L2
rm(list=ls())
sessionInfo()

library(xlsx)
library(gplots)
################################################################################

# Inputs

# Load gene expression table (counts from HTseq)
exDat1DF <- read.csv("../data/Exprs_HTSCexon_Run1.csv"
                     , header = TRUE, row.names = 1)
exDat2DF <- read.csv("../data/Exprs_HTSCexon_Run2.csv"
                     , header = TRUE, row.names = 1)
dim(exDat1DF)
dim(exDat2DF)
# Jason Stein's FPM gene expression table from Run 1 after filtering based on
# ERCCs, I'll use this to select for that subset of cells only
fpm1DF <- read.csv("../data/FPM.csv", header = TRUE)
# Metadata
metDatDF <- read.xlsx("../metadata/PercentageReadsMapping.xlsx", 1)
pctDupDF <- read.xlsx("../metadata/PercentDuplicates.xlsx", 1)
libBatch1DF <- read.xlsx("../metadata/LibraryBatch_Well_ID.xlsx"
                         , 1, header = FALSE)
libBatch2DF <- read.xlsx("../metadata/LibraryBatch_Well_ID.xlsx"
                         , 2, header = FALSE)
totReadsRun2dF <- read.table("../metadata/Total_Reads.txt", header = FALSE)
colnames(totReadsRun2dF) <- c("CellID", "NumReads2")
################################################################################

# Formatting and subsetting data

# Combine Run 1 and Run 2 counts
exDatDF <- Reduce('+', list(exDat2DF, exDat1DF))
# Remove ERCCs
exDatDF <- exDatDF[-c((nrow(exDatDF)-96):nrow(exDatDF)), ]
# (log2 + 1 counts)
exLg2p1dF <- log(exDatDF + 1, 2)

# Edit CellIDs to be same format
colnames(exDatDF) <- gsub("_.*$", "", colnames(exDatDF))
pctDupDF <- data.frame(PERCENT_DUPLICATION = t(pctDupDF)[-1, 1])
row.names(pctDupDF) <- gsub("\\..*\\.txt", "", row.names(pctDupDF))
# Add batch ID to metadata
batc2ID <- gsub("^ID:2-", "", as.vector(t(libBatch2DF)))
batc2ID <- gsub("^ID:1-", "1_", batc2ID)
batc2ID <- gsub("_[[:digit:]].*$", "", batc2ID)
batc1ID <- gsub("^ID:2-", "", as.vector(t(libBatch1DF)))
batc1ID <- gsub("^ID:1-", "1_", batc1ID)
batc1ID <- gsub("_[[:digit:]].*$", "", batc1ID)
batchIDdF <- rbind(data.frame(LibraryBatch = "1", CaptureWellID = batc1ID)
                   , data.frame(LibraryBatch = "2", CaptureWellID = batc2ID))
metDatDF <- merge(x = metDatDF, y = batchIDdF,
                  by.x = "CaptureWellID", by.y = "CaptureWellID")
# Add read depth from Run 2 to metadata
metDatDF <- merge(metDatDF, totReadsRun2dF, by.x = "CellID", by.y = "CellID")
metDatDF$TotalReads <- metDatDF$NumReads + metDatDF$NumReads2

# Subset cells to those in Jason's FPM gene expression table that was filtered
# based on ERCCs
fpmCells <- gsub("X", "Cell", colnames(fpm1DF))
exLg2p1FtDF <- exLg2p1dF[ , colnames(exLg2p1dF) %in% fpmCells]
metFtDF <- metDatDF[metDatDF$CellID %in% fpmCells, ]
# colnames() <- gsub("_.*", "", colnames(exDat1DF))
# ex1fTdF <- exDat1DF[ , colnames(exDat1DF) %in% fpmCells]
# ex2fTdF <- exDat2DF[ , colnames(exDat2DF) %in% fpmCells]

# Filter for only genes expressed >=1 count in >= 10% of cells
pres <- apply(exLg2p1FtDF >= 1, 1, sum) 
idx <- pres >= 0.1 * dim(exLg2p1FtDF)[2] ## exp > 10 in 80% of samples
exLg2p1FtErExDF = exLg2p1FtDF[idx,]


# Filter for only genes expressed >=1 count in >= 10% of cells
pres <- apply(exDatDF >= 1, 1, sum) 
idx <- pres >= 0.1 * dim(exDatDF)[2] ## exp > 10 in 80% of samples
exFtExDF = exDatDF[idx,]
# Normalize for read depth
tpmDF <- apply(exFtExDF, 1, function(gene) (gene/metDatDF$TotalReads)*10^6)
tpmDF <- data.frame(t(tpmDF))
# log(tpm + 1)
tpmLg2p1dF <- log(tpmDF + 1, 2)
# Subset cells to those in Jason's FPM gene expression table that was filtered
# based on ERCCs
fpmCells <- gsub("X", "Cell", colnames(fpm1DF))
tpmLg2p1fTdF <- tpmLg2p1dF[ , colnames(tpmLg2p1dF) %in% fpmCells]
metFtDF <- metDatDF[metDatDF$CellID %in% fpmCells, ]



################################################################################

MakeCorHeatmap <- function (metDat, metDatDF, corM, fileNameSfx
                            , corMethod, ftERCC) {
  # Set same sample correlation to NA to be able to see other heatmap colors
  diag(corM) <- NA
  jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan"
                                  , "#7FFF7F", "yellow", "#FF7F00", "red"
                                  , "#7F0000"))
  colsidecolorbar <- c("black", "blue", "green", "red", "grey", "yellow", "orange"
                       , "grey60", "grey20")
  colsidecolors <- colsidecolorbar[match(metDatDF[[metDat]]
                                         , unique(metDatDF[[metDat]]))]
  
  pdf(paste("../analysis/graphs/Heatmap_Clusters_", fileNameSfx, metDat
            , ".pdf", sep = ""), width = 16, height = 16)
  # par(mai = c(2,2,2,6))
  heatmap.2(corM, trace = "none", col = jet.colors(40)
          , ColSideColors = colsidecolors
          , RowSideColors = colsidecolors, cexRow = 0.2, cexCol = 0.2
          # , margins = c(16, 4)
          , xlab = "Sample ID", ylab = "Sample ID"
          , main = paste("Clustering Based on Inter-Sample ", corMethod
                         , " Coefficient"
                         , "\nlog2(counts + 1) ", ftERCC
                         , "\nColor bar indicates: ", metDat, sep = "")
          , keysize = 1)
  dev.off()
}

# Make a correlation matrix of all samples
corM <- cor(exLg2p1dF, method = "pearson")
MakeCorHeatmap("CapturePlateID", metDatDF, corM, "Pearson_", "Pearson"
               , "All Cells")
MakeCorHeatmap("VisualQC", metDatDF, corM, "Pearson_", "Pearson"
               , "All Cells")
MakeCorHeatmap("LibraryBatch", metDatDF, corM, "Pearson_", "Pearson"
               , "All Cells")

corM <- cor(exLg2p1FtDF, method = "pearson")
MakeCorHeatmap("CapturePlateID", metFtDF, corM, "Pearson_FtERCC_", "Pearson"
               , "Cells filtered based on ERCCs")
MakeCorHeatmap("VisualQC", metFtDF, corM, "Pearson_FtERCC_", "Pearson"
               , "Cells filtered based on ERCCs")
MakeCorHeatmap("LibraryBatch", metFtDF, corM, "Pearson_FtERCC_", "Pearson"
               , "Cells filtered based on ERCCs")

# log2(counts + 1) filtered for counts >= 1 in >= 10% of cells and cells filtered
# for those remaining after Jason's ERCC filtering
corM <- cor(exLg2p1FtErExDF, method = "pearson")
MakeCorHeatmap("CapturePlateID", metFtDF, corM, "Pearson_FtERCCg1c10c_", "Pearson"
               , "Cells filtered based on ERCCs")
MakeCorHeatmap("VisualQC", metFtDF, corM, "Pearson_FtERCCg1c10c_", "Pearson"
               , "Cells filtered based on ERCCs")
MakeCorHeatmap("LibraryBatch", metFtDF, corM, "Pearson_FtERCCg1c10c_", "Pearson"
               , "Cells filtered based on ERCCs")

# log2(TPM + 1) filtered for counts >= 1 in >= 10% of cells and cells filtered
# for those remaining after Jason's ERCC filtering
corM <- cor(tpmLg2p1fTdF, method = "pearson")
MakeCorHeatmap("CapturePlateID", metFtDF, corM, "Pearson_TPM_FtERCCg1c10c_", "Pearson"
               , "Cells filtered based on ERCCs")
MakeCorHeatmap("VisualQC", metFtDF, corM, "Pearson_TPM_FtERCCg1c10c_", "Pearson"
               , "Cells filtered based on ERCCs")
MakeCorHeatmap("LibraryBatch", metFtDF, corM, "Pearson_TPM_FtERCCg1c10c_", "Pearson"
               , "Cells filtered based on ERCCs")

corM <- cor(exLg2p1dF, method = "spearman")
MakeCorHeatmap("CapturePlateID", metDatDF, corM, "Spearman_", "Spearman"
               , "All Cells")
MakeCorHeatmap("VisualQC", metDatDF, corM, "Spearman_", "Spearman"
               , "All Cells")
MakeCorHeatmap("LibraryBatch", metDatDF, corM, "Spearman_", "Spearman"
               , "All Cells")

corM <- cor(exLg2p1FtDF, method = "spearman")
MakeCorHeatmap("CapturePlateID", metFtDF, corM, "Spearman_FtERCC_", "Spearman"
               , "Cells filtered based on ERCCs")
MakeCorHeatmap("VisualQC", metFtDF, corM, "Spearman_FtERCC_", "Spearman"
               , "Cells filtered based on ERCCs")
MakeCorHeatmap("LibraryBatch", metFtDF, corM, "Spearman_FtERCC_", "Spearman"
               , "Cells filtered based on ERCCs")
################################################################################


#Method to perform hierarchical clustering and cut at a certain height

# log2(counts + 1) filtered for counts >= 1 in >= 10% of cells and cells filtered
# for those remaining after Jason's ERCC filtering
corM <- cor(exLg2p1FtErExDF, method = "pearson")
diag(corM) <- NA
hm <- heatmap.2(corM, keysize = 1, trace = "none")
hc <- as.hclust( hm$rowDendrogram )
groups <- cutree( hc, k=10 )
jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan"
                                , "#7FFF7F", "yellow", "#FF7F00", "red"
                                , "#7F0000"))
pdf(paste("../analysis/graphs/Heatmap_Cut_Clusters_log2_counts.pdf", sep = "")
    , width = 16, height = 16)
heatmap.2(corM, keysize = 1, trace = "none", col = jet.colors(40)
          , ColSideColors = as.character(groups))
dev.off()
clusters <- split(as.data.frame(t(exLg2p1FtErExDF)), as.vector(groups))
clusters <- lapply(clusters, t)

# Write each cluster to .csv
for (i in 1:length(clusters)) {
  print(clusters[i])
  write.csv(clusters[i], file = paste("../data/clusters/Cluster_log2_counts_"
                                        , i, ".csv", sep = ""))
}

# log2(TPM + 1) filtered for counts >= 1 in >= 10% of cells and cells filtered
# for those remaining after Jason's ERCC filtering
corM <- cor(tpmLg2p1fTdF, method = "pearson")
diag(corM) <- NA
hm <- heatmap.2(corM, keysize = 1, trace = "none")
hc <- as.hclust( hm$rowDendrogram )
groups <- cutree( hc, k=10 )
jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan"
                                , "#7FFF7F", "yellow", "#FF7F00", "red"
                                , "#7F0000"))
pdf(paste("../analysis/graphs/Heatmap_Cut_Clusters_TPM.pdf", sep = "")
    , width = 16, height = 16)
heatmap.2(corM, keysize = 1, trace = "none", col = jet.colors(40)
          , ColSideColors = as.character(groups))
dev.off()
clusters <- split(as.data.frame(t(exLg2p1FtErExDF)), as.vector(groups))
clusters <- lapply(clusters, t)

# Write each cluster to .csv
for (i in 1:length(clusters)) {
  print(clusters[i])
  write.csv(clusters[i], file = paste("../data/clusters/Cluster_TPM_"
                                        , i, ".csv", sep = ""))
}
