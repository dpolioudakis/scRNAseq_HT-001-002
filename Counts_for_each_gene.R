# Damon Polioudakis
# 2016-01-07
# Graph of counts per gene  in 2015-12 scRNAseq run C196-001
# SxaQSEQsXbp060L2 and 2015-12 scRNAseq additional sequencing C196-002
rm(list=ls())
sessionInfo()

library(ggplot2)
library(reshape2)
library(xlsx)

# Input gene expression table (counts from HTseq)
exDat1DF <- read.csv("../data/Exprs_HTSCexon_Run1.csv", header = TRUE, row.names = 1)
exDat2DF <- read.csv("../data/Exprs_HTSCexon_Run2.csv", header = TRUE, row.names = 1)
fpm1DF <- read.csv("../data/FPM.csv", header = TRUE)
dim(exDat1DF)
dim(exDat2DF)
metDatDF <- read.xlsx("../metadata/PercentageReadsMapping.xlsx", 1)
totReadsRun2dF <- read.table("../metadata/Total_Reads.txt", header = FALSE)
colnames(totReadsRun2dF) <- c("CellID", "NumReads2")

# Add read depth from Run 2 to metadata
metDatDF <- merge(metDatDF, totReadsRun2dF, by.x = "CellID", by.y = "CellID")
metDatDF$TotalReads <- metDatDF$NumReads + metDatDF$NumReads2
metFtDF <- metDatDF

# Remove ERCCs
exDat1DF <- exDat1DF[-((nrow(exDat1DF)-96):nrow(exDat1DF)), ]
exDat2DF <- exDat2DF[-((nrow(exDat2DF)-96):nrow(exDat2DF)), ]

# Combine counts from scRNAseq Run1 and Run2
exDatDF <- Reduce('+', list(exDat2DF, exDat1DF))

# Subset to cells remaining after Jason's ERCC filters
fpmCells <- gsub("X", "Cell", colnames(fpm1DF))
exDatDF <- exDatDF[ , colnames(exDatDF) %in% fpmCells]
ex1fTdF <- exDat1DF[ , colnames(exDat1DF) %in% fpmCells]
ex2fTdF <- exDat2DF[ , colnames(exDat2DF) %in% fpmCells]
# Subset metadata
metFtDF <- metDatDF[metDatDF$CellID %in% fpmCells, ]

# Percentage of mapped reads made up by top X genes
nMapDF <- data.frame(apply(exDatDF, 2, sum))
topGenesDF <- data.frame(apply(exDatDF, 2, function(Cell) sum(sort(Cell, decreasing = TRUE)[1:10])))
pctTopDF <- topGenesDF/nMapDF[ ,1]
colnames(pctTopDF) <- "PercentTop"
pctTopDF$CellID <- row.names(pctTopDF)

# Number of reads mapped
nMapDF <- data.frame(apply(exDatDF, 2, sum))
colnames(nMapDF) <- "ReadsMapped"
nMapDF$CellID <- row.names(nMapDF)
# Format DF of gene expression for ggplot2
exDatFrmtDF <- exDatDF
exDatFrmtDF$Gene <- row.names(exDatDF)
exDatFrmtDF <- melt(exDatFrmtDF, id.vars = "Gene")
ggplot(exDatFrmtDF, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1) +
  geom_point(data = metFtDF, aes(x = CellID, y = TotalReads/100, color = "Total Reads / 100")) +
  geom_point(data = nMapDF, aes(x = CellID, y = ReadsMapped/100, color = "Reads Mapped to Exons / 100")) +
  geom_point(data = pctTopDF, aes(x = CellID, y = PercentTop*10000, color = "Percent Mapped to Top 10 Genes * 10^4")) +
  scale_colour_manual(values=c("red", "blue", "purple")) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.x=element_blank()) +
  xlab("Cell") +
  ylab("Counts") +
  ggsave("../analysis/graphs/Jackpotting_unfiltered.pdf")
  
  
ggplot(exDatFrmtDF, aes(x = variable, y = log(value + 1,2))) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1) +
  geom_point(data = metFtDF, aes(x = CellID, y = TotalReads/100000), color = "red") +
  geom_point(data = nMapDF, aes(x = CellID, y = ReadsMapped/100000), color = "blue") +
  geom_point(data = pctTopDF, aes(x = CellID, y = PercentTop*100), color = "purple")


# List of genes expressed >=1 count in >= 10% of cells
pres <- apply(exDatDF >= 1, 1, sum) 
idx <- pres >= 0.05 * dim(exDatDF)[2] ## exp >= 1 in 10% of samples

# Filter for genes not expressed >=1 count in >= 10% of cells
notPresDF <- exDatDF[!idx,]
# Number of reads mapped
nMapDF <- data.frame(apply(notPresDF, 2, sum))
colnames(nMapDF) <- "ReadsMapped"
nMapDF$CellID <- row.names(nMapDF)
# Format DF of gene expression for ggplot2
notPresDF$Gene <- row.names(notPresDF)
notPresDF <- melt(notPresDF, id.vars = "Gene")
ggplot(notPresDF, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1) +
  geom_point(data = metFtDF, aes(x = CellID, y = TotalReads/1000), color = "red") +
  geom_point(data = nMapDF, aes(x = CellID, y = ReadsMapped/100), color = "blue") +
  geom_point(data = pctTopDF, aes(x = CellID, y = PercentTop*10000), color = "purple")

# Filter for genes not expressed >=1 count in >= 10% of cells
exFtExDatDF <- exDatDF[idx,]
nMapDF <- data.frame(apply(notPresDF, 2, sum))
colnames(nMapDF) <- "ReadsMapped"
nMapDF$CellID <- row.names(nMapDF)
apply(exFtExDatDF, 2, mean)
boxplot(exFtExDatDF[ ,1:ncol(exFtExDatDF)])
exFtExDatDF$Gene <- row.names(exFtExDatDF)
exFtExDatDF <- melt(exFtExDatDF, id.vars = "Gene")
pdf("../analysis/graphs/tmp.pdf")
ggplot(exFtExDatDF, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1) +
  geom_point(data = metFtDF, aes(x = CellID, y = TotalReads/100), color = "red") +
  geom_point(data = nMapDF, aes(x = CellID, y = ReadsMapped/10), color = "blue") +
  geom_point(data = pctTopDF, aes(x = CellID, y = PercentTop*10000), color = "purple")
dev.off()







# Normalize for read depth
tpmDF <- apply(exDatDF, 1, function(gene) (gene/metDatDF$TotalReads)*10^6)
tpmDF <- data.frame(t(tpmDF))
lg2tpmDF <- log(tpmDF + 1, 2)

pres <- apply(lg2tpmDF >= 1, 1, sum) 
idx <- pres >= 0.1 * dim(lg2tpmDF)[2] ## exp >= 1 in 10% of samples

notPresDF <- lg2tpmDF[!idx,]
apply(notPresDF, 2, range)
apply(notPresDF, 2, mean)
boxplot(notPresDF[ ,1:ncol(notPresDF)])

tpmFtExDF <- lg2tpmDF[idx,]
apply(tpmFtExDF, 2, range)
apply(tpmFtExDF, 2, mean)
pdf(paste("../analysis/graphs/tmp_TPM.pdf", sep = ""), width = 16, height = 16)
boxplot(tpmFtExDF[ ,1:ncol(tpmFtExDF)])
dev.off()

# Check Run 1 Jason's TPMs
fpm1DF <- fpm1DF[ ,c(5:(ncol(fpm1DF)-1))]
pdf(paste("../analysis/graphs/tmp_JS.pdf", sep = ""), width = 16, height = 16)
boxplot(fpm1DF[ ,1:ncol(fpm1DF)])
dev.off()







d <- exDatDF[ ,2:5]
d$Genes <- row.names(d)
d <- melt(d, id.vars = "Genes", variable.name = "CellID", value.name = "Counts")

ggplot(d, aes(x = Genes, y = Counts)) +
  facet_wrap(~CellID) +
  geom_bar(stat = "identity")

ggplot(d, aes(x = Genes, y = Counts)) +
  geom_point()

ggplot(data.frame(exDatDF[ ,1]), aes(x = row.names(exDatDF), y = as.numeric(exDatDF[ ,1]))) +
  geom_point()

barplot(sort(exDatDF[ ,3]))
head(exDatDF)
