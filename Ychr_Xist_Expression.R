# Damon Polioudakis
# 2016-08-01
# Plot Y chromosome and Xist expression for Fluidigm HT-001-002 mouse wells
################################################################################

rm(list=ls())
sessionInfo()

require(ggplot2)
require(xlsx)
require(reshape2)
require(biomaRt)

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
# Reads aligned to mouse only (excluded those that also map to human)
nMapMouseDF <- read.table("../metadata/Reads_Aligned_Only_Mouse.txt"
                          , header = TRUE, sep = "\t")

# Visual QC
visualQCdF <- read.xlsx("../metadata/ExpID3-C1-HT-CaptureData-2016-02-02.xlsx", 1)

# HTseqID to SampleID
idMapDF <- read.table("../metadata/Compiled_Metadata_20160229.txt", header = TRUE)

# Picard Tools QC statistics
picardDF <- read.csv("../metadata/PicardToolsQC_Mouse.csv")

dir.create("../analysis/graphs/", recursive = TRUE)

outGraph <- "../analysis/graphs/Ychr_Xist_Expression_"
graphTitle <- "Ychr_Xist_Expression.R"

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


## Chomp . and everything after from ensembl IDs
# Mouse
row.names(exDatMmDF) <- gsub("\\..*", "", row.names(exDatMmDF))
# Human
row.names(exDatDF) <- gsub("\\..*", "", row.names(exDatDF))
################################################################################

### Filter to mouse or human capture sites

# Dataframe for ggplot2 of Number of Genes expressed > 1 CPM after different
# read mapped filters
mmFtExDF <- exDatMmDF[ ,nMapHumanDF$Uniquely_Mapped < 10^5
                          & nMapMouseDF$Uniquely_Mapped > 10^5]

hsFtExDF <- exDatDF[ ,nMapHumanDF$Uniquely_Mapped > 10^5
                       & nMapMouseDF$Uniquely_Mapped < 10^5]


## Convert Ensembl IDs to Gene Symbols

AddChromosomeMouse <- function (ensemblList) {
  ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  moduleGenes <- data.frame(ensemblList)
  # bioMart manual:
  #http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
  # Attribute: values to retrieve
  # Filters: input query type
  # Values: input query
  #ensembl <- useMart("ensembl")
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  ensembl <- useDataset("mmusculus_gene_ensembl", mart=ensembl)
  # Data frame of module Ensembl IDs and gene symbols
  ensemblGeneSymDF <- getBM(  attributes = c("ensembl_gene_id", "mgi_symbol", "chromosome_name")
                              , filters = "ensembl_gene_id"
                              , values = moduleGenes
                              , mart = ensembl
  )
  ensemblGeneSymDF
}

AddChromosomeHuman <- function (ensemblList) {
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  moduleGenes <- data.frame(ensemblList)
  # bioMart manual:
  #http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
  # Attribute: values to retrieve
  # Filters: input query type
  # Values: input query
  #ensembl <- useMart("ensembl")
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
  # Data frame of module Ensembl IDs and gene symbols
  ensemblGeneSymDF <- getBM(  attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name")
                              , filters = "ensembl_gene_id"
                              , values = moduleGenes
                              , mart = ensembl
  )
  ensemblGeneSymDF
}

# Dataframe of ensembl IDs and chromosome
# Mouse
mmEnsemblChrDF <- AddChromosomeMouse(row.names(mmFtExDF))
# Human
hsEnsemblChrDF <- AddChromosomeHuman(row.names(hsFtExDF))

# Add chromosome to expression dataframe
# Mouse
mmFtExDF <- merge(mmEnsemblChrDF, mmFtExDF, by.x = "ensembl_gene_id", by.y = "row.names")
# Human
hsFtExDF <- merge(hsEnsemblChrDF, hsFtExDF, by.x = "ensembl_gene_id", by.y = "row.names")

## Y expression
pdf(paste0(outGraph, "Mean_Ychr_Expression.pdf"))
# Mouse
yMmFtExDF <- mmFtExDF[mmFtExDF$chromosome_name == "Y", ]
barplot(sort(apply(yMmFtExDF[ ,-c(1:3)], 2, mean))
        , xaxt = "n", xlab = "Cells"
        , ylab = "Expression (counts)"
        , main = paste0(graphTitle
                        , "\n", "Mouse: Fluidigm HT mean Y chromosome expression"
                        , "\n", "> 10^5 mouse reads aligned, < 10^5 human reads aligned"))
# Human
yHsFtExDF <- hsFtExDF[hsFtExDF$chromosome_name == "Y", ]
barplot(sort(apply(yHsFtExDF[ ,-c(1:3)], 2, mean))
        , xaxt = "n", xlab = "Cells"
        , ylab = "Expression (counts)"
        , main = paste0(graphTitle
                        , "\n", "Human: Fluidigm HT mean Y chromosome expression"
                        , "\n", "> 10^5 human reads aligned, < 10^5 mouse reads aligned"))
dev.off()

## Xist expression
pdf(paste0(outGraph, "Xist_Expression.pdf"))
# Mouse
xistMmM <- as.matrix(mmFtExDF[mmFtExDF$mgi_symbol == "Xist", -c(1:3)])
barplot(xistMmM
        , xaxt = "n", xlab = "Cells"
        , ylab = "Expression (counts)"
        , main = paste0(graphTitle
                        , "\n", "Mouse: Fluidigm HT Xist expression"
                        , "\n", "> 10^5 mouse reads aligned, < 10^5 human reads aligned"))
# Human
xistHsM <- as.matrix(hsFtExDF[hsFtExDF$hgnc_symbol == "XIST", -c(1:3)])
barplot(xistHsM
        , xaxt = "n", xlab = "Cells"
        , ylab = "Expression (counts)"
        , main = paste0(graphTitle
                        , "\n", "Human: Fluidigm HT Xist expression"
                        , "\n", "> 10^5 human reads aligned, < 10^5 mouse reads aligned"))
dev.off()

## Log2 (expression + 1) by chromosome
# Mouse
ggDF <- melt(mmFtExDF)
ggDF <- ggDF[ggDF$chromosome_name != "MT", ]
ggDF$value <- log(ggDF$value + 1, 2)
aggregate(ggDF$value, list(chr = ggDF$chromosome_name), mean)
ggplot(ggDF, aes(y = value, x = chromosome_name)) +
  geom_jitter(alpha = 0.05) +
  ylab("Log2(counts + 1)") +
  xlab("Chromosome") + 
  ggtitle(paste0(graphTitle
                 , "\n", "Mouse: Fluidigm HT expression by chromosome"
                 , "\n", "> 10^5 mouse reads aligned, < 10^5 human reads aligned"
                 , "\n"))
ggsave(paste0(outGraph, "Expression_By_Chromosome_Mm.png"), width = 9, height = 6)
# Human
ggDF <- melt(hsFtExDF)
# Remove biomart extra chromosomes
rmvChr <- c(grep("CHR*", unique(ggDF$chromosome_name), value = TRUE)
            , grep("KI27*", unique(ggDF$chromosome_name), value = TRUE))
ggDF <- ggDF[! ggDF$chromosome_name %in% rmvChr, ]
ggDF$value <- log(ggDF$value + 1, 2)
aggregate(ggDF$value, list(chr = ggDF$chromosome_name), mean)
ggplot(ggDF, aes(y = value, x = chromosome_name)) +
  geom_jitter(alpha = 0.25) +
  ylab("Log2(counts + 1)") +
  xlab("Chromosome") + 
  ggtitle(paste0(graphTitle
         , "\n", "Human: Fluidigm HT expression by chromosome"
         , "\n", "> 10^5 human reads aligned, < 10^5 mouse reads aligned"
         , "\n"))
ggsave(paste0(outGraph, "Expression_By_Chromosome_Hs.png"), width = 9, height = 6)  
  
  



