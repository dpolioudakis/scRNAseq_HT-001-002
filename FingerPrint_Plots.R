# Damon Polioudakis
# 2016-05-23
# Heatmap and Histogram of variants detected for human cell subset from
# Fluidigm HT mouse human VZ and CP experiment
################################################################################

rm(list = ls())
sessionInfo()

require(reshape2)
require(ggplot2)

### Load data and assign variables

## Load data
vcfDF <- read.table("../data/vcf/Compiled_FingerPrint_VCFs.txt"
                    , comment.char = "", fill = TRUE, row.names = NULL)
colnames(vcfDF) <- c(colnames(vcfDF)[-1], "GENOTYPE_INFO")

## Variables
graphCodeTitle <- "FingerPrint_Plots.R"
outGraph <- "../analysis/graphs/FingerPrint_Plots_"

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text = element_text(size = 18)))
theme_update(plot.title = element_text(size = 14))
################################################################################

### Process and Plot Data

## Process

hsCellIDs <- unique(vcfDF$CellIDs)

# Filter SNPs GATK calls as Low Quality
vcfDF <- vcfDF[! vcfDF$FILTER == "LowQual",1:8]
# Remove unused factor levels
vcfDF$X.CHROM <- factor(vcfDF$X.CHROM)
# Check SNP is called same
split(vcfDF, vcfDF$X.CHROM)

# Make 0,1 matrix of called or not called in each sample
mx <- acast(vcfDF[ ,1:2], X.CHROM~CellIDs)
mx[mx > 1] <- 1

# Add Cell IDs for Cells with no variants detected to matrix
noVarCellIDs <- hsCellIDs[! hsCellIDs %in% colnames(mx)]
noVarMx <- matrix(data = 0, nrow(mx), length(noVarCellIDs))
colnames(noVarMx) <- noVarCellIDs
mx <- cbind(mx, noVarMx)

## Plot

# Heatmap of SNPs vs Samples called or not
ggDF <- melt(mx)

ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8)) +
  ylab("Variant location") +
  xlab("Cell ID") +
  ggtitle(paste0(graphCodeTitle
                 , "\nVariants called for each Cell: Fluidigm HT, Human VZ and CP"
                 , "\n"))
ggsave(paste0(outGraph, "Heatmap.pdf"), height = 12)


# Histogram of SNPs vs Samples called or not

ggDF <- data.frame(Count = rowSums(mx))

ggplot(ggDF, aes(x = Count)) +
  geom_histogram(binwidth = 1, origin = -0.5, col = "black") +
  scale_x_continuous(breaks = 0:12) +
  ylab("Count") +
  xlab("Number of cells variant is detected in") +
  ggtitle(paste0(graphCodeTitle
                 , "\nHistogram of number of cells in which each variant is called"
                 , "\n"))
ggsave(paste0(outGraph, "Histogram.pdf"))


# Number of variants detected per cell

ggDF <- data.frame(Number = colSums(mx))
ggDF$CellIDs <- row.names(ggDF)
ggplot(ggDF, aes(y = Number, x = CellIDs)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Variants detected") +
  xlab("Cell IDs") +
  ggtitle(paste0(graphCodeTitle
                 , "\nNumber of variants detected per cell"
                 , "\nMedian variants per cell: ", median(ggDF$Number)
                 , "\n"))
ggsave(paste0(outGraph, "Variants_Per_Cell.pdf"))
  


## If mean frequency of variants is 0.25% then
0.025 * 0.975 * 2
# = 0.04875
# Then 4.875% of variants would be ref in Donor 1 and alt in Donor 2
# Mean we need ~20-30 variants detected per cell to call Donor 1 vs Donor 2
