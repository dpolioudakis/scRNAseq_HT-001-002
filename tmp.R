# Damon Polioudakis
# 2016-05-23
# Heatmap and Histogram of SNPs detected
################################################################################

rm(list = ls())
sessionInfo()

require(reshape2)
require(ggplot2)

vcfDF <- read.table("../data/vcf/Compiled_FingerPrint_VCFs.txt"
                    , comment.char = "", fill = TRUE, row.names = NULL)
colnames(vcfDF) <- c(colnames(vcfDF)[-1], "GENOTYPE_INFO")

# Filter SNPs GATK calls as Low Quality
vcfDF <- vcfDF[! vcfDF$FILTER == "LowQual",1:8]
# Remove unused factor levels
vcfDF$X.CHROM <- factor(vcfDF$X.CHROM)
# Check SNP is called same
split(vcfDF, vcfDF$X.CHROM)

# Make 0,1 matrix of called or not called in each sample
mx <- acast(vcfDF[ ,1:2], X.CHROM~CellIDs)
mx[mx > 1] <- 1


# Heatmap of SNPs vs Samples called or not
ggDF <- melt(mx)

ggplot(ggDF, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Variant Location") +
  xlab("CellID")


# Histogram of SNPs vs Samples called or not

ggDF <- data.frame(Count = rowSums(mx))

ggplot(ggDF, aes(x = Count)) +
  geom_histogram(binwidth = 1, origin = -0.5, col = ) +
  scale_x_continuous(breaks = 0:12)


0.025 * 0.975 * 2
