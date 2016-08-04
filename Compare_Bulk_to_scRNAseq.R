
# Damon Polioudakis
# 2016-02-10
# Plot bulk RNAseq of VZ and CP from Luis and Jason's ATAC versus pooled
# scRNAseq VZ and CP

# Inputs
#   HTseq counts for bulk RNAseq VZ and CP from Luis and Jason's ATAC
#   HTseq counts for scRNAseq C196-001_002

# Outputs

################################################################################

rm(list=ls())
sessionInfo()

library(ggplot2)

# Load data and assign variables

buExDatDF <- read.csv("../data/htseq/bulk_VZ_CP_from_ATAC/Exprs_HTSCexon.csv"
                     , row.names = 1)

scExDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon.csv", row.names = 1)

# Picard Sequencing Statistics
picStatsDF <- read.csv("../metadata/PicardToolsQC.csv")

################################################################################

# Remove ERCCs
buExDatDF <- head(buExDatDF, -97)
scExDatDF <- head(scExDatDF, -97)
################################################################################

# Average counts for bulk for each transcript and counts of scRNAseq for each
# transript

# Pool scRNAseq
pScEx <- apply(scExDatDF, 1, sum)
head(pScEx)
tail(pScEx)

# Average expression for bulk for each gene
mBuEx <- apply(buExDatDF, 1, mean)

# Median for bulk for each gene
mdBuEx <- apply(buExDatDF, 1, median)
head(mdBuEx)
tail(mdBuEx)

cor(pScEx, mBuEx, method = "spearman")
cor(pScEx, mBuEx, method = "pearson")

cor(mdPdScEx, mBuEx, method = "spearman")
cor(mdPdScEx, mBuEx, method = "pearson")

cor(mnPdScEx, mBuEx, method = "spearman")
cor(mnPdScEx, mBuEx, method = "pearson")

# Log2 (counts + 1) Pooled scRNAseq vs bulk
plot(log(pScEx + 1, 2), log(mBuEx + 1,2))
plot(mdPdScEx, mBuEx)
plot(mnPdScEx, mBuEx, xlim = c(0,100000), ylim = c(0,50000))
################################################################################

# Read depth normalize bulk and scRNAseq

# Filter for number of mapped reads > 1.5*10^6
# Identify Cells IDs with number of mapped reads > 1.5*10^6
pfCells <- picStatsDF$X[picStatsDF$PF_READS_ALIGNED > 1.5*10^6]
# Filter expression dataframe for cell IDs with number of mapped reads > 1.5*10^6
ftScExDatDF <- scExDatDF[ ,colnames(scExDatDF) %in% pfCells]

# Pool scRNAseq
pScExDF <- data.frame(Counts = apply(ftScExDatDF, 1, sum))
head(pScExDF)

set.seed(11)
ncol(ftScExDatDF)
rNums <- sample(1:130, 130, replace = FALSE)
rndmGroups <- split(rNums, ceiling(seq_along(rNums) / (length(rNums) / 5)))
pScExDF <- data.frame(lapply(rndmGroups
                        , function (group) {apply(ftScExDatDF[ ,group], 1, sum)}))
head(pScExDF, 20)

# Read depth normalize pooled scRNAseq
rDep <- (apply(pScExDF, 2, sum) / 10^6)
pScRdnExDatDF <- pScExDF / rDep
head(pScRdnExDatDF)

# Read depth normalize bulk
rDep <- (apply(buExDatDF, 2, sum) / 10^6)
buRdNExDatDF <- buExDatDF / rDep

# Boxplot read depth normalized bulk log2 (counts)
boxplot(log(data.frame(buRdNExDatDF, pScRdnExDatDF) + 1, 2), range = 0)

# Mean counts for bulk RNAseq
mBuEx <- apply(buRdNExDatDF, 1, mean)

# Median counts for bulk RNAseq
mdBuEx <- apply(buRdNExDatDF, 1, median)

# Mean counts for pooled scRNAseq groups
mPdScEx <- apply(pScRdnExDatDF, 1, mean)

# Median counts for pooled scRNAseq groups
mdPdScEx <- apply(pScRdnExDatDF, 1, median)

cor(pScExDF, mBuEx, method = "spearman")
cor(pScExDF, mBuEx, method = "pearson")

plot(pScExDF[ ,1], mBuEx)

plot(log(mdPdScEx, 2), log(mdBuEx, 2))
cor(mdPdScEx, mdBuEx, method = "spearman")
cor(mdPdScEx, mdBuEx, method = "pearson")

plot(log(mPdScEx, 2), log(mBuEx, 2))
cor(mPdScEx, mBuEx, method = "spearman")
cor(mPdScEx, mBuEx, method = "pearson")

cor(data.frame(pScRdnExDatDF, buRdNExDatDF), method = "spearman")

######

# Plot read depth normalized pooled scRNAseq vs bulk

# Log2 bulk vs pooled scRNAseq
ggDF <- data.frame(scRNAseq = log(mPdScEx + 1, 2)
                   , Mean_Bulk_RNAseq = log(mBuEx + 1, 2))
ggplot(ggDF, aes(x = scRNAseq, y = Mean_Bulk_RNAseq)) +
  geom_point(shape = 1, alpha = 0.5) +
  stat_smooth(method = lm)
# + coord_cartesian(xlim = c(0, 500000), ylim = c(0, 2500) )    

# MA Plot
# Added +1 to all counts
# ***NOTE Inf values
ggDF <- data.frame(Avg = (0.5*log((pScExDF[ ,1] + 1) * (mBuEx + 1), 2))
                   , Log2Ratio = log(((pScExDF[ ,1] + 1) / (mBuEx + 1)), 2))
head(ggDF)
ggplot(ggDF, aes(x = Avg, y = Log2Ratio)) +
  geom_point(shape = 1, alpha = 0.5) +
  stat_smooth(method = lm)
# + coord_cartesian(xlim = c(0, 500000), ylim = c(0, 2500))

# MA Plot
# ***NOTE Inf values
ggDF <- data.frame(Avg = (0.5*log(mPdScEx * mBuEx, 2))
                   , Log2Ratio = log((mPdScEx / mBuEx), 2))
head(ggDF)
ggplot(ggDF, aes(x = Avg, y = Log2Ratio)) +
  geom_point(shape = 1, alpha = 0.5) +
  stat_smooth(method = lm)

# Filter MA plot for A < 0
Avg <- (0.5*log((mPdScEx) * (mBuEx), 2))
ftPdDF <- mPdScEx[Avg > 0]
ftBuDF <- mBuEx[Avg > 0]
ggDF <- data.frame(Avg = (0.5*log((ftPdDF + 1) * (ftBuDF + 1), 2))
                   , Log2Ratio = log(((ftPdDF + 1) / (ftBuDF + 1)), 2))
head(ggDF)
ggplot(ggDF, aes(x = Avg, y = Log2Ratio)) +
  geom_point(shape = 1, alpha = 0.5) +
  stat_smooth(method = lm)

log(4, 2)



# Correlation across individual bulk samples and pooled scRNAseq
df <- data.frame(buRdNExDatDF, pScRdnExDatDF)
cor(df, method = "spearman")




# Standard deviation for bulk RNAseq
sdBuEx <- apply(buRdNExDatDF, 1, sd)
