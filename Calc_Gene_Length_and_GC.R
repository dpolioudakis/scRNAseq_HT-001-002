rm(list=ls())

source("http://bioconductor.org/biocLite.R")
# biocLite("Genominator")
library(Genominator)

options(stringsAsFactors=FALSE)

exDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon.csv")

## Get Gencode 18 gtf file - this was cleaned by selecting only the columns containing the word "exon" and with the relevant information - chrnum, feature type, strand, start, end, and ensg and ense IDs separated by a semicolon
gtfinfo <- read.table("../../source/gencode.v19.annotation.gtf",sep="\t")
keep <- gtfinfo[,3]=="exon" ## Keep only the exon level features
gtfinfo <- gtfinfo[keep,]

genexoninfo <- unlist(strsplit(gtfinfo[,9],"[;]")) ## Split the semicolon separated information
gen.col<- which(regexpr("gene_id ", genexoninfo)>0)  ## finding which has gene_id for ensembl ID
getseq= genexoninfo[gen.col]
ENSGID <- substr(getseq,9,100)
length(unique(ENSGID)) ##62069

trans.col=which(regexpr("transcript_id ", genexoninfo)>0)
transeq = genexoninfo[trans.col]
ENSEID <- substr(transeq,16,100)
length(unique(ENSEID)) #213272
gtfinfo <- cbind(gtfinfo[,c(1:8)],ENSGID,ENSEID)

gene.col=which(regexpr("gene_name ", genexoninfo)>0)
geneseq = genexoninfo[gene.col]
GENEID <- substr(geneseq,12,100)
length(unique(GENEID)) #52775


gtfinfo <- cbind(gtfinfo[,c(1:8)],ENSGID,ENSEID)
gtfinfo= gtfinfo[,-c(6,8)] ## 6 and 8 columns are blank

## Keep only one copy of each ENSEID - the gtf file records one copy for each transcript id
keep <- match(unique(ENSEID),ENSEID)
gtfinfo1 <- gtfinfo[keep,]
##gtfinfo[,1] <- substr(gtfinfo[,1],4,10) ## 672406 exons is exactly what biomaRt contains

## Recode things for the Genominator package
chrnums <- gtfinfo1[,1] ## Using as.factor to coerce chromosome names can really botch things up... beware! So go ahead and convert MT, X, and Y to numbers throughout, unless necessary for other purposes
chrnums[chrnums=="MT"] <- "20"
chrnums[chrnums=="X"] <- "21"
chrnums[chrnums=="Y"] <- "22"
# rmChR.col1=which(regexpr("HG", chrnums)>0)
# rmChR.col2= which(regexpr("GL", chrnums)>0) ## removing Non-annotated(NT) chromosomes
# rmChR.col3= which(regexpr("HS", chrnums)>0)
# rmChR.col=c(rmChR.col1,rmChR.col2,rmChR.col3)
# gtfinfo1=gtfinfo1[-rmChR.col,]
chrnums=chrnums[-rmChR.col]
gtfinfo1[,1] <- chrnums ## Check here

gtfinfo=gtfinfo1

strinfo <- gtfinfo[,6]
strinfo[strinfo=="+"] <- 1L
strinfo[strinfo=="-"] <- -1L
gtfinfo[,6] <- strinfo


geneDat1=gtfinfo[,c(1,6,4,5,7,8)] ## chr integer, strand integer (-1L,0L,1L), start integer, end integer, ensg and transcript id
geneDat1 <- data.frame(as.numeric(chrnums),as.numeric(geneDat1[,2]),as.numeric(geneDat1[,3]),as.numeric(geneDat1[,4]),geneDat1[,5],geneDat1[,6])
names(geneDat1) <- c("chr","strand","start","end","ensembl_gene_id","ensembl_exon_id")

geneDatX <- geneDat1[order(geneDat1[,1],geneDat1[,3]),]
#DP edit - remove NAs from ERCC chromosomes
geneDatX <- geneDatX[complete.cases(geneDatX),]
validAnnotation(geneDatX) ## Have genominator check if this is a valid data object
geneDatX <- makeGeneRepresentation(annoData=geneDat1,type="Ugene",gene.id = "ensembl_gene_id", transcript.id = "ensembl_exon_id",verbose=TRUE) ##should take few minutes !!!

save(geneDatX,file="../../source/GenominatorUnionGeneModelsENSEMBLhg19.rda")
load(file="../../source/GenominatorUnionGeneModelsENSEMBLhg19.rda")

## Now use the genominator output to calculate GC content ### Use mac laptop ###
library(Repitools)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
geneDat2 <- cbind(geneDatX,geneDatX[,3]-geneDatX[,2])
geneDat2 <- geneDat2[order(geneDat2[,5]),]

## Change formatting again
chrnums <- geneDat2[,"chr"]
chrnums[chrnums=="20"] <- "M" ## important as UCSC codes MT as M
chrnums[chrnums=="21"] <- "X"
chrnums[chrnums=="22"] <- "Y"
stinfo <- geneDat2[,"strand"]
stinfo[stinfo==c(-1)] <- "-"
stinfo[stinfo==c(1)] <- "+"

## Calculate GC content from hg19 using the union exon ranges
gcQuery <- GRanges(paste("chr", chrnums,sep=""),IRanges(geneDat2[,2],geneDat2[,3]),strand=stinfo) ## Convert to a genomic ranges object
gcContent <- gcContentCalc(x=gcQuery,organism=Hsapiens)

## Take a length weighted average of GC content percentages to get the GC content for the union gene model
head(geneDat2)
geneDat2 <- cbind(geneDat2,gcContent)
geneDat2 <- cbind(geneDat2,gcContent*geneDat2[,6])
unionGenes <- by(geneDat2[,6],as.factor(geneDat2[,5]),sum)
unionGC <- by(geneDat2[,8],as.factor(geneDat2[,5]),sum)
geneDat3 <- cbind(unionGenes,unionGC/unionGenes)
colnames(geneDat3) <- c("UnionExonLength","UnionGCcontent")
ENSEMBLhg19.70UnionAnno <- geneDat3

## Save for further usage
save(ENSEMBLhg19.70UnionAnno, file="../../source/ENSEMBLhg19_UnionAnno.rda")
load("../../source/ENSEMBLhg19_UnionAnno.rda")

# DP edit - calculate length bias for each cell

exDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon_CellIDs.csv")

unionGenes <- data.frame(Length = ENSEMBLhg19.70UnionAnno[ ,1])
exLenDF <- merge(x = exDatDF, y = unionGenes, by.x = "X", by.y = "row.names" )
avgLength <- apply(exLenDF, 2
                 , function(counts) sum(as.numeric(counts) * exLenDF["Length"]) / 
                   sum(as.numeric(counts)))
avgLength <- tail(head(avgLength, -1), -1)
avgGCdF <- merge(x = exDatDF, y = ENSEMBLhg19.70UnionAnno, by.x = "X", by.y = "row.names" )
avgGCdF <- avgGCdF[complete.cases(avgGCdF), ]
avgGC <- apply(avgGCdF, 2
                 , function(counts) sum(as.numeric(counts) * avgGCdF["UnionGCcontent"]) / 
                   sum(as.numeric(counts)))
avgGC <- tail(head(avgGC, -2), -1)
save(avgLength, avgGC, file = "../analysis/tables/Avg_Gene_Length_and_GC.rda")
