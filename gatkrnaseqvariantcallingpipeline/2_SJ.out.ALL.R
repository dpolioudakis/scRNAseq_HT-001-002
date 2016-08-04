## This script creates a combined splice junction file for all samples for 2pass Star alignment

setwd("/geschwindlabshares/eQTL/Round2/2015-9086R/data/BAM") ##set directory to where the first STAR alignment occured
options(stringsAsFactors=FALSE)
##get all "SJ.out.tab files"
files <- list.files(pattern = "SJ.out.tab",recursive = TRUE)
##Make merged list of all novel splice junction 
#loop through all files, merge two files, then unique
file1=paste("./",files[1],sep="")
SJ.combine=read.delim(file1,header=FALSE,sep="\t",quote = "")
for (i in 2:length(files)){
  file2=paste("./",files[i],sep="")
  SJ=read.delim(file2,header=FALSE,sep="\t",quote = "")
  SJ.merge=rbind(SJ.combine,SJ)
  SJ.combine=unique(SJ.merge)
}
attach(SJ.combine)
new <- SJ.combine[order(V1),] 

setwd("/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code") ## set output directory
write.table(new,col.names=FALSE,row.names=FALSE,file="SJ.all.out.tab",sep="\t",quote=FALSE) ##name output file
