## This script creates a combined splice junction file for all samples for 2pass Star alignment

print("Starting script...")

outDir <- "/geschwindlabshares/RNAseq_singlecellfetal/HT-001-002_DP/data/bam/NoERCC/variant"

setwd("../data/bam/NoERCC") ##set directory to where the first STAR alignment occured
options(stringsAsFactors=FALSE)

##get all "SJ.out.tab files"
files <- Sys.glob('SxaQSEQsXap096L*/HT_*/SJ.out.tab')
# files <- list.files(pattern = "SJ.out.tab",recursive = TRUE)
print("Files to merge: ")
print(files)
##Make merged list of all novel splice junction 
#loop through all files, merge two files, then unique
file1 = paste("./", files[1], sep = "")
SJ.combine = read.delim(file1, header = FALSE, sep = "\t", quote = "")
print("Looping through files...")
for (i in 2:length(files)){
  # try(file2 = paste("./", files[i], sep=""))
 	file2 = files[i]
 	SJ = try(read.delim(file2, header = FALSE, sep = "\t", quote = ""))
	SJ.merge = rbind(SJ.combine,SJ)
    SJ.combine = unique(SJ.merge)
}
print("Done looping through files...")
attach(SJ.combine)
new <- SJ.combine[order(V1), ] 

print("Making outdir:")
print(outDir)
dir.create(outDir)
setwd(outDir) ## set output directory
write.table(new, col.names = FALSE, row.names = FALSE, file = "SJ.all.out.tab", sep = "\t", quote = FALSE) ##name output file
