## SumExpr.R
## Script to sum all expression data from cufflinks and htseq counts,
## as well as remove duplicates from cufflinks data

## To be used after the Orion RNA Seq Pipeline, before CD_RNASeq.R

####### Run this code on the server (recommend using a qsub, the for loop takes awhile)

### Use 6.5_Move.sh to move all individual files to compile folder

## Edit substr length to remove _(gene|exon)_union_count from file names and make column names

options(stringsAsFactors =F)
# setwd("../data/aligned/HTSCgene")

## Now HTSeq Counts, union exon

rm(list=ls())
setwd("../data/htseq/bulk_VZ_CP_from_ATAC")
algn = vector(mode = "list")
algn$Sample.ID = list.files(".","*exon_union_count")

Refrc=read.table(file=paste("./", algn$Sample.ID[1],sep=""),head=F)
ensembl =  Refrc[,1]
datExpr = data.frame(Sample1= Refrc[,2])
rownames(datExpr) = ensembl

nsample=length(algn$Sample.ID)
for (i in c(2:nsample)){
  temp=read.table(file=paste(algn$Sample.ID[i],sep=""),head=F)
  k=match(ensembl,temp[,1])
  print(paste(algn$Sample.ID[i],length(k[!(is.na(k))]),sep="       "))
  datExpr=cbind(datExpr,temp[k,2])
}

colnames(datExpr) <- gsub("_.*", "", as.character(algn$Sample.ID))
write.csv(datExpr, file="Exprs_HTSCexon.csv")

# ## Now HTSeq Counts, union gene

# rm(list=ls())
# algn = vector(mode = "list")
# algn$Sample.ID = list.files(".","*gene_union_count")

# Refrc=read.table(file=paste("./", algn$Sample.ID[1],sep=""),head=F)
# ensembl =  Refrc[,1]
# datExpr = data.frame(Sample1= Refrc[,2])
# rownames(datExpr) = ensembl

# nsample=length(algn$Sample.ID)
# for (i in c(2:nsample)){
  # temp=read.table(file=paste(algn$Sample.ID[i],sep=""),head=F)
  # k=match(ensembl,temp[,1])
  # print(paste(algn$Sample.ID[i],length(k[!(is.na(k))]),sep="       "))
  # datExpr=cbind(datExpr,temp[k,2])
# }

# colnames(datExpr)=substr(as.character(algn$Sample.ID),1,10)
# write.csv(datExpr, file="Exprs_HTSCgene.csv")


# ## First, Cufflinks
# 
# algn = vector(mode = "list")
# algn$Sample.ID = list.files(".","*.fpkm_tracking")
# 
# Refrc=read.table(file=paste("./", algn$Sample.ID[1],sep=""),head=T)
# ensembl =  Refrc$gene_id
# datExpr = data.frame(Sample1= Refrc[,10])
# rownames(datExpr) = make.names(ensembl, unique=TRUE)
# 
# nsample=length(algn$Sample.ID)
# for (i in c(2:nsample)){
#   temp=read.table(file=paste(algn$Sample.ID[i],sep=""),head=T)
#   temp_names=temp$gene_id
#   rownames(temp)=make.names(temp_names,unique=TRUE)
#   k=match(rownames(datExpr),rownames(temp))
#   print(paste(algn$Sample.ID[i],length(k[!(is.na(k))]),sep="       "))
#   datExpr=cbind(datExpr,temp[k,"FPKM"])
# }
# 
# colnames(datExpr)=substr(as.character(algn$Sample.ID),1,5)
# write.table(datExpr, "Wave3_FPKM.txt", sep="\t")
# 
# ## Next, remove duplicates from Cufflinks data
# 
# cuff = datExpr
# rm(datExpr)
# 
# match = cuff[grep("\\.1$",rownames(cuff)),]
# match_names = substr(as.character(rownames(match)),1,15)
# samples = colnames(cuff)
# 
# cuffsum = data.frame(matrix(NA,nrow=length(match_names),ncol=length(samples)))
# rownames(cuffsum)=match_names
# colnames(cuffsum)=samples
# 
# temp = data.frame()
# 
# ## Make a dataframe with sums of duplicate rows in cuff
# 
# for (i in 1:length(rownames(cuff)))
# {
#   for (j in 1:length(match_names))
#   {
#     if(rownames(cuff)[i]==match_names[j])
#     {
#       matchy_match = match_names[j]
#       temp=cuff[grep(as.character(matchy_match),rownames(cuff)),]
#       for (k in 1:length(samples))
#       {
#           temp_sum = sum(temp[,k])
#           cuffsum[j,k] = temp_sum
#       }
#     }
#   }
#   print(paste("On Row",i,"of",length(rownames(cuff)),sep=" "))
# }
# 
# ## Now substitute summed duplicate rows into rmdupcuff
# 
# rmdupcuff = cuff
# sub = which(rownames(rmdupcuff) %in% rownames(cuffsum))
# rmdupcuff[sub,] = cuffsum
# 
# ## Now remove duplicate rows from rmdupcuff
# 
# dup=cuff[grep("*[.]",rownames(cuff)),]
# length = length(rownames(rmdupcuff))
# 
# to_remove = which(rownames(rmdupcuff) %in% rownames(dup))
# rmdupcuff = rmdupcuff[-to_remove,]
# 
# ## rmdupcuff should have 63654 genes/rows left
# 
# rmdupcuff[grep("*[.]",rownames(rmdupcuff)),] ## Should be 0
# dim(rmdupcuff) ## Should be 63654 47
# 
# write.csv(rmdupcuff, file="Wave3exprs_Cuff.csv")





