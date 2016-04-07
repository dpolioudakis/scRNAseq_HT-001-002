
# Damon Polioudakis
# 2016-02-18
# Adopted from Luis Ubieta's script
# Check overlap of genes expressed in our Fluidigm medium HT VZ scRNA-seq
# samples to Kriegstein's 2015 scRNA-seq cell types from VZ
################################################################################

rm(list=ls())

library(xlsx)
library(WGCNA)

##Load Pollen et al, 2015 dataset (Single-cell transcriptomes of human fetal cortex)
#http://www.ncbi.nlm.nih.gov/pubmed/?term=26406371
PollenfnS1="../../kriegstein_2015/analysis/tables/TableS1.xlsx"
PollenfnS2="../../kriegstein_2015/analysis/tables/TableS2.xlsx"
PollenfnS3="../../kriegstein_2015/analysis/tables/TableS3.xlsx"
PollenfnS4="../../kriegstein_2015/analysis/tables/TableS4.csv"

PollenS3=read.xlsx(PollenfnS3,2) #oRG vs vRG genesets
PollenS4=read.csv(PollenfnS4)    #Correlations to classified cell types

##Convert PollenS3 to ENSG and output:

options(stringsAsFactors=FALSE);
library(biomaRt)
getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","strand","band","gene_biotype")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org") ## Using Gencode v19 annotations
geneAnno <- getBM(attributes = getinfo,filters=c("chromosome_name"),values=c(seq(1,22,by=1),"X","Y","M"),mart=mart)


PollenS3Genes=list(PollenS3$oRG,PollenS3$vRG,PollenS3$RG,PollenS3$UpORG,PollenS3$UpVRG);
names(PollenS3Genes)=c("oRG","vRG","RG","UpORG","UpVRG");
PollenS3GenesENSG=list();

for (i in 1:length(PollenS3Genes)){
    ind=match(PollenS3Genes[[i]],geneAnno$hgnc_symbol)
    PollenS3GenesENSG[[i]]=geneAnno[ind,1]
    PollenS3GenesENSG[[i]]=PollenS3GenesENSG[[i]][!is.na(PollenS3GenesENSG[[i]])]
    write.table(PollenS3GenesENSG[[i]],file=paste("TableS3",names(PollenS3Genes[i]),"ENSG.txt",sep="-"),row.names=FALSE, col.names=FALSE, quote=FALSE);
}

names(PollenS3GenesENSG)=c("oRG","vRG","RG","UpORG","UpVRG");
dir.create("../../kriegstein_2015/analysis/tables/TableS3-CellTypesGenes-ENSG")
save(PollenS3GenesENSG, file="../../kriegstein_2015/analysis/tables/TableS3-CellTypesGenes-ENSG/PollenS3GenesENSG.R")

##Convert correlations to pvalues in Pollen data.

setwd("../../kriegstein_2015/analysis/tables/")

#Remove negative correlations first

RGkeep= which(PollenS4$Correlation.RG>0);
RG=PollenS4[RGkeep,c(1,6)];

IPCkeep= which(PollenS4$Correlation.IPC>0);
IPC=PollenS4[IPCkeep,c(1,7)];

Neuronkeep= which(PollenS4$Correlation.Neuron>0);
Neuron=PollenS4[Neuronkeep,c(1,8)];

INkeep= which(PollenS4$Correlation.Interneuron>0);
Interneuron=PollenS4[INkeep,c(1,9)];

vRGkeep= which(PollenS4$Correlation.vRG>0);
vRG=PollenS4[vRGkeep,c(1,10)];

oRGkeep= which(PollenS4$Correlation.oRG>0);
oRG=PollenS4[oRGkeep,c(1,11)];

#Convert remaining correlations to pvalues in Pollen data.

RG$pvalue=corPvalueStudent(RG$Correlation.RG, 393);
IPC$pvalue=corPvalueStudent(IPC$Correlation.IPC, 393);
Neuron$pvalue=corPvalueStudent(Neuron$Correlation.Neuron, 393);
Interneuron$pvalue=corPvalueStudent(Interneuron$Correlation.Interneuron, 393);
vRG$pvalue=corPvalueStudent(vRG$Correlation.vRG, 393);
oRG$pvalue=corPvalueStudent(oRG$Correlation.oRG, 393);

#Calculate FDR adjusted pvalue

RG$FDRpvalue=p.adjust(RG$pvalue, method="fdr", n=length(RG$pvalue));
IPC$FDRpvalue=p.adjust(IPC$pvalue, method="fdr", n=length(IPC$pvalue));
Neuron$FDRpvalue=p.adjust(Neuron$pvalue, method="fdr", n=length(Neuron$pvalue));
Interneuron$FDRpvalue=p.adjust(Interneuron$pvalue, method="fdr", n=length(Interneuron$pvalue));
vRG$FDRpvalue=p.adjust(vRG$pvalue, method="fdr", n=length(vRG$pvalue));
oRG$FDRpvalue=p.adjust(oRG$pvalue, method="fdr", n=length(oRG$pvalue));

##Generate genelists FDR<0.05 and output

dir.create("../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes")

RGkeep= which(RG$FDRpvalue<=0.05); length(RGkeep)
RGf=RG[RGkeep,];
write.csv(RG,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-RG-Unfiltered.csv")
write.csv(RGf,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-RG-FDRFiltered.csv")

IPCkeep= which(IPC$FDRpvalue<=0.05); length(IPCkeep)
IPCf=IPC[IPCkeep,];
write.csv(IPC,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-IPC-Unfiltered.csv")
write.csv(IPCf,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-IPC-FDRFiltered.csv")

Neuronkeep= which(Neuron$FDRpvalue<=0.05); length(Neuronkeep)
Neuronf=Neuron[Neuronkeep,];
write.csv(Neuron,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-Neuron-Unfiltered.csv")
write.csv(Neuronf,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-Neuron-FDRFiltered.csv")

Interneuronkeep= which(Interneuron$FDRpvalue<=0.05); length(Interneuronkeep)
Interneuronf=Interneuron[Interneuronkeep,];
write.csv(Interneuron,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-Interneuron-Unfiltered.csv")
write.csv(Interneuronf,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-Interneuron-FDRFiltered.csv")

oRGkeep= which(oRG$FDRpvalue<=0.05); length(oRGkeep)
oRGf=oRG[oRGkeep,];
write.csv(oRG,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-oRG-Unfiltered.csv")
write.csv(oRGf,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-oRG-FDRFiltered.csv")

vRGkeep= which(vRG$FDRpvalue<=0.05); length(vRGkeep)
vRGf=vRG[vRGkeep,];
write.csv(vRG,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-vRG-Unfiltered.csv")
write.csv(vRGf,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-vRG-FDRFiltered.csv")

##Convert from Gene_ID to ENSG IDs, remove genes without ENSG


options(stringsAsFactors=FALSE);
library(biomaRt)
getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","strand","band","gene_biotype")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org") ## Using Gencode v19 annotations
geneAnno <- getBM(attributes = getinfo,filters=c("chromosome_name"),values=c(seq(1,22,by=1),"X","Y","M"),mart=mart)


ind=match(RGf$Gene_ID,geneAnno$hgnc_symbol);
RGf$ENSG_ID=geneAnno[ind,1];
RGf=RGf[!is.na(RGf$ENSG_ID),];
write.csv(RGf,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-RG-FDRFiltered-ENSG.csv")
      
ind=match(IPCf$Gene_ID,geneAnno$hgnc_symbol);
IPCf$ENSG_ID=geneAnno[ind,1];
IPCf=IPCf[!is.na(IPCf$ENSG_ID),];
write.csv(IPCf,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-IPC-FDRFiltered-ENSG.csv")      
      
ind=match(Neuronf$Gene_ID,geneAnno$hgnc_symbol);
Neuronf$ENSG_ID=geneAnno[ind,1];
Neuronf=Neuronf[!is.na(Neuronf$ENSG_ID),];
write.csv(Neuronf,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-Neuron-FDRFiltered-ENSG.csv")      
      
ind=match(Interneuronf$Gene_ID,geneAnno$hgnc_symbol);
Interneuronf$ENSG_ID=geneAnno[ind,1];
Interneuronf=Interneuronf[!is.na(Interneuronf$ENSG_ID),];
write.csv(Interneuronf,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-Interneuron-FDRFiltered-ENSG.csv")      

ind=match(oRGf$Gene_ID,geneAnno$hgnc_symbol);
oRGf$ENSG_ID=geneAnno[ind,1];
oRGf=oRGf[!is.na(oRGf$ENSG_ID),];
write.csv(oRGf,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-oRG-FDRFiltered-ENSG.csv")      

ind=match(vRGf$Gene_ID,geneAnno$hgnc_symbol);
vRGf$ENSG_ID=geneAnno[ind,1];
vRGf=vRGf[!is.na(vRGf$ENSG_ID),];
write.csv(vRGf,"../../kriegstein_2015/analysis/tables/TableS4-CellTypesGenes/TableS4-vRG-FDRFiltered-ENSG.csv")

##Create gene lists exlusive to each cell type; ie exclude genes overlapping any other cell type

oRGu=setdiff(oRGf$ENSG_ID,c(vRGf$ENSG_ID,Neuronf$ENSG_ID,Interneuronf$ENSG_ID,IPCf$ENSG_ID));
vRGu=setdiff(vRGf$ENSG_ID,c(oRGf$ENSG_ID,Neuronf$ENSG_ID,Interneuronf$ENSG_ID,IPCf$ENSG_ID));
IPCu=setdiff(IPCf$ENSG_ID,c(oRGf$ENSG_ID,vRGf$ENSG_ID,Neuronf$ENSG_ID,Interneuronf$ENSG_ID));
Neuronu=setdiff(Neuronf$ENSG_ID,c(vRGf$ENSG_ID,oRGf$ENSG_ID,Interneuronf$ENSG_ID,IPCf$ENSG_ID));
Interneuronu=setdiff(Interneuronf$ENSG_ID,c(vRGf$ENSG_ID,oRGf$ENSG_ID,Neuronf$ENSG_ID,IPCf$ENSG_ID));

oRGu=oRGf[match(oRGu,oRGf$ENSG_ID),]
vRGu=vRGf[match(vRGu,vRGf$ENSG_ID),]
IPCu=IPCf[match(IPCu,IPCf$ENSG_ID),]
Neuronu=Neuronf[match(Neuronu,Neuronf$ENSG_ID),]
Interneuronu=Interneuronf[match(Interneuronu,Interneuronf$ENSG_ID),]

#Calculate overlaps between PollenS3 and and Human-Specific Genes (Overlaps in Cell Type Genes allowed) 

setwd("../../kriegstein_2015/analysis/tables//TableS3-CellTypesGenes-ENSG");
load(file="../../kriegstein_2015/analysis/tables//TableS3-CellTypesGenes-ENSG/PollenS3GenesENSG.R");

CellTypesGenes=PollenS3GenesENSG
HumanEnhGenes=list(ATACSeqPC[[1]][,1],ATACSeqPC[[2]][,1],ATACSeqPC[[3]][,1],ATACSeqPC[[4]][,1],ATACSeqPC[[5]][,1],ATACSeqPC[[6]][,1],ATACSeqPC[[7]][,1],ATACSeqPC[[8]][,1],List2$V1,List3$V1,VZboth,CPboth,ATACSeqPC[[9]][,1],ATACSeqPC[[10]][,1],ATACSeqPC[[11]][,1],ATACSeqPC[[12]][,1],ASDWGS=List16PC[[1]]);
names(HumanEnhGenes)=c("ATAC.VZ.All","ATAC.CP.All","ATAC.VZ.Prox","ATAC.CP.Prox","ATAC.VZ.Dist","ATAC.CP.Dist","ATAC.Prox","ATAC.Dist","HiC.VZ","HiC.CP","VZboth","CPboth","VZboth.Dist-Strict","CPboth.Dist-Strict","VZboth.All-Strict","CPboth.All-Strict","ASDWGS");

pvalues=matrix(ncol=length(CellTypesGenes),nrow=length(HumanEnhGenes));
OR=matrix(ncol=length(CellTypesGenes),nrow=length(HumanEnhGenes));

for (j in 1:length(HumanEnhGenes)){ 
    for (i in 1:length(CellTypesGenes)){
        Q=length(intersect(HumanEnhGenes[[j]],CellTypesGenes[[i]]))
        M=length(HumanEnhGenes[[j]])
        K=length(CellTypesGenes[[i]])
        T=length(Background$SourceIdentifier)
        Fishers=fisher.test(matrix(c(Q,M-Q,K-Q,T-M-K+Q),2,2));
        pvalues[j,i]=Fishers$p.value
        OR[j,i]=Fishers$estimate
        write.table(intersect(HumanEnhGenes[[j]],CellTypesGenes[[i]]),file=paste(names(CellTypesGenes[i]),names(HumanEnhGenes[j]),"humanreg.txt",sep="-"),row.names=FALSE,col.names=FALSE,quote=FALSE)
    }
}

rownames(pvalues)=names(HumanEnhGenes)
colnames(pvalues)=names(CellTypesGenes)

rownames(OR)=names(HumanEnhGenes)
colnames(OR)=names(CellTypesGenes)

output=cbind(pvalues,OR)
outputnames=c(paste("pval.",colnames(OR),sep=""),paste("OR.",colnames(OR),sep=""))
colnames(output)=outputnames
write.csv(output,file="HumanReg-PollenS3-Overlaps-Fishers.csv")

##Calculate overlaps between Pollen et al and Human-Specific Genes (Overlaps in Cell Type Genes allowed) 

setwd("/geschwindlabshares/atacseq/HumanEnhancers/Pollen-2015-SingleCell/HumanRegGenesbyCellType/Non-Unique")
load(file="../../kriegstein_2015/analysis/tables/HumanRegGenesbyCellType/Non-Unique/PollenCellTypesGenes.R")      

#CellTypesGenes=list(RGf=RGf$ENSG_ID,IPCf=IPCf$ENSG_ID,Neuronf=Neuronf$ENSG_ID,Interneuronf=Interneuronf$ENSG_ID,vRGf=vRGf$ENSG_ID,oRGf=oRGf$ENSG_ID);
#save(CellTypesGenes, file="/geschwindlabshares/atacseq/HumanEnhancers/Pollen-2015-SingleCell/HumanRegGenesbyCellType/Non-Unique/PollenCellTypesGenes.R")

HumanEnhGenes=list(ATACSeqPC[[1]][,1],ATACSeqPC[[2]][,1],ATACSeqPC[[3]][,1],ATACSeqPC[[4]][,1],ATACSeqPC[[5]][,1],ATACSeqPC[[6]][,1],ATACSeqPC[[7]][,1],ATACSeqPC[[8]][,1],List2$V1,List3$V1,VZboth,CPboth,ATACSeqPC[[9]][,1],ATACSeqPC[[10]][,1],ATACSeqPC[[11]][,1],ATACSeqPC[[12]][,1],ASDWGS=List16PC[[1]]);
names(HumanEnhGenes)=c("ATAC.VZ.All","ATAC.CP.All","ATAC.VZ.Prox","ATAC.CP.Prox","ATAC.VZ.Dist","ATAC.CP.Dist","ATAC.Prox","ATAC.Dist","HiC.VZ","HiC.CP","VZboth","CPboth","VZboth.Dist-Strict","CPboth.Dist-Strict","VZboth.All-Strict","CPboth.All-Strict", "ASDWGS");

pvalues=matrix(ncol=length(CellTypesGenes),nrow=length(HumanEnhGenes));
OR=matrix(ncol=length(CellTypesGenes),nrow=length(HumanEnhGenes));

for (j in 1:length(HumanEnhGenes)){ 
    for (i in 1:length(CellTypesGenes)){
        Q=length(intersect(HumanEnhGenes[[j]],CellTypesGenes[[i]]))
        M=length(HumanEnhGenes[[j]])
        K=length(CellTypesGenes[[i]])
        T=length(Background$SourceIdentifier)
        Fishers=fisher.test(matrix(c(Q,M-Q,K-Q,T-M-K+Q),2,2));
        pvalues[j,i]=Fishers$p.value
        OR[j,i]=Fishers$estimate
        write.table(intersect(HumanEnhGenes[[j]],CellTypesGenes[[i]]),file=paste(names(CellTypesGenes[i]),names(HumanEnhGenes[j]),"humanreg.txt",sep="-"),row.names=FALSE,col.names=FALSE,quote=FALSE)
    }
}

rownames(pvalues)=names(HumanEnhGenes)
colnames(pvalues)=names(CellTypesGenes)

rownames(OR)=names(HumanEnhGenes)
colnames(OR)=names(CellTypesGenes)

output=cbind(pvalues,OR)
outputnames=c(paste("pval.",colnames(OR),sep=""),paste("OR.",colnames(OR),sep=""))
colnames(output)=outputnames
write.csv(output,file="HumanReg-Pollen-Overlaps-Fishers.csv")

##Calculate overlaps between Pollen et al and Human-Specific Genes (Cell Type unique list)

setwd("/geschwindlabshares/atacseq/HumanEnhancers/Pollen-2015-SingleCell/HumanRegGenesbyCellType/Unique")
load(file="/geschwindlabshares/atacseq/HumanEnhancers/Pollen-2015-SingleCell/HumanRegGenesbyCellType/Unique/PollenCellTypesGenes-Unique.R")      

#CellTypesGenes=list(IPC=IPCu$ENSG_ID,Neuron=Neuronu$ENSG_ID,Interneuron=Interneuronu$ENSG_ID,vRG=vRGu$ENSG_ID,oRG=oRGu$ENSG_ID);
#save(CellTypesGenes, file="/geschwindlabshares/atacseq/HumanEnhancers/Pollen-2015-SingleCell/HumanRegGenesbyCellType/Unique/PollenCellTypesGenes-Unique.R")

HumanEnhGenes=list(ATACSeqPC[[1]][,1],ATACSeqPC[[2]][,1],ATACSeqPC[[3]][,1],ATACSeqPC[[4]][,1],ATACSeqPC[[5]][,1],ATACSeqPC[[6]][,1],ATACSeqPC[[7]][,1],ATACSeqPC[[8]][,1],List2$V1,List3$V1,VZboth,CPboth,ATACSeqPC[[9]][,1],ATACSeqPC[[10]][,1],ATACSeqPC[[11]][,1],ATACSeqPC[[12]][,1]);
names(HumanEnhGenes)=c("ATAC.VZ.All","ATAC.CP.All","ATAC.VZ.Prox","ATAC.CP.Prox","ATAC.VZ.Dist","ATAC.CP.Dist","ATAC.Prox","ATAC.Dist","HiC.VZ","HiC.CP","VZboth","CPboth","VZboth.Dist-Strict","CPboth.Dist-Strict","VZboth.All-Strict","CPboth.All-Strict");

pvaluesu=matrix(ncol=length(CellTypesGenes),nrow=length(HumanEnhGenes))
ORu=matrix(ncol=length(CellTypesGenes),nrow=length(HumanEnhGenes))  

for (j in 1:length(HumanEnhGenes)){ 
    for (i in 1:length(CellTypesGenes)){
        Q=length(intersect(HumanEnhGenes[[j]],CellTypesGenes[[i]]))
        M=length(HumanEnhGenes[[j]])
        K=length(CellTypesGenes[[i]])
        T=length(Background$SourceIdentifier)
        Fishers=fisher.test(matrix(c(Q,M-Q,K-Q,T-M-K+Q),2,2));
        pvaluesu[j,i]=Fishers$p.value
        ORu[j,i]=Fishers$estimate
        write.table(intersect(HumanEnhGenes[[j]],CellTypesGenes[[i]]),file=paste(names(CellTypesGenes[i]),names(HumanEnhGenes[j]),"humanreg.txt",sep="-"),row.names=FALSE,col.names=FALSE,quote=FALSE)
    }
}

rownames(pvaluesu)=names(HumanEnhGenes)
colnames(pvaluesu)=names(CellTypesGenes)

rownames(ORu)=names(HumanEnhGenes)
colnames(ORu)=names(CellTypesGenes)

output=cbind(pvaluesu,ORu)
outputnames=c(paste("pval.",colnames(ORu),sep=""),paste("OR.",colnames(ORu),sep=""))
colnames(output)=outputnames
write.csv(output,file="HumanReg-PollenUnique-Overlaps-Fishers.csv")