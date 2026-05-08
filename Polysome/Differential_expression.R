library(tximport)
library(DESeq2)
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)
library(ggplot2)

rm(list=ls())
setwd(dir = "~")

#Loading the files
files <- paste("a",1:42,".genes.results",sep="")
files<-files[-which(grepl("a22",files))]
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
txi.rsem$length[txi.rsem$length == 0] <- 1

#Design
sample<-paste("a",1:42,sep="")
fraction<-factor(rep(c("F5","F10","F11","F12","F15","F16","F17"),6),
                  levels=c("F5","F10","F11","F12","F15","F16","F17"))

condition<-factor(c(rep("WT",14),rep("L100Q",14),rep("GFP",14)),levels=c("GFP","WT","L100Q"))
sampleInfo<-data.frame(sample,fraction,condition)
sampleInfo<-sampleInfo[which(sampleInfo$sample!="a22"),]
sampleInfo

#Reading into DESEQ2
deseqdata <- DESeqDataSetFromTximport(txi.rsem, colData = sampleInfo, design = ~fraction+condition:fraction)
deseqdata

#checking the factor levels
deseqdata$condition
deseqdata$fraction

# pre-filtering- Averaging at 5 per sample
deseqdata <- estimateSizeFactors(deseqdata)
nc <- counts(deseqdata, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 6
deseqdata <- deseqdata[filter,] #18293

#Running DESEQ
deseqdata <- estimateDispersions(deseqdata)
deseqdata <- nbinomWaldTest(deseqdata)

#Looking at the variability between replicates
vsd <- vst(deseqdata)
m<-plotPCA(vsd, "condition")
plotPCA(vsd, "fraction")
plotPCA(vsd, c("condition","fraction"))

### Comparing fold change in mutant vs wildtype
resultsNames(deseqdata)

##
rs_names<-resultsNames(deseqdata)
rs_names<-rs_names[grepl("condition",rs_names)]
Values<-list()
for(i in rs_names){
  df<-data.frame(results(deseqdata,name=i))
  df<-df[c("log2FoldChange","padj")]
  name<-gsub("fraction|condition","",i)
  names(df)<-paste(name,names(df),sep=".")
  df$Genes<-rownames(df)
  Values[[name]]<-df
}

LFC<-Values %>% purrr::reduce(left_join, by = "Genes")
rownames(LFC)<-LFC$Genes
LFC<-LFC[,-3]
write.table(LFC,"../LFC_by_conditions.txt",sep="\t",row.names = T,col.names = T,quote = F)

#output TPM
TPM<-data.frame(txi.rsem$abundance)
colnames(TPM)<-files
colnames(TPM)<-sapply(colnames(TPM),function(x){strsplit(x,"[.]")[[1]][1]})

#Getting the names
sample<-paste("a",1:42,sep="")
fraction<-factor(rep(c("F5","F10","F11","F12","F15","F16","F17"),6),
                 levels=c("F5","F10","F11","F12","F15","F16","F17"))

condition<-factor(c(rep("WT",14),rep("L100Q",14),rep("GFP",14)),levels=c("GFP","WT","L100Q"))
sampleInfo<-data.frame(sample,fraction,condition)
sampleInfo$ID<-paste(sampleInfo$fraction,sampleInfo$condition,sep="_")
sampleInfo$ID<-paste(sampleInfo$ID,rep(c(rep(1,7),rep(2,7)),3),sep="_")
sampleInfo<-sampleInfo[which(sampleInfo$sample!="a22"),]
sampleInfo

table(colnames(TPM)==sampleInfo$sample)
colnames(TPM)<-sampleInfo$ID
TPM$Genes<-rownames(TPM)
write.table(TPM,"../TPM.txt",sep="\t",row.names = F,col.names = T,quote = F)

#Normalized counts
Norm_counts<-counts(deseqdata, normalized=TRUE)
colnames(Norm_counts)<-files
colnames(Norm_counts)<-sapply(colnames(Norm_counts),function(x){strsplit(x,"[.]")[[1]][1]})

#Getting the names
sample<-paste("a",1:42,sep="")
fraction<-factor(rep(c("F5","F10","F11","F12","F15","F16","F17"),6),
                 levels=c("F5","F10","F11","F12","F15","F16","F17"))

condition<-factor(c(rep("WT",14),rep("L100Q",14),rep("GFP",14)),levels=c("GFP","WT","L100Q"))
sampleInfo<-data.frame(sample,fraction,condition)
sampleInfo$ID<-paste(sampleInfo$fraction,sampleInfo$condition,sep="_")
sampleInfo$ID<-paste(sampleInfo$ID,rep(c(rep(1,7),rep(2,7)),3),sep="_")
sampleInfo<-sampleInfo[which(sampleInfo$sample!="a22"),]
sampleInfo

table(colnames(Norm_counts)==sampleInfo$sample)
colnames(Norm_counts)<-sampleInfo$ID
Norm_counts<-data.frame(Norm_counts)
Norm_counts$Genes<-rownames(Norm_counts)
write.table(Norm_counts,"../Norm_counts.txt",sep="\t",row.names = F,col.names = T,quote = F)








