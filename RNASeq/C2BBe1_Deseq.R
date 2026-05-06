library(tximport)
library(DESeq2)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)

rm(list=ls())
setwd(dir = "~")

#Whole cell
files<-paste0(34:42,".genes.results")
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
txi.rsem$length[txi.rsem$length == 0] <- 1

#Design
sample<-sapply(files,function(x){strsplit(x,"[.]")[[1]][1]})
condition<-factor(c("WT","WT","WT","L100Q","L100Q","L100Q","GFP","GFP","GFP"),levels=c("GFP","WT","L100Q"))
sampleInfo<-data.frame(sample,condition)
sampleInfo

#Reading into DESEQ2
deseqdata <- DESeqDataSetFromTximport(txi.rsem, colData = sampleInfo, design = ~condition)
deseqdata

# pre-filtering- Averaging at 5 per sample
keep <- rowSums(counts(deseqdata)) >= 45 #average of 5 per sample 
deseqdata <- deseqdata[keep,] #18147

#checking the factor levels
deseqdata$condition

#Running DESEQ
deseqdata<-DESeq(deseqdata)

#Looking at the variability between replicates
vsd <- vst(deseqdata)
plotPCA(vsd, "condition")

#results
resultsNames(deseqdata)


#WT
wt<-data.frame(results(deseqdata,name="condition_WT_vs_GFP"))
wt$Genes<-rownames(wt)
wt <- wt %>% mutate(gene_type = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "up",
                                            log2FoldChange <= -1 & padj < 0.05 ~ "down",
                                            TRUE ~ "ns"))
wt$Gene_Sym<-mapIds(org.Hs.eg.db,keys = wt$Genes,column = "SYMBOL",keytype="ENSEMBL",multiVals = "first")
write.table(wt,"WT_vs_GFP_Whole.txt",sep="\t",col.names = T,row.names = F,quote=F)

#L100Q
l100q<-data.frame(results(deseqdata,name="condition_L100Q_vs_GFP"))
l100q$Genes<-rownames(l100q)
l100q <- l100q %>% mutate(gene_type = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "up",
                                          log2FoldChange <= -1 & padj < 0.05 ~ "down",
                                          TRUE ~ "ns"))
l100q$Gene_Sym<-mapIds(org.Hs.eg.db,keys = l100q$Genes,column = "SYMBOL",keytype="ENSEMBL",multiVals = "first")
write.table(l100q,"L100Q_vs_GFP_Whole.txt",sep="\t",col.names = T,row.names = F,quote=F)





