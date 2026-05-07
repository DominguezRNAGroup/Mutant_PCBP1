library(tximport)
library(DESeq2)
library(org.Hs.eg.db)
library(ggfortify)
library(dplyr)
library(ggrepel)
library(enrichplot)
library(clusterProfiler)

rm(list=ls())
setwd("RSEM/")

files <- list.files(pattern = "\\.genes.results$")
Counts<-list()

#Reading the files
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
txi.rsem$length[txi.rsem$length == 0] <- 1

##Design for DESeq2
sample<-sapply(files,function(x){strsplit(x,"[.]")[[1]][1]})
assay<-factor(sapply(sample,function(x){strsplit(x,"_|-")[[1]][3]}), levels=c("in","flagip"))
condition<-factor(sapply(sample,function(x){strsplit(x,"_|-")[[1]][2]}), levels=c("gfp","wt","mut"))
ind<-as.factor(c(rep(1,6),rep(2,6)))
sampleInfo<-data.frame(sample,condition,assay,ind) 
sampleInfo

#Reading into DESEQ2
deseqdata <- DESeqDataSetFromTximport(txi.rsem, colData = sampleInfo, design = ~condition+condition:ind+condition:assay)
deseqdata

#checking 
model.matrix(design(deseqdata),sampleInfo)

# pre-filtering- Averaging at 5 per sample
keep <- rowSums(counts(deseqdata)) >= 60 #18480 #Filter for more reads in RIP samples
deseqdata <- deseqdata[keep,]

#checking the factor levels
deseqdata$condition
deseqdata$assay

###Running DESEQ
deseqdata<-DESeq(deseqdata)

#looking at pca
vsd <- vst(deseqdata, blind=FALSE)
plotPCA(vsd, intgroup=c("assay"))
plotPCA(vsd, intgroup=c("ind"))
plotPCA(vsd, intgroup=c("condition","assay"))

#Looking at all the result outputs
resultsNames(deseqdata)

#WT vs GFP
res<-as.data.frame(results(deseqdata,contrast = list("conditionwt.assayflagip","conditiongfp.assayflagip")))
res$Genes<-rownames(res)
res <- res %>% mutate(gene_type = case_when(log2FoldChange >= (1) & padj < 0.05 ~ "up",
                                            log2FoldChange <= (-1) & padj < 0.05 ~ "down",
                                            TRUE ~ "ns"))
res$Gene_Sym<-mapIds(org.Hs.eg.db,keys = res$Genes,column = "SYMBOL",keytype="ENSEMBL",multiVals = "first")
res<-res[order(res$log2FoldChange,decreasing=T),]
write.table(res,"WT_vs_GFP_C2Bbe1.txt",sep="\t",col.names = T,row.names = F,quote=F)

#MUT vs GFP
res<-as.data.frame(results(deseqdata,contrast = list("conditionmut.assayflagip","conditiongfp.assayflagip")))
res$Genes<-rownames(res)
res <- res %>% mutate(gene_type = case_when(log2FoldChange >= (1) & padj < 0.05 ~ "up",
                                            log2FoldChange <= (-1) & padj < 0.05 ~ "down",
                                            TRUE ~ "ns"))
res$Gene_Sym<-mapIds(org.Hs.eg.db,keys = res$Genes,column = "SYMBOL",keytype="ENSEMBL",multiVals = "first")
res<-res[order(res$log2FoldChange,decreasing=T),]
write.table(res,"Mut_vs_GFP_C2Bbe1.txt",sep="\t",col.names = T,row.names = F,quote=F)
