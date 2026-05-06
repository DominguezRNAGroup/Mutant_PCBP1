library(tximport)
library(DESeq2)
library(org.Hs.eg.db)


rm(list=ls())
setwd(dir = "~/Desktop/Dominguez_Lab/PCBP1_FINAL/RNA-Seq/RSEM/HCT_new//")

#Loading the files
files <- list.files(pattern = "\\.genes.results$")
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
txi.rsem$length[txi.rsem$length == 0] <- 1

#Design
sample<-sapply(files,function(x){strsplit(x,"[.]")[[1]][1]})
condition<-factor(c("Parent","Parent","Parent",
                    "HDR1","HDR1","HDR1",
                    "HDR2","HDR2","HDR2",
                    "L100Q","L100Q","L100Q",
                    "F10","F10","F10",
                    "KO","KO","KO"),levels=c("Parent","HDR1","HDR2","L100Q","F10","KO"))
sampleInfo<-data.frame(sample,condition)
sampleInfo

#Reading into DESEQ2
deseqdata <- DESeqDataSetFromTximport(txi.rsem, colData = sampleInfo, design = ~condition)
deseqdata

# pre-filtering- Averaging at 5 per sample
keep <- rowSums(counts(deseqdata)) >= 90 #average of 5 per sample 
deseqdata <- deseqdata[keep,] #19432

#checking the factor levels
deseqdata$condition

#Running DESEQ
deseqdata<-DESeq(deseqdata)

#Looking at the variability between replicates
vsd <- vst(deseqdata)
plotPCA(vsd, "condition")

### Comparing fold change in mutant vs wildtype
resultsNames(deseqdata)

#HDR
res<-data.frame(results(deseqdata,name="condition_HDR1_vs_Parent"))
res$Genes<-rownames(res)
res <- res %>% mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                                            log2FoldChange <= -1 & padj <= 0.05 ~ "down",
                                            TRUE ~ "ns"))
res$Gene_Sym<-mapIds(org.Hs.eg.db,keys = res$Genes,column = "SYMBOL",keytype="ENSEMBL",multiVals = "first")
write.table(res,"HDR1_vs_WT_DESEQ.txt",sep="\t",col.names = T,row.names = F,quote=F)

#HDR2
res<-data.frame(results(deseqdata,name="condition_HDR2_vs_Parent"))
res$Genes<-rownames(res)
res <- res %>% mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                                            log2FoldChange <= -1 & padj <= 0.05 ~ "down",
                                            TRUE ~ "ns"))
res$Gene_Sym<-mapIds(org.Hs.eg.db,keys = res$Genes,column = "SYMBOL",keytype="ENSEMBL",multiVals = "first")
write.table(res,"HDR2_vs_WT_DESEQ.txt",sep="\t",col.names = T,row.names = F,quote=F)

#L100Q
res<-data.frame(results(deseqdata,name="condition_L100Q_vs_Parent"))
res$Genes<-rownames(res)
res <- res %>% mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                                            log2FoldChange <= -1 & padj <= 0.05 ~ "down",
                                            TRUE ~ "ns"))
res$Gene_Sym<-mapIds(org.Hs.eg.db,keys = res$Genes,column = "SYMBOL",keytype="ENSEMBL",multiVals = "first")
write.table(res,"L100Q_vs_WT_DESEQ.txt",sep="\t",col.names = T,row.names = F,quote=F)

#F10
res<-data.frame(results(deseqdata,name="condition_F10_vs_Parent"))
res$Genes<-rownames(res)
res <- res %>% mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                                            log2FoldChange <= -1 & padj <= 0.05 ~ "down",
                                            TRUE ~ "ns"))
res$Gene_Sym<-mapIds(org.Hs.eg.db,keys = res$Genes,column = "SYMBOL",keytype="ENSEMBL",multiVals = "first")
write.table(res,"F10_vs_WT_DESEQ.txt",sep="\t",col.names = T,row.names = F,quote=F)

#KO
res<-data.frame(results(deseqdata,name="condition_KO_vs_Parent"))
res$Genes<-rownames(res)
res <- res %>% mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                                            log2FoldChange <= -1 & padj <= 0.05 ~ "down",
                                            TRUE ~ "ns"))
res$Gene_Sym<-mapIds(org.Hs.eg.db,keys = res$Genes,column = "SYMBOL",keytype="ENSEMBL",multiVals = "first")
write.table(res,"KO_vs_WT_DESEQ.txt",sep="\t",col.names = T,row.names = F,quote=F)






