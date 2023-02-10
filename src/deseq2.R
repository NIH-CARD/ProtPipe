#R
#deseq2 in iDNI proteomics project__________________________________________________
library(DESeq2)
library(ggplot2)
library(ggrepel)
#read files_________________________________________________________________________
setwd("/Users/liz36/Documents/INDI/proteomics/")
pro_nor_counts <- read.csv("iNDI_15group_normarlized.csv")

#deseq2_____________________________________________________________________________
samplename <-unique((gsub("_.","",colnames(pro_nor_counts)[4:ncol(pro_nor_counts)])))
samplename <- samplename[-grep("NT",samplename)]

for (i in samplename) {
  tmp1 <- pro_nor_counts[,c(grep(i,colnames(pro_nor_counts)))]
  count <- cbind(tmp1,pro_nor_counts[,c(grep("NT",colnames(pro_nor_counts)))])
  rownames(count) <- pro_nor_counts$PG.ProteinGroups
  sampleTable<- data.frame(sampleName = colnames(count),
                           sampleCondition = c(rep('treat', times=ncol(tmp1)),rep('control', times=12)))
  dds <- DESeqDataSetFromMatrix(countData = count,
                                colData=sampleTable,
                                design=~sampleCondition)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c('sampleCondition', 'treat', 'control'))
  dds_df  <- as.data.frame(subset(res))
  tmp2<- merge(count,dds_df,by=0,all=F)
  tmp3 <-merge(pro_nor_counts[,1:3],tmp2,by=1,all=F)
  filename <- paste0("deseq2_result/",i,"_NT_dds.csv")
  write.csv(tmp3,file = paste0("deseq2_results/",i,"_NT_dds.csv",row.names = F))
  }

