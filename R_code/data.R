#R
#proteomics data modified

#read files------------------------------------------------------------------------------
setwd("/Users/liz36/Documents/INDI/proteomics/")
iNDI <- read.csv("20220630_130647_20220629_mini_iNDI_normarlized.csv")

#modify-----------------------------------------------------------------------------------
#integer & NA
iNDI[,4:ncol(iNDI)] <- as.data.frame(apply(iNDI[,4:ncol(iNDI)],2,as.integer))
iNDI[is.na(iNDI)] <- 0
#change name
iNDI <- iNDI[,-grep("T1A1_4",colnames(iNDI))] 
colnames(iNDI) <- gsub('X.*35V_','',colnames(iNDI))
colnames(iNDI) <- gsub('.raw.*','',colnames(iNDI))
colnames(iNDI) <- gsub('TARDBP_14','TARDBP_4',colnames(iNDI))
colnames(iNDI) <- gsub('HNRN2B','HNRNPA2B1',colnames(iNDI))
colnames(iNDI) <- gsub('EWSR','EWSR1',colnames(iNDI))
colnames(iNDI) <- gsub('MART','MATR3',colnames(iNDI))
colnames(iNDI) <- gsub('T1A1','TIA1',colnames(iNDI))
colnames(iNDI) <- gsub('TBR','TBK1',colnames(iNDI))
colnames(iNDI) <- gsub('TMEM','TMEM106B',colnames(iNDI))
#order
iNDI <- iNDI[,c(colnames(iNDI)[1:3],sort(colnames(iNDI[,4:ncol(iNDI)])))]
iNDI[,4:ncol(iNDI)] <- as.data.frame(apply(iNDI[,4:ncol(iNDI)],2,as.integer))
write.csv(iNDI,"iNDI_all_normarlized.csv",row.names = F)
#log2
log2_iNDI <- log2(iNDI[,4:ncol(iNDI)]+1)
log2_iNDI <- cbind(iNDI[1:3],log2_iNDI)
write.csv(log2_iNDI,"log2_iNDI_all_normarlized.csv",row.names = F)
#NT 2 sample
iNDI_2nt <- iNDI[,-grep('NT76', colnames(iNDI))]
write.csv(iNDI_2nt,"iNDI_2nt_normarlized.csv",row.names = F)
log2_iNDI_nt <- log2(iNDI_2nt[,4:ncol(iNDI_2nt)]+1)
log2_iNDI_nt <- cbind(iNDI_2nt[1:3],log2_iNDI_nt)
write.csv(log2_iNDI_nt,"log2_iNDI_2nt_normarlized.csv",row.names = F)

#KD efficiency--------------------------------------------------------------------------
library(reshape)
KO_sheet <- data.frame()
KOsample <- c(unique(gsub("_.",'',colnames(iNDI_2nt[,4:ncol(iNDI_2nt)]))))
KOsample <- KOsample[-(grep("NT.*",KOsample))]
for (i in KOsample) {
  tmp1 <- iNDI_2nt[which(iNDI_2nt$PG.Genes == i),grep(i, colnames(iNDI_2nt)),]
  tmp2 <- iNDI_2nt[which(iNDI_2nt$PG.Genes == i),grep('NT', colnames(iNDI_2nt)),]
  tmp_i <- cbind(tmp1,tmp2)
  tmp3<-melt(tmp_i)
  tmp3$gene <- i
  KO_sheet <- rbind(KO_sheet,tmp3)
}
library(ggplot2)
library(ggpubr)#compare_means(): easy to use solution to performs one and multiple mean comparisons.
library(scales)# Log2 scaling of the y axis (with visually-equal spacing)
KO_sheet$sample <- 'KO'
KO_sheet$sample[grep('NT24',KO_sheet$variable)] <- "NT24"
KO_sheet$sample[grep('NT94',KO_sheet$variable)] <- "NT94"
KO_sheet$number <- KO_sheet$value
KO_sheet$number[which(KO_sheet$number==0)]  <- 1 
# Global test
compare_means(number ~ gene,  data = KO_sheet, method = "anova")

write.csv(KO_sheet,"KO_sheet.csv",row.names = F)
pdf("KO_efficiency.pdf",width = 16,height = 8)
ggplot(KO_sheet, aes(x=gene, y=number, fill=sample)) +
  geom_boxplot()+
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))+
  labs(y = "Normalized Quantification of Protein Level",x= "Gene")
dev.off()

pdf("KO_efficiency.pdf",width = 14,height = 8)

ggplot(KO_sheet, aes(x=gene, y=number, fill=sample)) +
  geom_boxplot()+
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  labs(y = "Normalized Quantification of Protein Level",x= "Gene")
dev.off()

#sample down___________________________________________________________________________________
#'CHMP2B',"OPTN","SETX","VCP"    
iNDI_15group <- iNDI_2nt[,-grep('CHMP2B|OPTN|SETX|VCP', colnames(iNDI_2nt))]
iNDI_15group <- iNDI_15group[-which(rowMeans(iNDI_15group[,4:ncol(iNDI_15group)])==0),]
write.csv(iNDI_15group,"iNDI_15group_normarlized.csv",row.names = F)
log2_iNDI_15group <- log2(iNDI_15group[,4:ncol(iNDI_15group)]+1)
log2_iNDI_15group <- cbind(iNDI_15group[,1:3],log2_iNDI_15group)
write.csv(log2_iNDI_15group,"log2_iNDI_15group_normarlized.csv",row.names = F)




