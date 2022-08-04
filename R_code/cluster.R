#R
#proteomics data cluster

library(ggplot2)
#read files_____________________________________________________________________________
setwd("/Users/liz36/Documents/INDI/proteomics/")
iNDI <- read.csv("iNDI_all_normarlized.csv")
iNDI_2nt <- read.csv("iNDI_2nt_normarlized.csv")
log2_iNDI <- read.csv("log2_iNDI_all_normarlized.csv")
iNDI_15group <- read.csv("iNDI_15group_normarlized.csv")
log2_iNDI_15group <- read.csv("log2_iNDI_15group_normarlized.csv")
#PCA and plot___________________________________________________________________________
#all sample
pca <- prcomp(t(log2_iNDI[,4:119]), center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
pca_df <- as.data.frame(pca$x)
pca_df$Condition <- gsub('_.*','',rownames(pca_df))
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500",
            "#9370DB","#98FB98","#F08080","#1E90FF",
            "#7CFC00","#FFFF00","#808000","#FF00FF",
            "#FA8072","#7B68EE","#9400D3","#800080",
            "#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD",
            "#008B8B","#FFE4B5","#8A2BE2","#228B22",
            "#E9967A","#4682B4","#32CD32","#F0E68C",
            "#FFFFE0","#EE82EE","#FF6347","#6A5ACD",
            "#9932CC","#8B008B","#8B4513","#DEB887")
pdf(file = 'pca_log2_iNDI.pdf',width = 8,height = 6)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size=4)+
  ggtitle("PCA plot") +
  scale_color_manual(values = allcolour)
dev.off()
#the samples we used
pca <- prcomp(t(log2_iNDI_15group[,4:ncol(log2_iNDI_15group)]), center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
pca_df <- as.data.frame(pca$x)
pca_df$Condition <- gsub('_.*','',rownames(pca_df))
pdf(file = 'pca_log2_iNDI_15group.pdf',width = 8,height = 6)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size=4)+
  ggtitle("PCA plot") 
dev.off()


#hc_cluster________________________________________________________________________________
dist_mat <- dist(t(log2_iNDI[,4:119])) #
hc_cluster <- hclust(dist_mat,method = "complete")
samplesname<- colnames(log2_iNDI[,4:119])

pdf(file = 'hc_cluster_all_log2_iND.pdf',width = 18,height =8 )
plot(hc_cluster,cex=0.8,col="dark red",labels = samplesname,main="HC Cluster")
dev.off()

#mean to one sample_______________________________________________________________________
#all samples
log2_iNDI_in1 <- data.frame(row.names = rownames(log2_iNDI))
for (i in unique(gsub("_.",'',colnames(log2_iNDI[,4:ncol(log2_iNDI)])))) {
  tmp <- log2_iNDI[,grep(i, colnames(log2_iNDI)),]
  tmp_i <- data.frame(rowMeans(tmp),row.names = rownames(log2_iNDI))
  colnames(tmp_i) <-i
  log2_iNDI_in1 <- cbind(log2_iNDI_in1,tmp_i)
}
dist_mat <- dist(t(log2_iNDI_in1)) #
hc_cluster <- hclust(dist_mat,method = "complete")
samplesname<- colnames(log2_iNDI_in1)
#hc cluster
pdf(file = 'hc_cluster_log2_iNDI_mean_all.pdf')
plot(hc_cluster,cex=0.8,col="dark red",labels = samplesname,main="HC Cluster")
dev.off()
#the sample we used
log2_iNDI_15group_in1 <- data.frame(row.names = rownames(log2_iNDI_15group))
for (i in unique(gsub("_.",'',colnames(log2_iNDI_15group[,4:ncol(log2_iNDI_15group)])))) {
  tmp <- log2_iNDI_15group[,grep(i, colnames(log2_iNDI_15group)),]
  tmp_i <- data.frame(rowMeans(tmp),row.names = rownames(log2_iNDI_15group))
  colnames(tmp_i) <-i
  log2_iNDI_15group_in1 <- cbind(log2_iNDI_15group_in1,tmp_i)
}
log2_iNDI_15group_in1$NT <- rowMeans(log2_iNDI_15group_in1[,grep("NT",colnames(log2_iNDI_15group_in1))])
log2_iNDI_15group_in1 <- log2_iNDI_15group_in1[,-grep("NT.",colnames(log2_iNDI_15group_in1))]
dist_mat <- dist(t(log2_iNDI_15group_in1)) #
hc_cluster <- hclust(dist_mat,method = "complete")
samplesname<- colnames(log2_iNDI_15group_in1)
#hc cluster
pdf(file = 'hc_cluster_log2_iNDI_15group_mean_all.pdf')
plot(hc_cluster,cex=0.8,col="dark red",labels = samplesname,main="HC Cluster")
dev.off()


##nt cluster_______________________________________________________________________________
log2_iNDI_NT <- log2_iNDI[,c(1:3,grep("NT",colnames(log2_iNDI)))]
dist_mat <- dist(t(log2_iNDI_NT[,4:21])) #
hc_cluster <- hclust(dist_mat,method = "complete")
samplesname<- colnames(log2_iNDI_NT[,4:21])
pdf(file = 'hc_cluster_nt.pdf')
plot(hc_cluster,cex=0.8,col="dark red",labels = samplesname,main="HC Cluster")
dev.off()
pca <- prcomp(t(log2_iNDI_NT[which(rowSums(log2_iNDI_NT[,4:21])>0),4:21]), center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
pca_df <- as.data.frame(pca$x)
pca_df$Condition <- gsub('_.*','',rownames(pca_df))
pdf(file = 'pca_log2_iNDI_nt.pdf')
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size=4)+
  ggtitle("PCA plot") +stat_ellipse()
dev.off()

#correlation______________________________________________________________________________
cor_matrix <- cor (as.matrix(log2_iNDI[,4:119]))
library(reshape2)
melted_cormat <- melt(cor_matrix)
pdf('cor_log2normailized_INDI.pdf',width = 15,height = 15)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 8, hjust = 1))
dev.off()

#umap_________________________________________________________________________________________
library(umap)
library(tidyverse)
umap_data <- t(log2_iNDI[,4:ncol(log2_iNDI)])
set.seed(100)
umap_out <- umap(umap_data)
umap_df <- umap_out$layout %>% 
  as.data.frame() %>%
  rename(umap1="V1",
         umap2="V2") %>%
  mutate(ID=row.names(umap_data))
umap_df$Condition <- gsub('_.*','',umap_df$ID)
umap_df$col <- 'KD_samples'
umap_df$col[grep("NT",umap_df$Condition)] <- "NT_samples"

pdf("umap_log2indi.pdf",width = 8,height = 6)
ggplot(umap_df) +
  geom_point(aes(x=umap1, y=umap2, color=Condition),size=4)+
  scale_color_manual(values = allcolour)
dev.off()

set.seed(100)
umap_data <- t(log2_iNDI_15group[,4:ncol(log2_iNDI_15group)])
umap_out <- umap(umap_data)
umap_df <- umap_out$layout %>% 
  as.data.frame()
colnames(umap_df) <-c("UMAP1","UMAP2")
umap_df$Condition <- gsub('_.*','',row.names(umap_data))
umap_df$col <- 'KD_samples'
umap_df$col[grep("NT",umap_df$Condition)] <- "NT_samples"

pdf("umap_log2_iNDI_15group.pdf",width = 8,height = 6)
ggplot(umap_df) +
  geom_point(aes(x=UMAP1, y=UMAP2, color=Condition),size=4)+
  scale_color_manual(values = allcolour)
dev.off()

