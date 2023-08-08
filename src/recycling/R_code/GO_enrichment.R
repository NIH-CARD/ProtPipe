#R
#GO enrichment in iDNI proteomics project

## --- install and load  the package  manager--------------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#bio_pkgs = c("org.Hs.eg.db", "clusterProfiler", "ReactomePA")
#if (!requireNamespace(bio_pkgs, quietly = TRUE))
#  BiocManager::install(bio_pkgs)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
source('/Users/liz36/Documents/INDI/proteomics/go_enrichment/enrichment_clusterProfiler.R')
## --- parameters -------------------------------------------------------------
fdr_cutoff=0.05
lfc_cutoff=0.585
enrich_pvalue=0.05
gsea_fdr_cutoff=0.05


## --- input and output -------------------------------------------------------
out_dir <- "/Users/liz36/Documents/INDI/proteomics/go_enrichment/CT_desq2_fdr0.05_lgfc0.585_go0.05_gsea0.05/"
MSigDb_gmt_dir='/Users/liz36/Documents/INDI/proteomics/go_enrichment/entrez/'
input='/Users/liz36/Documents/INDI/proteomics/deseq2_results/'
setwd(out_dir)

## ------- preprocess data ---------------------------------------------
files=list.files(input)
for (i in files){
  if (!dir.exists(i)){
    dir.create(i,recursive = T)
  }
  
  deseq_res <- read.csv(paste0(input,i),row.names = 1)
  diff_genes <- deseq_res[,c((1:3),((ncol(deseq_res)-5):ncol(deseq_res)))]
  entrizid = data.frame(bitr(diff_genes$PG.Genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop=T))
  df_all=merge(entrizid,diff_genes,by.x='SYMBOL',by.y='PG.Genes')
  df_all=df_all[order(df_all$log2FoldChange,decreasing = T),]
  df_all=na.omit(df_all)
  all_gene_vector=df_all$log2FoldChange
  names(all_gene_vector)=df_all$ENTREZID
  ## up and down regulated genes 
  up_genes=df_all[which(df_all$log2FoldChange>=lfc_cutoff&df_all$padj<=fdr_cutoff),]
  down_genes=df_all[which(df_all$log2FoldChange<=(-lfc_cutoff)&df_all$padj<=fdr_cutoff),]
  ## ------- Enrichment -------
  
  if (nrow(up_genes)>0){
    enrichAll(gene_id=up_genes$ENTREZID, all_gene_vector=all_gene_vector, MSigDb_gmt_dir = MSigDb_gmt_dir, out_prefix = 'up',outdir=i,width = 12,height = 8,enrich_pvalue=enrich_pvalue)
  }
  
  if (nrow(down_genes)>0){
    enrichAll(gene_id=down_genes$ENTREZID, all_gene_vector=all_gene_vector, MSigDb_gmt_dir = MSigDb_gmt_dir, out_prefix = 'down',outdir=i,width = 12,height = 8,enrich_pvalue=enrich_pvalue)
  }
  ## ------- barplot of GSEA results -------
  gsea_res_path=paste0(out_dir,i)
  
  gsea_res=file.path(gsea_res_path,'/enrich_MSigDb_GSEA.txt')
  
  if (file.size(gsea_res) > 1L){
    gsea_res=read.delim(gsea_res)
    gsea_res_sig=gsea_res[which(gsea_res$p.adjust<=gsea_fdr_cutoff),]
    
    for (j in unique(gsea_res_sig$group)){
      tmp_gsea_res=gsea_res_sig[which(gsea_res_sig$group==j),]
      if (nrow(tmp_gsea_res)>=20){
        tmp_gsea_res_top=tmp_gsea_res[order(tmp_gsea_res$p.adjust,decreasing = F),]
        tmp_gsea_res_top=tmp_gsea_res_top[1:20,]
        
        
      }else{
        tmp_gsea_res_top=tmp_gsea_res
      }
      tmp_gsea_res_top=tmp_gsea_res_top[order(tmp_gsea_res_top$NES,decreasing = F),]
      tmp_gsea_res_top$Description=factor(tmp_gsea_res_top$Description,tmp_gsea_res_top$Description)
      
      # barplot
      p<-ggplot(data=tmp_gsea_res_top, aes(x=Description, y=NES,fill=NES)) +
        geom_bar(stat="identity")+scale_fill_gradient(low = "blue", high = "red")+
        coord_flip()+theme_bw()
      p
      ggsave(file.path(gsea_res_path,paste0('GSEA_',j,'.pdf')),width = 12,height = (2+nrow(tmp_gsea_res_top)*0.25))
    }
  }
}




