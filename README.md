
## Proteomics data anylysis 

This workflow performs the following tasks:
- [Spectronaut](https://biognosys.com/resources/spectronaut-the-deepest-proteome-coverage-available/) or [DIA-NN](https://github.com/vdemichev/DiaNN) for Protein Analysis and Quantification
- Data process and Cluster
- Differential expression analysis(DESeq2)
- Downstream functional enrichment analysis

## System requirements
- Windows for Spectronaut
- R 


## Spectronaut for Protein Analysis and Quantification
We uses the Spectronaut software for protein analysis and quantification in Windows systerm.

## Data process and Cluster
We used the R to process the MS data and visualize the samples cluster(including HC-cluster, PCA, UMAP,Correlation)
``` bash
Rscript data.R
Rscript cluster.R
``` 

## Differential expression analysis(DESeq2)
DE analysis is done using the DESeq2 Bioconductor package. It takes the merged raw read counts (from Spectronaut) as an input:
```
Rscript ~rcode/deseq2.R
```

## Downstream functional enrichment analysis
Downstream functional enrichment analysis for Differential expression gene.
```
Rscript ~rcode/GO_enrichment.R
```
