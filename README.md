
## Proteomics data anylysis 

This workflow performs the following tasks:
- [Spectronaut](https://biognosys.com/resources/spectronaut-the-deepest-proteome-coverage-available/) or [DIA-NN](https://github.com/vdemichev/DiaNN) for Protein Analysis and Quantification
- Data process and Cluster
- Differential expression analysis(DESeq2)
- Downstream functional enrichment analysis

## System requirements
- Windows for Spectronaut
- Linux for DIA-NN
- R 


## Spectronaut for Protein Analysis and Quantification
We uses the Spectronaut software for protein analysis and quantification in Windows systerm.

## DIA-NN for Protein Analysis and Quantification
We uses the DIA-NN software for protein analysis and quantification on biowulf.
``` bash
module load diann
diann --f ../20210208_KLOF_DIA_FAIMS_35V_d0_1.mzML   --lib  --threads 24 --verbose 1 --out ./report.tsv --qvalue 0.01 --matrices --out-lib ./report-lib.tsv --gen-spec-lib --predictor --fasta ../uniprot-proteome_Human_UP000005640_20191105.fasta --fasta-search --min-fr-mz 200 --max-fr-mz 2000 --met-excision --cut K*,R* --missed-cleavages 2 --min-pep-len 7 --max-pep-len 52 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4 --unimod4 --var-mods 5 --var-mod UniMod:35,15.994915,M --var-mod UniMod:1,42.010565,*n --monitor-mod UniMod:1 --reanalyse --relaxed-prot-inf --smart-profiling --peak-center --no-ifs-removal 
``` 

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

## Bash comand line
Downstream functional enrichment analysis for Differential expression gene.
```
bash pro.bash
```

