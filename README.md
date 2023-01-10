
## Proteomics data anylysis 

This workflow performs the following tasks:
- [Spectronaut](https://biognosys.com/resources/spectronaut-the-deepest-proteome-coverage-available/) or [DIA-NN](https://github.com/vdemichev/DiaNN) for Protein Analysis and Quantification
- Quality control and data filltering
- Normalization and imputation(OPTIONAL)
- Data Cluster
- Differential expression analysis and visualization
- Downstream functional enrichment analysis
- Run all steps from start to finish

## System requirements
- Windows for Spectronaut
- Linux for DIA-NN
- R 


## Spectronaut for Protein Analysis and Quantification
We uses the Spectronaut software for protein analysis and quantification in Windows systerm.

## DIA-NN for Protein Analysis and Quantification
We uses the DIA-NN software for protein analysis and quantification on biowulf or Windows systerm.

Use the following commands to run DIA-NN on biowulf：
``` bash
module load diann
diann --f ../20210208_KLOF_DIA_FAIMS_35V_d0_1.mzML   --lib  --threads 24 --verbose 1 --out ./report.tsv --qvalue 0.01 --matrices --out-lib ./report-lib.tsv --gen-spec-lib --predictor --fasta ../uniprot-proteome_Human_UP000005640_20191105.fasta --fasta-search --min-fr-mz 200 --max-fr-mz 2000 --met-excision --cut K*,R* --missed-cleavages 2 --min-pep-len 7 --max-pep-len 52 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4 --unimod4 --var-mods 5 --var-mod UniMod:35,15.994915,M --var-mod UniMod:1,42.010565,*n --monitor-mod UniMod:1 --reanalyse --relaxed-prot-inf --smart-profiling --peak-center --no-ifs-removal 
``` 

## Quality control and data filltering
We used the R to process the MS data, visualize the samples quality and fillter some sample with poor quality. This step will provide the figures including the total number of identified and quantified protein groups, the distribution of protein intensity, and correlation among protein abundance of the biological replicates.

Runing the code:
``` bash
Rscript QC.R --pro_input $pro_input --pep_input $pep_input -p $project_name -o $outdir
``` 

## Data Cluster
The data cluster include HC-cluster, PCA, and UMAP.

``` bash
Rscript cluster_plot.R -i $pro_input -p $project_name -o $outdir
``` 


## Differential expression analysis(T-test) and downstream functional enrichment analysis
DE analysis is done using the t-test. Downstream functional enrichment analysis for Differential expression gene. It takes the MS quatification from Spectronaut as an input:
```
Rscript DE_enrichment.R -i $pro_input -c $control -o $outdir 
```

## Bash comand line

```
bash spe_pro.bash
```

