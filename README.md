
## Proteomics data anylysis 

This workflow performs the following tasks:
- Spectronaut for Protein Analysis and Quantification
- Data process and Cluster
- Differential expression analysis(DESeq2)
- Downstream functional enrichment analysis

## System requirements
- Windows for Spectronaut
- Python
- R 


## Installation
We uses the Miniconda3 package management system to harmonize all of the software packages. 
Use the following commands to install Minicoda3：
``` bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
### Create an isolated environment for CRISPR screen analysis
``` bash
conda create -n CRSIPR
conda activate CRSIPR
``` 

### Install tools
Tools needed for this analysis are: R, MAGeCK, libxml2, MAGeCK-VISPR,MAGeCKFlute 
~~~
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels r
conda install -c anaconda libxml2
conda install -c bioconda mageck
conda install -c bioconda -c conda-forge mageck-vispr
~~~

### Install MAGeCKFlute using R
~~~
R

> install.packages(c("devtools", "BiocManager"), repos = "https://cloud.r-project.org")
> BiocManager::install(c("pathview", "biomaRt", "msigdbr", "dendextend", "pheatmap", "sva", "ggrepel", "knitr", "clusterProfiler", "depmap"))
> BiocManager::install("MAGeCKFlute") # Released version
# Or
> devtools::install_github("liulab-dfci/MAGeCKFlute")
~~~

## Processing of CRISPR screen data with MAGeCK or MAGeCK-VISPR
In a typical use case, CRISPR screen data are processed with MAGeCK (option A) step by step. If users want to perform QC and visualize the results, we recommend MAGeCK-VISPR (option B) instead. 

### A.Process CRISPR screen data step by step with MAGeCK ● Timing 1.5 h
```
#Activate the CRSIPR environment 
conda activate CRSIPR

# Download and unzip the test data for both datasets, using the following commands:
wget http://cistrome.org/MAGeCKFlute/demo.tar.gz 
tar zxvf demo.tar.gz
cd demo_data

# Generate a count table for Dataset 1 with the mageck count function, by first changing the working directory to a directory that contains raw .fastq data and is able to store the output of mageck count as follows:
cd path/to/demo_data/mageck_count

# To run the mageck count on Dataset 1, type the following command:
mageck count -l library.csv -n GSC_0131 --sample-label day0_r1, day0_r2,day23_r1,day23_r2 --fastq GSC_0131_Day0_Rep1.fastq.gz GSC_0131_Day0_Rep2.fastq.gz GSC_0131_Day23_Rep1.fastq.gz GSC_0131_ Day23_Rep2.fastq.gz

# Identify screen hits using MAGeCK RRA(MAGeCK RRA for comparison between two conditions, such as an initial condition versus cells cultured for a period of time. )
$mageck test -k GSC_0131.count.txt -t day23_r1,day23_r2 -c day0_r1,day0_r2 -n GSC_0131_rra --remove-zero both --remove- zero-threshold 0

# Identify screen hits using MAGeCK MLE. (If an experiment contains more than two conditions, for example, a three-condition design: day 0, drug treatment and DMSO treatment, we recommend using MAGeCK MLE)
cd path/to/demo_data/mageck_mle
mageck mle --count-table rawcount.txt --design-matrix designmatrix. txt --norm-method control --control-sgrna nonessential_ctrl_sgrna_ list.txt --output-prefix braf.mle ##modify designmatrix table as your expriemnt design
```

### (B) Process CRISPR screen data with MAGeCK-VISPR ● Timing 1.5 h
```
#Activate the CRSIPR environment 
conda activate CRSIPR

#Choose a workflow directory and initialize the workflow with the .fastq or .fastq.gz files that contain the raw reads
mageck-vispr init workflow --reads path/to/file/*.fastq*

#Configure the workflow
cd workflow

# To check whether the ‘config.yaml’ files have been configured correctly,enter the following command line into the terminal:
$snakemake –n

# Execute the workflow.
snakemake --cores 8

# Visualize the results with VISPR.
vispr server results/*.vispr.yaml
```

## MAGeCKFlute 
This package implements methods to perform quality control (QC), normalization, batch effect removal, gene hit identification and downstream functional enrichment analysis for CRISPR screens. 

### Functional analysis for MAGeCK RRA results
```
#Activate the CRSIPR environment 
conda activate CRSIPR

#use the R interactive shell
R

 >library(MAGeCKFlute)
 >setwd(‘path/to/file/’)
 
# To perform functional analysis for MAGeCK RRA results 
>FluteRRA(gene_summary = "path/to/file/rra.gene_summary.txt", prefix="FluteRRA", organism="hsa")
 
```

### Functional analysis for MAGeCK RRA results
```
#Activate the CRSIPR environment 
conda activate CRSIPR

#use the R interactive shell
R

 >library(MAGeCKFlute)
 >setwd(‘path/to/file/’)
 
# To perform functional analysis for MAGeCK MLE  results 
>FluteMLE(gene_summary="path/to/file/mle.gene_summary.txt", ctrlname="dmso", treatname="plx", organism="hsa", prefix="FluteMLE", -pathway_limit = c(3,50))


```

### (Optional) Batch effect removal
```
R
> library(MAGeCKFlute)
> BatchRemove(mat = "rawcount.txt", batchMat = "BatchMatrix. txt", prefix = "BatchCorrect", -pca = T, -cluster = T, -outdir = ".")
```

### (Optional) Correct copy-number bias. 
MAGeCK RRA and MAGeCK MLE contain an optional method to correct copy-number biases in the calculated RRA scores and beta scores, respectively. We recommend that users perform copy-number bias correction if the CNV information is available for the cell line. 

```
# To perform MAGeCK RRA with copy-number bias correction, type the following:
mageck test -k rawcount.txt -t HL60.final -c HL60.initial -n rra_ cnv --cnv-norm cnv_data.txt –cell-line HL60_HAEMATOPOIETIC_AND_ LYMPHOID_TISSUE

# To perform MAGeCK MLE with copy-number bias correction, type the following:
mageck mle --count-table rawcount.txt --design-matrix designma- trix.txt --cnv-norm cnv_data.txt

```
