# Protpipe Local downstream analysis

## General Overview
If you have the output file from the search enegine, you can run this downstream analysis without the sigularity installed.


## Installing

Opean your R and install the dependencies:

```R
bioc_deps <- c("clusterProfiler", "org.Hs.eg.db", "limma" )
cran_deps <- c("corrplot", "data.table", "ggplot2", "umap","ggbeeswarm","ggrepel",'ggdendro', 
               "pheatmap", "reshape2", "rlang", "magick",'ecodist',
               "ggthemes", "dplyr", "tidyr", "foreach") 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(bioc_deps) 
install.packages(cran_deps)
```

## Basic MS data analysis: `./protpipe.sh basic`

For performing QC and running differential abundance or enrichment analysis for typical mass spec data. The required inputs are
- protein intensity estimates from DIA-NN or Spectronaut
- experimental design matrix csv file
Open your ternimal and cd the dir of Protpipe:
```bash
Rscript src/basic_analysis.R \
    --pgfile EXAMPLES/DIFF_ABUNDANCE/iPSC.csv \
    --design EXAMPLES/DIFF_ABUNDANCE/design_matrix_iPSC.csv \
    --out EXAMPLES/DIFF_ABUNDANCE/
```

## Affinity Purification Mass Spec analysis: `./protpipe.sh APMS`

Similar to `basic` but for affinity purification mass spec. Requires the user to specify which protein was used for pulldown (`--ip`)
Open your ternimal and cd the dir of Protpipe:
```bash
Rscript src/APMS.R \
    --pgfile EXAMPLES/APMS/APMS.csv \
    --design EXAMPLES/APMS/design_matrix_APMS.csv \
    --ip UNC13A \
    --out EXAMPLES/APMS/
```

## Immunopeptidome analysis: `./protpipe.sh immuno`

Requires the csv or tsv output from `FragPipe` and a `csv` specifying HLA typing.
Open your ternimal and cd the dir of Protpipe:
```bash
Rscript src/immunopeptidome.R \
    --pepfile EXAMPLES/IMMUNOPEPTIDOME/combined_peptide.tsv \
    --out EXAMPLES/IMMUNOPEPTIDOME/ \
    --hla EXAMPLES/IMMUNOPEPTIDOME/HLA_typing.csv
```
