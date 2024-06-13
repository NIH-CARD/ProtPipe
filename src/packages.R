bioc_deps <- c("clusterProfiler", "org.Hs.eg.db", "limma" )
cran_deps <- c("corrplot", "data.table", "ggplot2", "umap","ggbeeswarm","ggrepel",'ggdendro', 
               "pheatmap", "reshape2", "rlang", "magick",'ecodist',
               "ggthemes", "dplyr", "tidyr", "foreach") 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(bioc_deps) 
install.packages(cran_deps)