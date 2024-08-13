#########################
#phospho analysis
#July-2024
#Ziyi Li
#CARD NIH
#########################

#### ARG PARSING ###################################################################################
library(optparse)
pwd = getwd()
optparse_indent = '\n                '
option_list = list( 
  make_option(
    "--pepfile",
    default=NULL,
    help=paste(
      'Input file of Peptide Intensity (from DIA-NN or Spectronaut)',
      'Required.',
      sep=optparse_indent
    )
  ),
  make_option(
    "--PTM075",
    default=NULL,
    help=paste(
      'Input report uses the PTM Localization probability cutoff of 0.75 from Spectronaut',
      'Required.',
      sep=optparse_indent
    )
  ),
  make_option(
    "--PTM000",
    default=NULL,
    help=paste(
      'Input report uses the PTM Localization probability cutoff of 0 from Spectronaut',
      'Required.',
      sep=optparse_indent
    )
  ),
  make_option(
    "--out",
    dest="outdir",
    default='output', 
    help=paste(
      'Directory to direct all output. Directory will be created if does not exist).',
      'Defaults to the current working directory:',
      pwd,
      sep=optparse_indent
    )
  ),
  make_option(
    "--labelgene",
    dest="labelgene",
    default=NULL, 
    help='Gene to always label in output plots'
  ),
  make_option(
    "--heatmap",
    dest="heatmap",
    default=NULL, 
    help='Gene to plot heatmap'
  ),
  make_option(
    "--base",
    dest="log_base",
    default=10, 
    help='Base for log transformation of intensity data. Default: 10'
  ),
  make_option(
    "--normalize",
    default='none',
    type='character',
    help=paste(
      'shift: adjust sample intensities to match global median by adding a constant',
      'scale: adjust sample intensities to match global median by multiplicative scaling',
      'none: do not normalize',
      sep=optparse_indent
    )
  ),
  make_option(
    "--exclude",
    default=NULL,
    type='character',
    help=paste(
      'semicolon-separated string of files to exclude from analysis'
    )
  ),
  make_option(
    "--sds",
    dest = 'sds',
    default=3,
    type='numeric',
    help=paste(
      'Filter out samples with protein group counts > N standard deviations from the mean.',
      'Increase to higher values for greater tolerance of variance in protein group counts.',
      'Default: 3',
      sep=optparse_indent
    )
  ),
  make_option(
    "--minintensity",
    dest = 'minintensity',
    default=500,
    type='numeric',
    help='Minimum LINEAR (not log) intensity. Default: 500'
  ),
  make_option(
    "--fdr",
    dest = 'fdr_threshold',
    default=0.01,
    type='numeric',
    help=paste(
      'False Discovery Rate threshold for differential abundance analysis.',
      'Default: 0.01',
      sep=optparse_indent
    )
  ),
  make_option(
    "--lfc_threshold",
    dest = 'lfc_threshold',
    default=1,
    type='numeric',
    help=paste(
      'Minimum LINEAR fold change [NOT log, as log base can be modified] for labeling',
      'protein groups in differential abundance analysis. Default: 2 (equivalent to',
      'log2 fold-change threshold of 1)',
      sep=optparse_indent
    )
  ),
  make_option(
    "--imputation",
    action = 'store_true',
    default=FALSE, 
    type='logical',
    help=paste(
      'Applies data imputation. Not yet implimented.',
      sep=optparse_indent
    )
  ),
  make_option(
    "--design",
    default=NULL,  
    help=paste(
      'Comma- or tab-delimited, three-column text file specifying the experimental design.',
      'File should contain headers. Header names do not matter; column order DOES matter.',
      'Columns order: <sample_name> <condition> <control>',
      sep=optparse_indent
    )
  ),
  make_option(
    "--DE_method",
    dest="DE_method",
    default='ttest',
    type='character',
    help=paste(
      'ttest',
      'limma',
      sep=optparse_indent
    )
  ),
  make_option(
    "--neighbors",
    default=15,
    type='numeric',
    help=paste(
      'N Neighbors to use for UMAP. Default: 15',
      sep=optparse_indent
    )
  ),
  make_option(
    "--dry",
    action = 'store_true',
    default=FALSE, 
    type='logical',
    help=paste(
      'Applies data imputation. Not yet implimented.',
      sep=optparse_indent
    )
  ),
  make_option(
    "--enrich",
    dest = 'enrich_pvalue',
    default=0.01,
    type='numeric',
    help=paste(
      'The cutoff of p-value for gene enrichment analysis.',
      'Default: 0.01',
      sep=optparse_indent
    )
  ),
  make_option(
    "--gsea",
    dest = 'gsea_fdr_cutoff',
    default=0.01,
    type='numeric',
    help=paste(
      'The cutoff False Discovery Rate of gsea analysis.',
      'Default: 0.01',
      sep=optparse_indent
    )
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$PTM075)  && is.null(opt$PTM000) && is.null(opt$pepfile)) {
  cat("ERROR: --pepfile <file> or --PTM075 <file> and --PTM000<file> must be provided\n")
  badargs <- TRUE
}

#### Source the function ###
source('src/functions.R')
#### PACKAGES ######################################################################################
package_list = c('plyr','dplyr','ggplot2','ggridges', 'data.table', 'corrplot', 
                 'umap', 'magick', 'ggdendro', 'ecodist','ggbeeswarm', 
                 'ggrepel', 'ggthemes', 'foreach','reshape2','org.Hs.eg.db',
                 'clusterProfiler','pheatmap','limma','UpSetR','maSigPro','KSEAapp','DOSE')
cat("INFO: Loading required packages\n      ")
cat(paste(package_list, collapse='\n      ')); cat('\n')

defaultW <- getOption("warn"); options(warn = -1)   # Temporarily disable warnings for quiet loading
if(all((lapply(package_list, require, character.only=TRUE)))) {
  cat("INFO: All packages successfully loaded\n")
} else {
  cat("ERROR: One or more packages not available. Are you running this within the container?\n")
}
options(warn = defaultW)    # Turn warnings back on


#### IMPORT AND FORMAT DATA#########################################################################
##peptide file
if (!is.null(opt$pepfile)) {
  #### IMPORT AND FORMAT DATA#########################################################################
  tryTo(paste0('INFO: Reading input file ', opt$pepfile),{
    dat <- fread(opt$pepfile)
  }, paste0('ERROR: problem trying to load ', opt$pepfile, ', does it exist?'))
  
  tryTo(paste0('INFO: Massaging data from ', opt$pepfile, ' into a common style format for processing'), {
    dat <- phospho_standardize_format(dat)
  }, 'ERROR: failed! Check for missing/corrupt headers?')
}

##PTM000 and PTM075
if (!is.null(opt$PTM075)  && !is.null(opt$PTM000)) {
  #### IMPORT AND FORMAT DATA#########################################################################
  tryTo(paste0('INFO: Reading input file ', opt$PTM075,' and ',opt$PTM000),{
    data_ptm075=fread(opt$PTM075)%>%
      phospho_standardize_format()
    data_ptm000=fread(opt$PTM000)%>%
      phospho_standardize_format()
  }, paste0('ERROR: problem trying to load ', opt$pepfile, ', does it exist?'))
  
  tryTo(paste0('INFO: Fillter the “PTM000” report to only contain precursors (EG.PrecursorId) exported in the “PTM075” '), {
    dat=data_ptm000 %>%
      filter(Precursor %in% data_ptm075$Precursor)
  }, 'ERROR: failed! Check for missing/corrupt headers?')
  
}else{
  cat("ERROR: Both --PTM075 <file> and --PTM000<file> must be provided\n")
  badargs <- TRUE
}

#### DATA FILTER###############################################
tryTo(paste0('INFO: Applying Filter Intensity > ',opt$minintensity),{
  dat <- dat%>%
    mutate(across(where(is.numeric), ~ ifelse(. < opt$minintensity, NA, .)))
}, 'ERROR: failed!')

tryTo(paste0('INFO: Identify unique PTM with the least missing values and highest median intensity'), {
  phospho_dat <- dat %>%
    #Filter for phospho sites
    filter(grepl('Phospho', Modified_Sequence))%>%
    # Generate PTM identifiers
    mutate(PTM = gsub(',*C[0-9]+,*|,*M[0-9]+,*', '', paste0(Genes, PTM_Location)))%>%
    # Calculate missing values 
    mutate(missing_value = rowSums(is.na(select(., -contains("Protein_Group|Genes|PTM|Precursor|Modified_Sequence")))))%>% 
    # Calculate median values
    mutate(median = rowMedians(as.matrix(select(., -contains("Protein_Group|Genes|PTM|Precursor|Modified_Sequence|missing_value")) %>% 
                                           select_if(is.numeric)), na.rm = TRUE)) %>%
    #Identify unique PTM with the least missing values and highest median intensity
    group_by(PTM) %>%
    filter(missing_value == min(missing_value)) %>%
    slice(which.max(median)) %>%
    ungroup() %>%
    select(-c(missing_value, median)) 
  ezwrite(phospho_dat,output_dir = opt$outdir,output_filename ='phospho_data.tsv' )
}, 'ERROR: failed! Check for colnames?')

#Converting to long format
tryTo(paste0('INFO: Converting to long format'), {
  dat_long <- melt_intensity_table(dat)
  phospho_dat_long=melt_intensity_table(phospho_dat)
}, 'ERROR: failed! Check for missing/corrupt headers?')

tryTo('INFO: Excluding all unquantified or zero intensities', {
  dat_long <- dat_long[! is.na(Intensity)][Intensity != 0]
  phospho_dat_long <- phospho_dat_long[! is.na(Intensity)][Intensity != 0]
}, 'ERROR: failed!')

#### QC ############################################################################################
## Make QC dir
QC_dir <- paste0(opt$outdir, '/QC/')
if(! dir.exists(QC_dir)){
  dir.create(QC_dir, recursive = T)
}

## Plotting intensity distribution
tryTo('INFO: Plotting Precursor counts',{
  plot_phospho_site_counts(DT_phospho_long =phospho_dat_long,
                           output_dir = QC_dir,
                           height = 6, 
                           width = 12)
   }, 'ERROR: failed!')
tryTo('INFO: Plotting intensity distribution',{
  plot_phospho_site_intensity(DT_phospho_long =phospho_dat_long,
                              output_dir = QC_dir)
}, 'ERROR: failed!')

#### CLUSTERING ####################################################################################
cluster_dir <- paste0(opt$outdir, '/Clustering/')
if(! dir.exists(cluster_dir)){
  dir.create(cluster_dir, recursive = T)
}
# PCA
tryTo('INFO: running PCA and plotting first two components',{
  pca <- get_PCs(phospho_dat)
  ezwrite(pca$components, cluster_dir, 'PCA.tsv')
  ezwrite(pca$summary, cluster_dir, 'PCA_summary.tsv')
  plot_PCs(pca, cluster_dir, 'PCA.pdf')
}, 'ERROR: failed!')

# Hierarchical Clustering
tryTo('INFO: running Hierarchical Clustering',{
  plot_hierarchical_cluster(phospho_dat, cluster_dir)
}, 'ERROR: failed!')

# UMAP
if ((ncol(dat)-3)>opt$neighbors) {
  tryTo('INFO: running UMAP',{
    umap <- get_umap(dat, opt$neighbors)
    ezwrite(umap, cluster_dir, 'UMAP.tsv')
    plot_umap(umap, cluster_dir, 'UMAP.pdf')
  }, 'ERROR: failed!')
}



#### DE and pathway analysis##########################
if (!is.null(opt$design)) {
  tryTo('INFO: Importing experimental design',{
    design <- fread(opt$design, header=TRUE)
    setnames(design, c('sample_name', 'condition', 'control'))
  }, 'ERROR: failed!')
  tryTo('INFO: Validating experimental design',{
    print(design[])
    cat('\n')
    conditions <- unique(design$condition)
    for (condition.i in conditions) {
      samples <- design[condition == condition.i, sample_name]
      control <- unique(design[condition == condition.i, control])
      if (length(control) != 1) {
        cat(paste0('ERROR: condition ', condition.i, ' maps to multiple controls: ', control, '\n'))
        cat(paste0('       Check the design matrix and esure no more than one control label per condition\n'))
        quit(exit=1)
      } else {
        cat(paste0('INFO: condition ', condition.i, ' maps to control ', control, '\n'))
      }
    }
    cat(paste0('INFO: all conditions pass check (i.e. map to one control condition)\n'))
  }, 'ERROR: failed!')
  ##DE ttest
  if (opt$DE_method == 'ttest') {
    DE_dir <- paste0(opt$outdir, '/ttest/Differential_Intensity/')
    if(! dir.exists(DE_dir)){
      dir.create(DE_dir, recursive = T)
    }
    EA_dir <- paste0(opt$outdir, '/ttest/Enrichiment_Analysis/')
    if(! dir.exists(EA_dir)){
      dir.create(EA_dir, recursive = T)
    }
    tryTo('INFO: Running differential intensity t-tests and pathway analysis',{
      phospho_ttest(DT = phospho_dat,
                    design_matrix = design,
                    DE_dir =DE_dir ,
                    EA_dir = EA_dir)
      
    }, 'ERROR: DE ttest failed!')
  }
  ##DE limma
  else if (opt$DE_method == 'limma') {
    DE_dir <- paste0(opt$outdir, '/limma/Differential_Intensity/')
    if(! dir.exists(DE_dir)){
      dir.create(DE_dir, recursive = T)
      }
    EA_dir <- paste0(opt$outdir, '/limma/Enrichiment_Analysis/')
    if(! dir.exists(EA_dir)){
      dir.create(EA_dir, recursive = T)
      }
    Log2_phospho_dat=as.data.frame(phospho_dat) %>%
      mutate_if(is.numeric, ~ log2(. ))
    phospho_limma(Log2_DT = Log2_phospho_dat,
             design_matrix = design,
             DE_dir = DE_dir,
             EA_dir=EA_dir)
  }
} 
#### HEATMAP ########################################################################
plot_heatmap(dat)
if (!is.null(opt$heatmap)) {
  plot_heatmap_subset(dat,opt$heatmap)
}

quit()

