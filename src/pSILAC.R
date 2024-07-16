#########################
#pSALIC analysis
#June-2024
#Ziyi Li
#CARD NIH
#########################

#### ARG PARSING ###################################################################################
library(optparse)
pwd = getwd()
optparse_indent = '\n                '
option_list = list( 
  make_option(
    "--pgfile",
    default=NULL,
    help=paste(
      'Input file of Protein Group Intensity (from DIA-NN or Spectronaut)',
      'Required.',
      sep=optparse_indent
    )
  ),
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
    default=0,
    type='numeric',
    help='Minimum LINEAR (not log) intensity. Default: 0'
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
    "--foldchange",
    dest = 'foldchange',
    default=2,
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



if (is.null(opt$pgfile) && is.null(opt$pepfile)) {
  cat("ERROR: --pgfile <file> or --pepfile  <file> must be provided\n")
  badargs <- TRUE
}

#### Source the function ###
source('src/functions.R')
#### PACKAGES ######################################################################################
package_list = c('ggplot2', 'data.table', 'corrplot', 'umap', 
                 'magick', 'ggdendro', 'ecodist','ggbeeswarm',
                 'ggrepel', 'ggthemes', 'foreach','reshape2',
                 'org.Hs.eg.db','clusterProfiler','pheatmap',
                 'limma','ggridges')
cat("INFO: Loading required packages\n      ")
cat(paste(package_list, collapse='\n      ')); cat('\n')

defaultW <- getOption("warn"); options(warn = -1)   # Temporarily disable warnings for quiet loading
if(all((lapply(package_list, require, character.only=TRUE)))) {
  cat("INFO: All packages successfully loaded\n")
} else {
  cat("ERROR: One or more packages not available. Are you running this within the container?\n")
}
options(warn = defaultW)    # Turn warnings back on

##pepfile########  
if (!is.null(opt$pepfile)) {
  #### IMPORT AND FORMAT DATA#########################################################################
  tryTo(paste0('INFO: Reading input file ', opt$pepfile),{
    dat <- fread(opt$pepfile)
  }, paste0('ERROR: problem trying to load ', opt$pepfile, ', does it exist?'))
  
  tryTo(paste0('INFO: Massaging data from ', opt$pepfile, ' into a common style format for processing'), {
    dat <- psilac_standardize_format(dat)
  }, 'ERROR: failed! Check for missing/corrupt headers?')
  
  #Converting to long format
  tryTo(paste0('INFO: Converting to long format'), {
    dat.long <- melt_intensity_table(dat)
  }, 'ERROR: failed! Check for missing/corrupt headers?')
  
  tryTo('INFO: Excluding all unquantified or zero intensities', {
    dat.long <- dat.long[! is.na(Intensity)][Intensity != 0]
  }, 'ERROR: failed!')
  
  tryTo(paste0('INFO: Applying Filter Intensity > ',opt$minintensity),{
    dat.long <- dat.long[Intensity > opt$minintensity]
  }, 'ERROR: failed!')
  
  #### QC ############################################################################################
  ## Make QC dir
  QC_dir <- paste0(opt$outdir, '/QC/')
  if(! dir.exists(QC_dir)){
    dir.create(QC_dir, recursive = T)
  }
  ## Plotting intensity distribution
  
  tryTo('INFO: Plotting counts',{
    plot_silac_pep_counts(DT.long = dat.long,
                          output_dir = QC_dir,
                          height = 6, 
                          width = 12)
    }, 'ERROR: failed!')
  
  tryTo('INFO: Plotting intensity distribution',{
    plot_silac_pep_intensity(DT.long = dat.long,
                             output_dir = QC_dir,
                             height = 6, 
                             width = 12)
    }, 'ERROR: failed!')
}

##pgfile########  
if (!is.null(opt$pgfile)) {
  #### IMPORT AND FORMAT DATA#########################################################################
  tryTo(paste0('INFO: Reading input file ', opt$pgfile),{
    dat <- fread(opt$pgfile)
  }, paste0('ERROR: problem trying to load ', opt$pgfile, ', does it exist?'))
  
  tryTo(paste0('INFO: Massaging data from ', opt$pgfile, ' into a common style format for processing'), {
    dat <- psilac_standardize_format(dat)
  }, 'ERROR: failed! Check for missing/corrupt headers?')
  
  #Converting to long format
  tryTo(paste0('INFO: Converting to long format'), {
    dat.long <- melt_intensity_table(dat)
  }, 'ERROR: failed! Check for missing/corrupt headers?')
  
  tryTo('INFO: Excluding all unquantified or zero intensities', {
    dat.long <- dat.long[! is.na(Intensity)][Intensity != 0]
  }, 'ERROR: failed!')
  
  tryTo(paste0('INFO: Applying Filter Intensity > ',opt$minintensity),{
    dat.long <- dat.long[Intensity > opt$minintensity]
  }, 'ERROR: failed!')
  
  #### QC ############################################################################################
  ## Make QC dir
  QC_dir <- paste0(opt$outdir, '/QC/')
  if(! dir.exists(QC_dir)){
    dir.create(QC_dir, recursive = T)
  }
  ## Plotting intensity distribution
  
  tryTo('INFO: Plotting counts',{
    plot_silac_pg_counts(DT.long = dat.long,
                          output_dir = QC_dir,
                          height = 6, 
                          width = 12)
  }, 'ERROR: failed!')
  
  tryTo('INFO: Plotting intensity distribution',{
    plot_silac_pg_intensity(DT.long = dat.long,
                             output_dir = QC_dir,
                             height = 6, 
                             width = 12)
  }, 'ERROR: failed!')
}


quit()

