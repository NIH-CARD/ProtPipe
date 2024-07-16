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
package_list = c('plyr','dplyr','ggplot2','ggridges', 'data.table', 'corrplot', 
                 'umap', 'magick', 'ggdendro', 'ecodist','ggbeeswarm', 
                 'ggrepel', 'ggthemes', 'foreach','reshape2','org.Hs.eg.db',
                 'clusterProfiler','pheatmap','limma','UpSetR','maSigPro','KSEAapp')
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
tryTo(paste0('INFO: Reading input file ', opt$pepfile),{
  dat <- fread(opt$pepfile)
  }, paste0('ERROR: problem trying to load ', opt$pepfile, ', does it exist?'))
  
tryTo(paste0('INFO: Massaging data from ', opt$pepfile, ' into a common style format for processing'), {
  dat <- phospho_standardize_format(dat)
  }, 'ERROR: failed! Check for missing/corrupt headers?')

#Identify unique PTM with the least missing values and highest median intensity
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
  
tryTo(paste0('INFO: Applying Filter Intensity > ',opt$minintensity),{
  dat_long <- dat_long[Intensity > opt$minintensity]
  phospho_dat_long <- phospho_dat_long[Intensity > opt$minintensity]
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


####DE##########################
if (!is.null(opt$design)) {
  DE_dir <- paste0(opt$outdir, '/Differential_Intensity/')
  if(! dir.exists(DE_dir)){
    dir.create(DI_dir, recursive = T)
  }
  EA_dir <- paste0(opt$outdir, '/Enrichiment_Analysis/')
  if(! dir.exists(EA_dir)){
    dir.create(EA_dir, recursive = T)
  }
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
    tryTo('INFO: Running differential intensity t-tests and pathway analysis',{
    }, 'ERROR: DE ttest failed!')
  }
  ##DE limma
  else if (opt$DE_method == 'limma') {
    Log2_phospho_dat=as.data.frame(phospho_dat) %>%
      mutate_if(is.numeric, ~ log2(. ))
    do_limma(Log2_DT = Log2_phospho_dat,
             design_matrix = design,
             DE_dir = DI_dir,
             EA_dir=EA_dir,
             lfc_threshold =opt$lfc_threshold,
             fdr_threshold =opt$fdr_threshold,
             enrich_pvalue = opt$enrich_pvalue )
  }
}
  



quit()

