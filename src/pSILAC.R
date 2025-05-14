#!/usr/bin/env Rscript
# R/4
# Ziyi Li @ NIH-CARD
# June-2024
# pSILAC analysis for Spectronaut quantity estimates

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
    "--normalize",
    default='none',
    type='character',
    help=paste('median:adjust sample intensities to match global median',
               'mean:adjust sample intensities to match global mean',
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

if (is.null(opt$pepfile)) {
  cat("ERROR: --pepfile  <file> must be provided\n")
  badargs <- TRUE
}

#### Source the function ###
source('src/functions.R')
#### PACKAGES ######################################################################################
package_list = c('ggplot2','ggridges', 'ggdendro','ggbeeswarm','ggrepel', 'ggthemes', 
                 'data.table','foreach','reshape2','corrplot', 'umap', 
                 'magick',  'ecodist',
                 'org.Hs.eg.db','clusterProfiler','DOSE','enrichplot',
                 'dplyr' )
cat("INFO: Loading required packages\n      ")
cat(paste(package_list, collapse='\n      ')); cat('\n')

defaultW <- getOption("warn"); options(warn = -1)   # Temporarily disable warnings for quiet loading
if(all(lapply(package_list, require, character.only=TRUE))) {
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
  
  tryTo(paste0('INFO: Standardize into a common style format for processing'), {
    psilac_dat <- psilac_standardize_format(dat)
  }, 'ERROR: failed! Check for missing/corrupt headers?')
  
  tryTo(paste0('INFO: Identify unique Peptide with the least missing values and highest median intensity'), {
    psilac_dat <- unique_data(DT = psilac_dat,
                              col = 'Peptide_Sequence',
                              output_dir = opt$outdir,
                              output_filename = 'psilac_data_output.tsv')
  }, 'ERROR: failed! Check for colnames?')
  #Converting to long format
  tryTo(paste0('INFO: Converting to long format'), {
    dat.long <- melt_intensity_table(psilac_dat)
  }, 'ERROR: failed! Check for missing/corrupt headers?')
  
  tryTo('INFO: Impute the NA as 0', {
    dat.long[is.na(dat.long)] = 0
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
  
  ## Plotting counts
  tryTo('INFO: Plotting counts',{
    pep_counts=plot_silac_counts(DT.long = dat.long,
                               output_dir = QC_dir,
                               intensity_cutoff = 0,
                               name = '',
                               height = 6, 
                               width = 12)
  }, 'ERROR: failed!')
  
  ## Exclude samples with N protein groups < opt$sds away from mean
  ## Default value: 3 standard deviations, modifiable with --sds [N]
  tryTo('INFO: Identifying samples with protein group count outliers',{
    cat(paste0('INFO: defining outliers as samples with [N protein groups] > ', opt$sds, ' standard deviations from the mean\n'))
    stdev <- sd(pgcounts$counts)
    mean_count <- mean(pgcounts$counts)
    min_protein_groups <- floor(mean_count - (opt$sds * stdev))
    max_protein_groups <- ceiling(mean_count + (opt$sds * stdev))
    cat(paste0('INFO: Tolerating protein group counts in the range [', min_protein_groups,',',max_protein_groups,']'))
    low_count_samples <- as.character(pgcounts$Sample[pgcounts$counts < min_protein_groups])
    high_count_samples <- as.character(pgcounts$Sample[pgcounts$counts > max_protein_groups])
    if(length(low_count_samples)==0) {
      cat('\nINFO: No low group count samples to remove\n')
    } else {
      cat(paste0('\nINFO: Runing low-count outlier ', low_count_samples))
      cat('\n\n')
      print(pgcounts[which(pgcounts$Sample %in% low_count_samples),])
      cat('\n')
      psilac_dat[, c(low_count_samples) ] = NULL   # remove sample columns from wide table
      dat.long <- dat.long[! (Sample %in% low_count_samples)] # remove rows from long table
    }
    if(length(high_count_samples)==0) {
      cat('INFO: No high group count samples to remove\n')
    } else {
      cat(paste0('\nINFO: Pruning high-count outlier ', high_count_samples))
      cat('\n')
      print(pgcounts[pgcounts$Sample %in% high_count_samples,])
      psilac_dat[, c(high_count_samples) ] = NULL   # remove sample columns from wide table
      dat.long <- dat.long[! (Sample %in% high_count_samples)] # remove rows from long table
    }
    if(length(high_count_samples) + length(low_count_samples) >0){
      ezwrite(x = psilac_dat,output_dir = opt$outdir ,output_filename = 'psilac_data_output.tsv')
      tryTo('INFO: Plotting counts of the filltered files',{
        plot_silac_counts(DT.long = dat.long,
                          output_dir = QC_dir,
                          intensity_cutoff = 0,
                          name = 'filtered_',
                          height = 6, 
                          width = 12)
      }, 'ERROR: failed!')
    }
  }, 'ERROR: failed!')

  ## Plotting intensity distribution
  tryTo('INFO: Plotting intensity distribution',{
    plot_silac_pep_intensity(DT.long = dat.long,
                             output_dir = QC_dir,
                             name = '',
                             height = 10, 
                             width = 12)
  }, 'ERROR: failed!')
  
  #### Cluster analysis ############################################
  cluster_dir <- paste0(opt$outdir, '/Clustering/')
  if(! dir.exists(cluster_dir)){
    dir.create(cluster_dir, recursive = T)
  }
  
  ## PCA
  tryTo('INFO: Plotting PCA',{
    psilac_pca(DT = psilac_dat,
               output_dir = cluster_dir,
               output_filename = 'PCA')
  }, 'ERROR: failed!')
  
  ## UMAP
  tryTo('INFO: Plotting UMAP',{
    psilac_umap(DT = psilac_dat,
               output_dir = cluster_dir,
               output_filename = 'UMAP')
  }, 'ERROR: failed!')
  
  
  #### Caculate half life #############################
  HL_dir <- paste0(opt$outdir, '/HALF_Life/')
  if(! dir.exists(HL_dir)){
    dir.create(HL_dir, recursive = T)
  }
  ## Caculate the remaining light peptide by L/(L+H)
  degradation_ratio=L_remain_channel(DT = psilac_dat,
                                     output_dir = HL_dir ,
                                     output_filename ='remaining_light_peptide_ratio' )
  ## Caculate the remaining heavy peptide by H/(L+H)
  synthesis_ratio=H_remain_channel(DT = psilac_dat,
                                   output_dir = HL_dir ,
                                   output_filename ='remaining_heavy_peptide_ratio' )
  
  ## Caculate the remaining heavy peptide and light peptide ration
  
  
}

if (!is.null(opt$pgfile)) {
  #### IMPORT AND FORMAT DATA#########################################################################
  tryTo(paste0('INFO: Reading input file ', opt$pepfile),{
    dat <- fread(opt$pgfile)
  }, paste0('ERROR: problem trying to load ', opt$pepfile, ', does it exist?'))
  
  tryTo(paste0('INFO: Standardize into a common style format for processing'), {
    psilac_dat <- psilac_standardize_format(dat)
  }, 'ERROR: failed! Check for missing/corrupt headers?')
  
  tryTo(paste0('INFO: Identify unique genes with the least missing values and highest median intensity'), {
    psilac_dat <- unique_data(DT = psilac_dat,
                              col = 'Genes',
                              output_dir = opt$outdir,
                              output_filename = 'psilac_pro_data_output.tsv')
  }, 'ERROR: failed! Check for colnames?')
  #Converting to long format
  tryTo(paste0('INFO: Converting to long format'), {
    dat.long <- melt_intensity_table(psilac_dat)
  }, 'ERROR: failed! Check for missing/corrupt headers?')
  
  tryTo('INFO: Impute the NA as 0', {
    dat.long[is.na(dat.long)] = 0
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
  
  ## Plotting counts
  tryTo('INFO: Plotting counts',{
    pgcounts=plot_silac_pg_counts(DT.long = dat.long,
                               output_dir = QC_dir,
                               intensity_cutoff = 0,
                               name = 'pro',
                               height = 6, 
                               width = 12)
  }, 'ERROR: failed!')
  
  ## Exclude samples with N protein groups < opt$sds away from mean
  ## Default value: 3 standard deviations, modifiable with --sds [N]
  tryTo('INFO: Identifying samples with protein group count outliers',{
    cat(paste0('INFO: defining outliers as samples with [N protein groups] > ', opt$sds, ' standard deviations from the mean\n'))
    stdev <- sd(pgcounts$counts)
    mean_count <- mean(pgcounts$counts)
    min_protein_groups <- floor(mean_count - (opt$sds * stdev))
    max_protein_groups <- ceiling(mean_count + (opt$sds * stdev))
    cat(paste0('INFO: Tolerating protein group counts in the range [', min_protein_groups,',',max_protein_groups,']'))
    low_count_samples <- as.character(pgcounts$Sample[pgcounts$counts < min_protein_groups])
    high_count_samples <- as.character(pgcounts$Sample[pgcounts$counts > max_protein_groups])
    if(length(low_count_samples)==0) {
      cat('\nINFO: No low group count samples to remove\n')
    } else {
      cat(paste0('\nINFO: Runing low-count outlier ', low_count_samples))
      cat('\n\n')
      print(pgcounts[which(pgcounts$Sample %in% low_count_samples),])
      cat('\n')
      psilac_dat[, c(low_count_samples) ] = NULL   # remove sample columns from wide table
      dat.long <- dat.long[! (Sample %in% low_count_samples)] # remove rows from long table
    }
    if(length(high_count_samples)==0) {
      cat('INFO: No high group count samples to remove\n')
    } else {
      cat(paste0('\nINFO: Pruning high-count outlier ', high_count_samples))
      cat('\n')
      print(pgcounts[pgcounts$Sample %in% high_count_samples,])
      psilac_dat[, c(high_count_samples) ] = NULL   # remove sample columns from wide table
      dat.long <- dat.long[! (Sample %in% high_count_samples)] # remove rows from long table
    }
    if(length(high_count_samples) + length(low_count_samples) >0){
      ezwrite(x = psilac_dat,output_dir = opt$outdir ,output_filename = 'psilac_data_output.tsv')
      tryTo('INFO: Plotting counts of the filltered files',{
        plot_silac_counts(DT.long = dat.long,
                          output_dir = QC_dir,
                          intensity_cutoff = 0,
                          name = 'filtered_',
                          height = 6, 
                          width = 12)
      }, 'ERROR: failed!')
    }
  }, 'ERROR: failed!')
  
  ## Plotting intensity distribution
  tryTo('INFO: Plotting intensity distribution',{
    plot_silac_pg_intensity(DT.long = dat.long,
                             output_dir = QC_dir,
                             height = 10, 
                             width = 12)
  }, 'ERROR: failed!')
  
  #### Cluster analysis ############################################
  cluster_dir <- paste0(opt$outdir, '/Clustering/')
  if(! dir.exists(cluster_dir)){
    dir.create(cluster_dir, recursive = T)
  }
  
  ## PCA
  tryTo('INFO: Plotting PCA',{
    psilac_pca(DT = psilac_dat,
               output_dir = cluster_dir,
               output_filename = 'PCA')
  }, 'ERROR: failed!')
  
  ## UMAP
  tryTo('INFO: Plotting UMAP',{
    psilac_umap(DT = psilac_dat,
                output_dir = cluster_dir,
                output_filename = 'UMAP')
  }, 'ERROR: failed!')
  
  
  #### Caculate half life #############################
  HL_dir <- paste0(opt$outdir, '/HALF_Life/')
  if(! dir.exists(HL_dir)){
    dir.create(HL_dir, recursive = T)
  }
  ## Caculate the remaining light peptide by L/(L+H)
  degradation_ratio=L_remain_channel(DT = psilac_dat,
                                     output_dir = HL_dir ,
                                     output_filename ='remaining_light_peptide_ratio' )
  ## Caculate the remaining heavy peptide by H/(L+H)
  synthesis_ratio=H_remain_channel(DT = psilac_dat,
                                   output_dir = HL_dir ,
                                   output_filename ='remaining_heavy_peptide_ratio' )
  
  ## Caculate the remaining heavy peptide and light peptide ration
  
  
}


quit()
