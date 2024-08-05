#!/usr/bin/env Rscript
# R/4
#proteomics analysis for DIA-NN and Spectronaut quantity estimates

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


if(! opt$normalize %in% c('none','mean','median')) {
    cat("ERROR: --normalize must be 'mean','median', or 'none'\n")
    badargs <- TRUE
}

if (is.null(opt$pgfile) && is.null(opt$pepfile)) {
    cat("ERROR: --pgfile <file> or --pepfile  <file> must be provided\n")
    badargs <- TRUE
}

#### Source the function ###
source('src/functions.R')
#### PACKAGES ######################################################################################
package_list = c('ggplot2', 'data.table','dplyr' ,'corrplot', 'umap', 
                 'magick', 'ggdendro', 'ecodist','ggbeeswarm',
                 'ggrepel', 'ggthemes', 'foreach','reshape2',
                 'org.Hs.eg.db','clusterProfiler','pheatmap','limma','DOSE')
cat("INFO: Loading required packages\n      ")
cat(paste(package_list, collapse='\n      ')); cat('\n')

defaultW <- getOption("warn"); options(warn = -1)   # Temporarily disable warnings for quiet loading
if(all((lapply(package_list, require, character.only=TRUE)))) {
    cat("INFO: All packages successfully loaded\n")
} else {
    cat("ERROR: One or more packages not available. Are you running this within the container?\n")
}
options(warn = defaultW)    # Turn warnings back on


###pgfile ############# 
if (!is.null(opt$pgfile)) {
  #### IMPORT AND FORMAT DATA#########################################################################
  tryTo(paste0('INFO: Reading input file ', opt$pgfile),{
    dat <- fread(opt$pgfile)
  }, paste0('ERROR: problem trying to load ', opt$pgfile, ', does it exist?'))
  
  tryTo(paste0('INFO: Massaging data from ', opt$pgfile, ' into a common style format for processing'), {
    dat <- standardize_format(dat)
  }, 'ERROR: failed! Check for missing/corrupt headers?')
  
  tryTo(paste0('INFO: Trimming extraneous column name info'), {
    setnames(dat, trim_colnames(dat))
    #set column order
    col_order=c(colnames(dat)[1:2],sort(colnames(dat)[3:ncol(dat)]))
    setcolorder(dat,col_order)
  }, 'ERROR: failed! Check for missing/corrupt headers?')
  
  # exclude samples in opt$exclude
  if (! is.null(opt$exclude)) {
    opt$exclude <- strsplit(opt$exclude, split=';')[[1]]
    tryTo(paste('INFO: excluding samples', opt$exclude),{
      for(i in opt$exclude) {
        dat[, (i) := NULL]
      }
    }, paste0('ERROR: problem trying to load ', opt$pgfile, ', does it exist?'))
  }
  
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
  QC_dir <- paste0(opt$outdir, '/QC/')
  if(! dir.exists(QC_dir)){
    dir.create(QC_dir, recursive = T)
  }
  ## Plotting intensity distribution
  
  tryTo('INFO: Plotting intensity distribution',{
    plot_pg_intensities(dat.long, QC_dir, 'intensities.pdf')
    # plot_density(dat.long, QC_dir, 'intensity_density.pdf')
    # plot_density(dat.long.normalized, QC_dir, 'intensity_density_normalized.pdf')
  }, 'ERROR: failed!')
  
  ## Normalization takes place by default, and can be modified with the --normalize flag. See opts.
  if (opt$normalize == 'none') {
    cat('Skipping median-normalization due to --normalize none\n')
  } else {
    if (opt$normalize == 'mean') {
      tryTo('INFO: Global Normalization by Mean',{
        dat.long <- mean_normalize_intensity(dat.long)
      }, 'ERROR: failed!')
      
      tryTo('INFO: Plotting mean-normalized intensity distributions',{
        plot_pg_intensities(dat.long, QC_dir, 'intensities_mean_normalized.pdf')
      }, 'ERROR: failed!')
    } else if (opt$normalize == 'median') {
      
      tryTo('INFO: Global Normalization by Median',{
        dat.long <- median_normalize_intensity(dat.long)
      }, 'ERROR: failed!')
      
      tryTo('INFO: Plotting median-normalized intensity distributions',{
        plot_pg_intensities(dat.long, QC_dir, 'intensities_median_normalized.pdf')
      }, 'ERROR: failed!')
    }
    
    tryTo('INFO: Re-generating wide table with normalized intensities',{
      original_colorder <- colnames(dat)
      dat <- dcast(dat.long, Protein_Group+Genes~Sample, value.var='Intensity')
      setcolorder(dat, original_colorder)
    }, 'ERROR: failed!')
  }
  
  # pgcounts represents the distribution of Protein Groups with Intensity > 0
  # Visually, it is represented as a bar plot with x=sample, y=N, ordered by descending N
  # Get counts of [N=unique gene groups with `Intensity` > 0]
  tryTo('INFO: Plotting protein group counts',{
    pgcounts=plot_pg_counts(DT.long = dat.long,
                   output_dir = QC_dir
                  )
  }, 'ERROR: failed!')
  
  tryTo('INFO: Plotting sample intensity correlations',{
    dat.correlations <- get_spearman(dat)
    ezwrite(dat.correlations, QC_dir, 'sample_correlation.tsv')
    plot_correlation_heatmap(dat.correlations, QC_dir, 'sample_correlation.pdf')
  }, 'ERROR: failed!')
  
  if (!is.null(opt$exclude)) {
    tryTo(paste('INFO: excluding samples', opt$exclude),{
      design <- design[! sample_name %in% opt$exclude]
    }, 'ERROR: failed!')
  }
  
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
  }
  
  
  ## Exclude samples with N protein groups < opt$sds away from mean
  ## Default value: 3 standard deviations, modifiable with --sds [N]
  tryTo('INFO: Identifying samples with protein group count outliers',{
    cat(paste0('INFO: defining outliers as samples with [N protein groups] > ', opt$sds, ' standard deviations from the mean\n'))
    stdev <- sd(pgcounts[,N])
    mean_count <- mean(pgcounts[,N])
    min_protein_groups <- floor(mean_count - (opt$sds * stdev))
    max_protein_groups <- ceiling(mean_count + (opt$sds * stdev))
    cat(paste0('INFO: Tolerating protein group counts in the range [', min_protein_groups,',',max_protein_groups,']'))
    low_count_samples <- as.character(pgcounts[N < min_protein_groups, Sample])
    high_count_samples <- as.character(pgcounts[N > max_protein_groups, Sample])
    if(length(low_count_samples)==0) {
      cat('\nINFO: No low group count samples to remove\n')
    } else {
      cat(paste0('\nINFO: Pruning low-count outlier ', low_count_samples))
      cat('\n\n')
      print(pgcounts[Sample %in% low_count_samples])
      cat('\n')
      dat[, c(low_count_samples) := NULL]    # remove sample columns from wide table
      dat.long <- dat.long[! (Sample %in% low_count_samples)] # remove rows from long table
    }
    if(length(high_count_samples)==0) {
      cat('INFO: No high group count samples to remove\n')
    } else {
      cat(paste0('\nINFO: Pruning high-count outlier ', high_count_samples))
      cat('\n')
      print(pgcounts[Sample %in% high_count_samples])
      dat[, c(high_count_samples) := NULL]    # remove sample columns from wide table
      dat.long <- dat.long[! (Sample %in% high_count_samples)] # remove rows from long table
    }
  }, 'ERROR: failed!')
  
  #### CLUSTERING ####################################################################################
  cluster_dir <- paste0(opt$outdir, '/Clustering/')
  if(! dir.exists(cluster_dir)){
    dir.create(cluster_dir, recursive = T)
  }
  # PCA
  tryTo('INFO: running PCA and plotting first two components',{
    pca <- get_PCs(dat)
    ezwrite(pca$components, cluster_dir, 'PCA.tsv')
    ezwrite(pca$summary, cluster_dir, 'PCA_summary.tsv')
    plot_PCs(pca, cluster_dir, 'PCA.pdf')
  }, 'ERROR: failed!')
  
  # Hierarchical Clustering
  tryTo('INFO: running Hierarchical Clustering',{
    plot_hierarchical_cluster(dat, cluster_dir)
  }, 'ERROR: failed!')
  
  # UMAP
  if ((ncol(dat)-3)>opt$neighbors) {
    tryTo('INFO: running UMAP',{
      umap <- get_umap(dat, opt$neighbors)
      ezwrite(umap, cluster_dir, 'UMAP.tsv')
      plot_umap(umap, cluster_dir, 'UMAP.pdf')
    }, 'ERROR: failed!')
  }
  
  #### DIFFERENTIAL INTENSITY########################################################################
  if (!is.null(opt$design)) {
    if (opt$DE_method == 'ttest') {
      DI_dir <- paste0(opt$outdir, '/ttest/Differential_Intensity/')
      if(! dir.exists(DI_dir)){
        dir.create(DI_dir, recursive = T)
      }
      
      EA_dir <- paste0(opt$outdir, '/ttest/Enrichiment_Analysis/')
      if(! dir.exists(EA_dir)){
        dir.create(EA_dir, recursive = T)
      }
      tryTo('INFO: Running differential intensity t-tests and pathway analysis',{
        for (treatment in conditions) {
          print(paste0(treatment, ' vs ', control, ' DE analysis'))
          treatment_sample_names <- intersect(colnames(dat), design[condition == treatment, sample_name])
          if (length(treatment_sample_names)==0) {next}
          control <- unique(design[condition == treatment, control])
          control_sample_names <- colnames(dat)[colnames(dat) %like% control]
          if (length(control_sample_names)==0) {next}
          if(treatment != control) {
            t_test <- do_t_test(DT = dat, treatment_samples = treatment_sample_names,control_samples =  control_sample_names)
            ezwrite(t_test[order(adj.P.Val)], DI_dir, paste0(treatment, '_vs_', control, '.tsv'))
            plot_volcano(DT.original = t_test, 
                         out_dir = DI_dir,
                         output_filename = paste0(treatment, '_vs_', control, '_ttest'),
                         label_col = 'Genes',
                         lfc_threshold = opt$lfc_threshold,
                         fdr_threshold = opt$fdr_threshold,
                         labelgene =  opt$labelgene)
            enrich_pathway(t_test, treatment, control, EA_dir, opt$lfc_threshold, opt$fdr_threshold,opt$enrich_pvalue)
          }
        }
      }, 'ERROR: failed!')
      
    }
    else if (opt$DE_method == 'limma') {
      DI_dir <- paste0(opt$outdir, '/limma/Differential_Intensity/')
      if(! dir.exists(DI_dir)){
        dir.create(DI_dir, recursive = T)
      }
      
      EA_dir <- paste0(opt$outdir, '/limma/Enrichiment_Analysis/')
      if(! dir.exists(EA_dir)){
        dir.create(EA_dir, recursive = T)
      }
      Log2_dat=dat
      Log2_dat=as.data.frame(Log2_dat)
      Log2_dat[, 3:ncol(Log2_dat)]=log2(Log2_dat[, 3:ncol(Log2_dat)]+1)
      
      do_limma(Log2_DT = Log2_dat,
               design_matrix = design,
               DE_dir = DI_dir,
               EA_dir=EA_dir,
               lfc_threshold =opt$lfc_threshold,
               fdr_threshold =opt$fdr_threshold,
               enrich_pvalue = opt$enrich_pvalue )
    }
  }
  
  #### HEATMAP ########################################################################
  plot_heatmap(dat)
  if (!is.null(opt$heatmap)) {
    plot_heatmap_subset(dat,opt$heatmap)
  }
}

##pepfile########  
if (!is.null(opt$pepfile)) {
  #### IMPORT AND FORMAT DATA#########################################################################
  tryTo(paste0('INFO: Reading input file ', opt$pepfile),{
    dat <- fread(opt$pepfile)
  }, paste0('ERROR: problem trying to load ', opt$pepfile, ', does it exist?'))
  
  tryTo(paste0('INFO: Massaging data from ', opt$pepfile, ' into a common style format for processing'), {
    dat <- standardize_format(dat)
  }, 'ERROR: failed! Check for missing/corrupt headers?')
  
  tryTo(paste0('INFO: Trimming extraneous column name info'), {
    setnames(dat, trim_colnames(dat))
    #set column order
    col_order=c(colnames(dat)[1:2],sort(colnames(dat)[3:ncol(dat)]))
    setcolorder(dat,col_order)
  }, 'ERROR: failed! Check for missing/corrupt headers?')
  
  # exclude samples in opt$exclude
  if (! is.null(opt$exclude)) {
    opt$exclude <- strsplit(opt$exclude, split=';')[[1]]
    tryTo(paste('INFO: excluding samples', opt$exclude),{
      for(i in opt$exclude) {
        dat[, (i) := NULL]
      }
    }, paste0('ERROR: problem trying to load ', opt$pepfile, ', does it exist?'))
  }
  
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
  QC_dir <- paste0(opt$outdir, '/QC/')
  if(! dir.exists(QC_dir)){
    dir.create(QC_dir, recursive = T)
  }
  ## Plotting intensity distribution
  
  tryTo('INFO: Plotting intensity distribution',{
    plot_pep_intensities(dat.long, QC_dir, 'intensities.pdf')
    # plot_density(dat.long, QC_dir, 'intensity_density.pdf')
    # plot_density(dat.long.normalized, QC_dir, 'intensity_density_normalized.pdf')
  }, 'ERROR: failed!')
  
  ## Normalization takes place by default, and can be modified with the --normalize flag. See opts.
  if (opt$normalize == 'none') {
    cat('Skipping median-normalization due to --normalize none\n')
  } else {
    if (opt$normalize == 'shift') {
      tryTo('INFO: Calculating median-normalized intensities by shifting sample intensities',{
        dat.long <- shift_normalize_intensity(dat.long)
      }, 'ERROR: failed!')
      
      tryTo('INFO: Plotting shift-normalized intensity distributions',{
        plot_pg_intensities(dat.long, QC_dir, 'intensities_shift_normalized.pdf')
      }, 'ERROR: failed!')
    } else if (opt$normalize == 'scale') {
      
      tryTo('INFO: Calculating median-normalized intensities by scaling sample intensities',{
        dat.long <- scale_normalize_intensity(dat.long)
      }, 'ERROR: failed!')
      
      tryTo('INFO: Plotting scale-normalized intensity distributions',{
        plot_pg_intensities(dat.long, QC_dir, 'intensities_scale_normalized.pdf')
      }, 'ERROR: failed!')
    }
    
    tryTo('INFO: Re-generating wide table with normalized intensities',{
      original_colorder <- colnames(dat)
      dat <- dcast(dat.long, Protein_Group+Genes+First_Protein_Description~Sample, value.var='Intensity')
      setcolorder(dat, original_colorder)
    }, 'ERROR: failed!')
  }
  
  # pgcounts represents the distribution of Protein Groups with Intensity > 0
  # Visually, it is represented as a bar plot with x=sample, y=N, ordered by descending N
  # Get counts of [N=unique gene groups with `Intensity` > 0]
  tryTo('INFO: Tabulating protein group counts',{
    pgcounts <- dat.long[, .N, by=Sample]
    # Order samples by ascending counts
    ezwrite(pgcounts, QC_dir, 'protein_group_nonzero_counts.tsv')
    plot_pep_counts(pgcounts, QC_dir, 'protein_group_nonzero_counts.pdf')
  }, 'ERROR: failed!')
  
  tryTo('INFO: Plotting sample intensity correlations',{
    dat.correlations <- get_spearman(dat)
    ezwrite(dat.correlations, QC_dir, 'sample_correlation.tsv')
    plot_correlation_heatmap(dat.correlations, QC_dir, 'sample_correlation.pdf')
  }, 'ERROR: failed!')
  
  if (!is.null(opt$exclude)) {
    tryTo(paste('INFO: excluding samples', opt$exclude),{
      design <- design[! sample_name %in% opt$exclude]
    }, 'ERROR: failed!')
  }
  
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
  }
  
  
  ## Exclude samples with N protein groups < opt$sds away from mean
  ## Default value: 3 standard deviations, modifiable with --sds [N]
  tryTo('INFO: Identifying samples with protein group count outliers',{
    cat(paste0('INFO: defining outliers as samples with [N protein groups] > ', opt$sds, ' standard deviations from the mean\n'))
    stdev <- sd(pgcounts[,N])
    mean_count <- mean(pgcounts[,N])
    min_protein_groups <- floor(mean_count - (opt$sds * stdev))
    max_protein_groups <- ceiling(mean_count + (opt$sds * stdev))
    cat(paste0('INFO: Tolerating protein group counts in the range [', min_protein_groups,',',max_protein_groups,']'))
    low_count_samples <- as.character(pgcounts[N < min_protein_groups, Sample])
    high_count_samples <- as.character(pgcounts[N > max_protein_groups, Sample])
    if(length(low_count_samples)==0) {
      cat('\nINFO: No low group count samples to remove\n')
    } else {
      cat(paste0('\nINFO: Pruning low-count outlier ', low_count_samples))
      cat('\n\n')
      print(pgcounts[Sample %in% low_count_samples])
      cat('\n')
      dat[, c(low_count_samples) := NULL]    # remove sample columns from wide table
      dat.long <- dat.long[! (Sample %in% low_count_samples)] # remove rows from long table
    }
    if(length(high_count_samples)==0) {
      cat('INFO: No high group count samples to remove\n')
    } else {
      cat(paste0('\nINFO: Pruning high-count outlier ', high_count_samples))
      cat('\n')
      print(pgcounts[Sample %in% high_count_samples])
      dat[, c(high_count_samples) := NULL]    # remove sample columns from wide table
      dat.long <- dat.long[! (Sample %in% high_count_samples)] # remove rows from long table
    }
  }, 'ERROR: failed!')
  
  #### CLUSTERING ####################################################################################
  cluster_dir <- paste0(opt$outdir, '/Clustering/')
  if(! dir.exists(cluster_dir)){
    dir.create(cluster_dir, recursive = T)
  }
  
  # PCA
  tryTo('INFO: running PCA and plotting first two components',{
    pca <- get_PCs(dat)
    ezwrite(pca$components, cluster_dir, 'PCA.tsv')
    ezwrite(pca$summary, cluster_dir, 'PCA_summary.tsv')
    plot_PCs(pca, cluster_dir, 'PCA.pdf')
  }, 'ERROR: failed!')
  
  # Hierarchical Clustering
  tryTo('INFO: running Hierarchical Clustering',{
    plot_hierarchical_cluster(dat, cluster_dir)
  }, 'ERROR: failed!')
  
  # UMAP
  if ((ncol(dat)-3)>opt$neighbors) {
    tryTo('INFO: running UMAP',{
      umap <- get_umap(dat, opt$neighbors)
      ezwrite(umap, cluster_dir, 'UMAP.tsv')
      plot_umap(umap, cluster_dir, 'UMAP.pdf')
    }, 'ERROR: failed!')
  }
  
  #### DIFFERENTIAL INTENSITY########################################################################
  if (!is.null(opt$design)) {
    if (opt$DE_method == 'ttest') {
      DI_dir <- paste0(opt$outdir, '/ttest/Differential_Intensity/')
      if(! dir.exists(DI_dir)){
        dir.create(DI_dir, recursive = T)
      }
      
      EA_dir <- paste0(opt$outdir, '/ttest/Enrichiment_Analysis/')
      if(! dir.exists(EA_dir)){
        dir.create(EA_dir, recursive = T)
      }
      tryTo('INFO: Running differential intensity t-tests and pathway analysis',{
        for (treatment in conditions) {
          print(paste0(treatment, ' vs ', control, ' DE analysis'))
          treatment_sample_names <- intersect(colnames(dat), design[condition == treatment, sample_name])
          if (length(treatment_sample_names)==0) {next}
          control <- unique(design[condition == treatment, control])
          control_sample_names <- colnames(dat)[colnames(dat) %like% control]
          if (length(control_sample_names)==0) {next}
          if(treatment != control) {
            t_test <- do_t_test(dat, treatment_sample_names, control_sample_names)
            ezwrite(t_test[order(p.adj)], DI_dir, paste0(treatment, '_vs_', control, '.tsv'))
            plot_volcano(DT.original = t_test,
                         out_dir = DI_dir,
                         output_filename = paste0(treatment, '_vs_', control, '_ttest'),
                         label_col = 'Peptide_Sequence',
                         lfc_threshold =opt$lfc_threshold, 
                         fdr_threshold = opt$fdr_threshold, 
                         labelgene =  opt$labelgene)
            enrich_pathway(t_test, treatment, control, EA_dir, opt$lfc_threshold, opt$fdr_threshold,opt$enrich_pvalue)
          }
        }
      }, 'ERROR: failed!')
    }
    else if (opt$DE_method == 'limma') {
      DI_dir <- paste0(opt$outdir, '/limma/Differential_Intensity/')
      if(! dir.exists(DI_dir)){
        dir.create(DI_dir, recursive = T)
      }
      
      EA_dir <- paste0(opt$outdir, '/limma/Enrichiment_Analysis/')
      if(! dir.exists(EA_dir)){
        dir.create(EA_dir, recursive = T)
      }
      Log2_dat=dat
      Log2_dat=as.data.frame(Log2_dat)
      Log2_dat[, 3:ncol(Log2_dat)]=log2(Log2_dat[, 3:ncol(Log2_dat)]+1)
      
      do_limma(Log2_DT = Log2_dat,
               design_matrix = opt$design,
               DE_dir = DI_dir,
               lfc_threshold =opt$lfc_threshold,
               fdr_threshold =opt$fdr_threshold,
               enrich_pvalue = opt$enrich_pvalue )
    }
    
    
  }
  
  #### HEATMAP ########################################################################
  plot_heatmap(dat)
  if (!is.null(opt$heatmap)) {
    plot_heatmap_subset(dat,opt$heatmap)
  }
  
}

##Done######## 
quit()
