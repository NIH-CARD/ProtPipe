#!/usr/bin/env Rscript
# R/4
#Peptidome analysis for FragPipe quantity estimates

#### ARG PARSING ###################################################################################
library(optparse)
pwd = getwd()
optparse_indent = '\n                '
option_list = list( 
    make_option(
        "--pepfile",
        default=NULL,
        help=paste(
            'Input file of Peptide Intensity (from FragPipe)',
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
            'Filter out samples with peptide counts > N standard deviations from the mean.',
            'Increase to higher values for greater tolerance of variance in peptide group counts.',
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
            'peptide groups in differential abundance analysis. Default: 2 (equivalent to',
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
    ),
    
    make_option(
      "--hla",
      default=NULL,
      help=paste(
        'The HLA typing information.',
        sep=optparse_indent
      )
    )
)



usage_string <- "Rscript %prog --pepfile [filename] --design [filename] [other options] "
opt <- parse_args(OptionParser(usage = usage_string, option_list))


source('src/functions.R')



# Set to TRUE when running interactively for debugging, to set test opts

badargs <- FALSE

if(opt$dry) {
    cat("INFO: Quitting due to --dry run\n")
    quit(status=0)
}

if(! opt$normalize %in% c('shift','scale','none')) {
    cat("ERROR: --normalize must be 'shift' 'scale' or 'none'\n")
    badargs <- TRUE
}

if (is.null(opt$pepfile)) {
    cat("ERROR: --pepfile <file> must be provided\n")
    badargs <- TRUE
}



if (badargs == TRUE) {
    quit(status=1)
}


#### PACKAGES ######################################################################################
package_list = c('ggplot2', 'data.table', 'corrplot', 'umap', 'magick', 'ggdendro', 'ecodist',
                 'ggbeeswarm', 'ggrepel', 'ggthemes', 'foreach','reshape2','org.Hs.eg.db',
                 'clusterProfiler','STRINGdb','eulerr')
cat("INFO: Loading required packages\n      ")
cat(paste(package_list, collapse='\n      ')); cat('\n')

defaultW <- getOption("warn"); options(warn = -1)   # Temporarily disable warnings for quiet loading
if(all(suppressWarnings(suppressMessages((lapply(package_list, require, character.only=TRUE)))))) {
    cat("INFO: All packages successfully loaded\n")
} else {
    cat("ERROR: One or more packages not available. Are you running this within the container?\n")
}
options(warn = defaultW)    # Turn warnings back on



opt$lfc_threshold <- log(opt$foldchange, base=2)
cat(paste0('INFO: LFC threshold of log2(Intensity) > ', opt$lfc_threshold, '\n'))
cat(paste0('INFO: FDR threshold of ', opt$fdr_threshold, '\n'))


#### MAKE DIRS #####################################################################################

QC_dir <- paste0(opt$outdir, '/QC/')
if(! dir.exists(QC_dir)){
    dir.create(QC_dir, recursive = T)
}


cluster_dir <- paste0(opt$outdir, '/Clustering/')
if(! dir.exists(cluster_dir)){
    dir.create(cluster_dir, recursive = T)
}

if (!is.null(opt$design)) {
  DI_dir <- paste0(opt$outdir, '/Differential_Intensity/')
  if(! dir.exists(DI_dir)){
    dir.create(DI_dir, recursive = T)
    }
  EA_dir <- paste0(opt$outdir, '/Enrichiment_Analysis/')
  if(! dir.exists(EA_dir)){
    dir.create(EA_dir, recursive = T)
  }
}

if (!is.null(opt$hla)) {
  MHC_dir <- paste0(opt$outdir, '/MHC_Prediction/')
  if(! dir.exists(MHC_dir)){
    dir.create(MHC_dir, recursive = T)
  }
}


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

## Plotting intensity distribution

tryTo('INFO: Plotting intensity distribution',{
  plot_pep_intensities(dat.long, QC_dir, 'intensities.pdf', plot_title='Normalized intensities')
      # plot_density(dat.long, QC_dir, 'intensity_density.pdf')
    # plot_density(dat.long.normalized, QC_dir, 'intensity_density_normalized.pdf')
}, 'ERROR: failed!')


## Normalization does not take place by default, and can be modified with the --normalize flag. See opts.
if (opt$normalize == 'none') {
    cat('Skipping median-normalization due to --normalize none\n')
} else {
    if (opt$normalize == 'shift') {

        tryTo('INFO: Calculating median-normalized intensities by shifting sample intensities',{
            dat.long <- shift_normalize_intensity(dat.long)
        }, 'ERROR: failed!')

        tryTo('INFO: Plotting shift-normalized intensity distributions',{
          plot_pg_intensities(dat.long, QC_dir, 'intensities_shift_normalized.pdf', plot_title='shift-normalized intensities')
        }, 'ERROR: failed!')

    } else if (opt$normalize == 'scale') {

        tryTo('INFO: Calculating median-normalized intensities by scaling sample intensities',{
            dat.long <- scale_normalize_intensity(dat.long)
        }, 'ERROR: failed!')

        tryTo('INFO: Plotting scale-normalized intensity distributions',{
            plot_pg_intensities(dat.long, QC_dir, 'intensities_scale_normalized.pdf', plot_title='scale-normalized intensities')
        }, 'ERROR: failed!')
    }

    tryTo('INFO: Re-generating wide table with normalized intensities',{
        original_colorder <- colnames(dat)
        dat <- dcast(dat.long, Protein_Group+Genes+First_Protein_Description~Sample, value.var='Intensity')
        setcolorder(dat, original_colorder)
    }, 'ERROR: failed!')
}


# Peptide_counts represents the distribution of peptide with Intensity > 0
# Visually, it is represented as a bar plot with x=sample, y=N, ordered by descending N
# Get counts of [N=unique gene groups with `Intensity` > 0]
tryTo('INFO: Cabulating peptide counts',{
  pep_counts <- dat.long[, .N, by=Sample]
  # Order samples by ascending counts
  ezwrite(pep_counts, QC_dir, 'peptide_group_nonzero_counts.tsv')
  plot_pep_counts(pep_counts, QC_dir, 'peptide_nonzero_counts.pdf')
}, 'ERROR: failed!')


tryTo('INFO: Plotting sample intensity correlations',{
    dat.correlations <- get_spearman(dat)
    ezwrite(dat.correlations, QC_dir, 'sample_correlation.tsv')
    plot_correlation_heatmap(dat.correlations, QC_dir, 'sample_correlation.pdf')
}, 'ERROR: failed!')

if (!is.null(opt$design)) {
  tryTo('INFO: Importing experimental design',{
    design <- fread(opt$design, header=TRUE)
    setnames(design, c('sample_name', 'condition', 'control'))
    }, 'ERROR: failed!')
}

if (!is.null(opt$exclude)) {
    tryTo(paste('INFO: excluding samples', opt$exclude),{
        design <- design[! sample_name %in% opt$exclude]
    }, 'ERROR: failed!')
}


#### CLUSTERING ####################################################################################


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




#### DIFFERENTIAL INTENSITY ########################################################################
if (!is.null(opt$design)) {
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


  tryTo('INFO: Running differential intensity t-tests and pathway analysis',{
    for (treatment in conditions) {
      treatment_sample_names <- intersect(colnames(dat), design[condition == treatment, sample_name])
        if (length(treatment_sample_names)==0) {next}
        control <- unique(design[condition == treatment, control])
        control_sample_names <- colnames(dat)[colnames(dat) %like% control]
        if (length(control_sample_names)==0) {next}
        if(treatment != control) {
            t_test <- do_t_test_APMS(dat, treatment_sample_names, control_sample_names)
            ezwrite(t_test[order(p.adj)], DI_dir, paste0(treatment, '_vs_', control, '.tsv'))
            plot_volcano(t_test, treatment, control, opt$log_base, opt$lfc_threshold, opt$fdr_threshold, DI_dir, opt$labelgene)
            enrich_pathway(t_test, treatment, control, EA_dir, opt$lfc_threshold, opt$fdr_threshold,opt$enrich_pvalue)
        }
    }
}, 'ERROR: failed!')
  
  

}

#### MHC BINDING PREDICTION########################################################################
if (!is.null(opt$hla)) {
  tryTo('INFO:Peptide/MHC I binding affinity prediction',{
    tryTo('INFO: Importing HLA_typing file',{
      hla_typing=fread(opt$hla, header=TRUE)
    }, 'ERROR: failed!')
    print(hla_typing[])
    
    for (sample in unique(hla_typing$sample_name)) {
      #Generate the MHCflurry input file
      get_mhcflurry_input(hla_typing,dat.long,sample,MHC_dir)
      # Using the system() function to run the mhcflurry
      runing_mhcflurry(sample)
      #Plot for the Peptide/MHC I binding affinity prediction
      plot_mhc_affinity(sample)
    }
  }, 'ERROR: failed!')
}




quit()

## TODO:
# Evaluate normalization with samples from very different batches?
# Calculate differential intensities using DESeq2?
# Function to compare q-values from spectronaut VS DIA-NN?
