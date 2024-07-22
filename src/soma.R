#################
#somalogic data 
#Ziyi Li
#2024-05-10
########################

#package and source code
package_list = c('SomaDataIO', 'dplyr', 'data.table', 
                 'reshape2','ggplot2','umap','ggdendro')
cat("INFO: Loading required packages\n      ")
cat(paste(package_list, collapse='\n      ')); cat('\n')
defaultW=getOption("warn"); options(warn = -1)   # Temporarily disable warnings for quiet loading
if(all((lapply(package_list, require, character.only=TRUE)))) {
  cat("INFO: All packages successfully loaded\n")
} else {
  cat("ERROR: One or more packages not available. Are you running this within the container?\n")
  quit()
}
source('src/soma_function.R')
#
#### ARG PARSING ###################################################################################
library(optparse)
pwd = getwd()
optparse_indent = '\n                '
option_list = list( 
  make_option(
    "--input",
    default=NULL,
    help=paste(
      'Input file of Protein Group Intensity (Somalogic)',
      'Required.',
      sep=optparse_indent
    )
  ),
  make_option(
    "--n",
    dest='out_name',
    default='out',
    help=paste(
      'name of outfile',
      'Defaults:out.tsv',
      sep=optparse_indent
    )
  ),
  make_option(
    "--condition",
    dest='condition',
    default=NULL,
    help=paste(
      'Sample condition',
      sep=optparse_indent
    )
  ),
  make_option(
    "--out",
    dest="outdir",
    default='output', 
    help=paste(
      'Directory to direct all output. Directory will be created if does not exist.',
      'Defaults to the current working directory:',
      pwd,
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

usage_string="Rscript %prog --input [filename] -n [filename] --design [filename] --condition [other options] "
opt=parse_args(OptionParser(usage = usage_string, option_list))


#### IMPORT AND FORMAT DATA#########################################################################
tryTo(paste0('INFO: Reading input file ', opt$input),{
  my_adat=read_adat(opt$input)
}, paste0('ERROR: problem trying to load ', opt$input, ', does it exist?'))

tryTo(paste0('INFO: Reading condition file ',opt$condition),{
  condition_file=read.csv(opt$condition)
}, paste0('ERROR: problem trying to load ', opt$condition, ', does it exist?'))

## MAKE DIRS
if(! dir.exists(opt$outdir)){
  dir.create(opt$outdir, recursive = T)
}

#Writing the standard CSV output 
tryTo(paste0('INFO: Writing the standard CSV output '), {
  dat_all=soma_all_output(DT = my_adat,
                          output_dir = opt$outdir)
  dat=soma_sample_out(DT = my_adat)
  dat_fillter=Buffer_filter(DT = dat_all,output_dir = opt$outdir)
  condition=condition=my_adat%>%
    select(-contains("seq."))
}, 'ERROR: failed! Check for missing/corrupt headers?')

#Converting to long format
tryTo(paste0('INFO: Converting to long format'), {
  dat_long=melt_intensity_table(dat)
  dat_fillter_long=melt_intensity_table(dat_fillter)
  dat_all_long=melt_intensity_table(dat_all)
}, 'ERROR: failed! Check for missing/corrupt headers?')


#### QC ############################################################################################
## MAKE DIRS
QC_dir=paste0(opt$outdir, '/QC/')
if(! dir.exists(QC_dir)){
  dir.create(QC_dir, recursive = T)
}
## Get counts of [N=unique gene groups with `Intensity` > 0]
tryTo('INFO: Tabulating protein group counts',{
  pgcounts=dat_fillter_long[, .N, by=Sample]
  # Order samples by ascending counts
  ezwrite(x = pgcounts, 
          output_dir = QC_dir, 
          output_filename = 'protein_group_counts_over_buffer.tsv')
  plot_pg_counts(DT = pgcounts, 
                 output_dir = QC_dir, 
                 output_filename = 'protein_group_counts_over_buffer.pdf')
  soma_plot_counts(counts = pgcounts,
                   condition_file = condition,
                   output_dir = QC_dir)
}, 'ERROR: failed!')

## Plotting intensity distribution
tryTo('INFO: Plotting intensity distribution',{
  plot_pg_intensities(dat_long, QC_dir, 'intensities.pdf')
  plot_pg_intensities(dat_all_long, QC_dir, 'all_intensities.pdf')
  
}, 'ERROR: failed!')

#### CLUSTERING ####################################################################################
cluster_dir=paste0(opt$outdir, '/Clustering/')
if(! dir.exists(cluster_dir)){
  dir.create(cluster_dir, recursive = T)
}

# Hierarchical Clustering
tryTo('INFO: running Hierarchical Clustering',{
  plot_hierarchical_cluster(dat, cluster_dir)
}, 'ERROR: failed!')

# PCA
tryTo('INFO: running PCA and plotting first two components',{
  pca=soma_get_PCs(dat,condition)
  pca$summary=merge(pca$summary,condition,by.x='Sample',by.y='SampleId')
  ezwrite(pca$components, cluster_dir, 'PCA.tsv')
  ezwrite(pca$summary, cluster_dir, 'PCA_summary.tsv')
  soma_plot_PCs(pca, cluster_dir, 'PCA.pdf')
}, 'ERROR: failed!')


# UMAP
if ((ncol(dat)-3)>opt$neighbors) {
  tryTo('INFO: running UMAP',{
    umap=get_umap(dat)
    ezwrite(umap, cluster_dir, 'UMAP.tsv')
    plot_umap(umap, cluster_dir, 'UMAP.pdf')
  }, 'ERROR: failed!')
}

#### DIFFERENTIAL INTENSITY########################################################################
#### DIFFERENTIAL INTENSITY########################################################################
if (!is.null(opt$design)) {
  if (opt$DE_method == 'ttest') {
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
          plot_volcano_pep(t_test, treatment, control, opt$log_base, opt$lfc_threshold, opt$fdr_threshold, DI_dir, opt$labelgene)
          enrich_pathway(t_test, treatment, control, EA_dir, opt$lfc_threshold, opt$fdr_threshold,opt$enrich_pvalue)
        }
      }
    }, 'ERROR: failed!')
  }else if (opt$DE_method == 'limma') {
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

quit()
