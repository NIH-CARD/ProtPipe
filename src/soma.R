#################
#somalogic data 
#Ziyi Li
#2024-05-10
########################
#### ARG PARSING ###################################################################################
library(optparse)
pwd = getwd()
optparse_indent = '\n                '
option_list = list( 
  make_option(
    "--input",
    default=NULL,
    help=paste(
      'Input file of protein group intensity (Somalogic)',
      'Required.',
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
    "--lfc_threshold",
    dest = 'lfc_threshold',
    default=0,
    type='numeric',
    help=paste(
      'The cutoff of log2 fold change for DE analysis.',
      'Default: 0',
      sep=optparse_indent
    )
  ),
  make_option(
    "--fdr_threshold",
    dest = 'fdr_threshold',
    default=0.01,
    type='numeric',
    help=paste(
      'The cutoff of p-value for DE analysis.',
      'Default: 0.01',
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

usage_string="Rscript %prog --input [filename]"
opt=parse_args(OptionParser(usage = usage_string, option_list))


#Package and Source Code#######################################################################
source('src/functions.R')
package_list = c('SomaDataIO' ,'data.table', 
                 'reshape2','ggplot2','umap','ggdendro','ggrepel',
                 'clusterProfiler','org.Hs.eg.db','DOSE','enrichplot',
                 'dplyr',
                 'pheatmap')
cat("INFO: Loading required packages\n      ")
cat(paste(package_list, collapse='\n      ')); cat('\n')
defaultW=getOption("warn"); options(warn = -1)   # Temporarily disable warnings for quiet loading
if(all((lapply(package_list, require, character.only=TRUE)))) {
  cat("INFO: All packages successfully loaded\n")
} else {
  cat("ERROR: One or more packages not available. Are you running this within the container?\n")
  quit()
}


#### IMPORT AND FORMAT DATA#########################################################################
tryTo(paste0('INFO: Reading input file ', opt$input),{
  my_adat=read_adat(opt$input)
}, paste0('ERROR: problem trying to load ', opt$input, ', does it exist?'))

if (!is.null(opt$condition)) {
  tryTo(paste0('INFO: Reading condition file ', opt$condition),{
    condition=fread(opt$condition)
  }, paste0('ERROR: problem trying to load ', opt$input, ', does it exist?'))
}else{
  condition <- my_adat %>%
    select_if(~ !is.numeric(.))
}

## MAKE DIRS
if(! dir.exists(opt$outdir)){
  dir.create(opt$outdir, recursive = T)
}

#Writing the standard TSV output 
tryTo(paste0('INFO: Writing the standard CSV output '), {
  dat=soma_sample_out(DT = my_adat)
  dat_all=soma_all_output(DT = my_adat,
                          output_dir = opt$outdir)
  dat_filter=Buffer_filter(DT = dat_all,output_dir = opt$outdir)
  if (!is.null(opt$condition)){
    tryTo(paste0('INFO: Reading condition file ',opt$condition),{
      condition_file=read.csv(opt$condition)
    }, paste0('ERROR: problem trying to load ', opt$condition, ', does it exist?'))
  }else{condition=condition=my_adat%>%
    select(-contains("seq."))}
  
}, 'ERROR: failed! Check for missing/corrupt headers?')

#Converting to long format
tryTo(paste0('INFO: Converting to long format'), {
  dat_long=melt_intensity_table(dat)
  dat_filter_long=melt_intensity_table(dat_filter)
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
  pgcounts=dat_filter_long[, .N, by=Sample]
  # Order samples by ascending counts
  ezwrite(x = pgcounts, 
          output_dir = QC_dir, 
          output_filename = 'protein_group_counts_over_buffer.tsv')
  soma_plot_counts(counts = pgcounts,
                   condition_file = condition,
                   output_dir = QC_dir)
 }, 'ERROR: failed!')

## Plotting intensity distribution
tryTo('INFO: Plotting intensity distribution',{
  plot_pg_intensities(dat_long, QC_dir, 'intensities.pdf')
  plot_pg_intensities(dat_all_long, QC_dir, 'all_intensities.pdf')
}, 'ERROR: failed!')

## Missing value check and plot
tryTo('INFO: Missing value check and plot',{
  miss_value_plot(DT = dat_filter,outdir = QC_dir)
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
  soma_plot_PCs(pca, cluster_dir, 'PCA.pdf')
}, 'ERROR: failed!')


# UMAP
tryTo('INFO: running UMAP',{
  umap=soma_get_umap(dat,condition)
  ezwrite(umap, cluster_dir, 'UMAP.tsv')
  soma_plot_umap(umap, cluster_dir, 'UMAP.pdf')
}, 'ERROR: failed!')

#### DIFFERENTIAL INTENSITY########################################################################
DI_dir <- paste0(opt$outdir, '/Differential_Intensity/')
if(! dir.exists(DI_dir)){
  dir.create(DI_dir, recursive = T)
}

EA_dir <- paste0(opt$outdir, '/Enrichiment_Analysis/')
if(! dir.exists(EA_dir)){
  dir.create(EA_dir, recursive = T)
}
tryTo('INFO: Running differential intensity t-tests and pathway analysis',{
  conditions=na.omit(unique(condition$SampleGroup))
  for (i in 1:(length(conditions) - 1)) {
    for (j in (i + 1):length(conditions)) {
      treat <- conditions[i]
      control <- conditions[j]
      print(paste0(treat, ' vs ', control, ' DE analysis'))
      treatment_sample_names <- intersect(colnames(dat), condition$SampleId[condition$SampleGroup == treat])
      if (length(treatment_sample_names)==0) {next}
      control_sample_names <- intersect(colnames(dat), condition$SampleId[condition$SampleGroup == control])
      t_test <- do_t_test(DT = dat, 
                          treatment_samples = treatment_sample_names, 
                          control_samples = control_sample_names)
      dat_Buffer_Calibrator <- dat_all %>% 
        select(AptName, Buffer, Calibrator)
      t_test <- merge(dat_Buffer_Calibrator, t_test, by = "AptName")%>%
        select(AptName, TargetFullName, Protein_Group, Genes, Buffer, Calibrator, 
               treatment_estimate, control_estimate, logFC, P.Value, adj.P.Val)
      # Rename columns
      colnames(t_test)[colnames(t_test) == "treatment_estimate"] <- paste0(treat,'_estimate')
      colnames(t_test)[colnames(t_test) == "control_estimate"] <- paste0(control,'_estimate')
      
      ezwrite(t_test[order(t_test$adj.P.Val),], DI_dir, paste0(treat, '_vs_', control, '_ttest.tsv'))
      soma_plot_volcano(DT.original = t_test,
                        out_dir =DI_dir ,
                        output_filename = paste0(treat, '_vs_', control),
                        lfc_threshold = opt$lfc_threshold,
                        fdr_threshold = opt$fdr_threshold,
                        labelgene =opt$labelgene )
      soma_box_plot(soma_adat = my_adat,
                    out_dir = DI_dir,
                    gene = opt$labelgene,
                    t_test = t_test,
                    levels = NULL,
                    color = NULL)
     
      enrich_pathway(DT.original = t_test, 
                     treatment = treat, 
                     control = control, 
                     outdir = EA_dir, 
                     lfc_threshold = opt$lfc_threshold,
                     fdr_threshold = opt$fdr_threshold,
                     enrich_pvalue = opt$enrich_pvalue)
      
    }
  }
}, 'ERROR: failed!')
#### HEATMAP ########################################################################
tryTo('INFO: Plot heatmap',{
  soma_plot_heatmap(DT_heatmap = dat,condition_file = condition)
  if (!is.null(opt$heatmap)) {
    soma_plot_heatmap_subset(DT_heatmap = dat,
                             condition_file = condition,
                             gene_subset = opt$heatmap)
  }
}, 'ERROR: failed!')


quit()


