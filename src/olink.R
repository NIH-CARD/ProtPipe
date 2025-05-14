##########################
#Olink analysis
#August-2024
#Ziyi Li
#CARD NIH
#########################

#### ARG PARSING ###################################################################################
library(optparse)
pwd = getwd()
optparse_indent = '\n                '
option_list = list( 
  make_option(
    "--input",
    default=NULL,
    help=paste(
      'Input file of NPX data (Olink)',
      'Required.',
      sep=optparse_indent
    )
  ),
  make_option(
    "--condition",
    dest='condition',
    default=NULL,
    help=paste(
      'Comma- or tab-delimited, three-column text file specifying the experimental design.',
      'File should contain headers. Header names do not matter; column order DOES matter.',
      'Columns order: <sample_name> <condition>',
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
      'Bridging normalization:  One of the dataframes is adjusted to another using overlapping samples (bridge samples)',
      'Subset normalization: A subset of samples is used to normalize two dataframes, one of which is used as a reference_project',
      'Reference median normalization: Works only on one dataframe.',
      'none: do not normalize',
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


#### PACKAGES PACKAGES AND SOURCE CODE ######################################################################################
package_list = c('OlinkAnalyze','ggplot2','ggdendro', 'reshape2','data.table',
                 'clusterProfiler','DOSE','enrichplot','org.Hs.eg.db',
                 'dplyr')
cat("INFO: Loading required packages\n      ")
cat(paste(package_list, collapse='\n      ')); cat('\n')

defaultW <- getOption("warn"); options(warn = -1)   # Temporarily disable warnings for quiet loading
if(all((lapply(package_list, require, character.only=TRUE)))) {
  cat("INFO: All packages successfully loaded\n")
} else {
  cat("ERROR: One or more packages not available. Are you running this within the container?\n")
}
options(warn = defaultW)    # Turn warnings back on

source('src/functions.R')

#### IMPORT AND FORMAT DATA#########################################################################

tryTo(paste0('INFO: Reading input file ', opt$input),{
  my_npx=read_NPX(opt$input)
}, paste0('ERROR: problem trying to load ', opt$input, ', does it exist?'))

if (!is.null(opt$condition)) {
  tryTo(paste0('INFO: Reading condition file ', opt$condition),{
    condition=read.csv(opt$condition)
  }, paste0('ERROR: problem trying to load ', opt$condition, ', does it exist?'))
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
  npx_wide <- my_npx |>
    dplyr::filter(SampleType == "SAMPLE") |>
    dplyr::filter(AssayType == "assay") |>
    dplyr::select(SampleID, UniProt,Assay, OlinkID, NPX) |>
    tidyr::pivot_wider(names_from = SampleID, values_from = NPX) |>
    dplyr::rename(Protein_Group = UniProt, Genes = Assay)
  npx_wide_all<- my_npx |>
    dplyr::filter(AssayType == "assay") |>
    dplyr::select(SampleID, UniProt,Assay, OlinkID, NPX) |>
    tidyr::pivot_wider(names_from = SampleID, values_from = NPX) |>
    dplyr::rename(Protein_Group = UniProt, Genes = Assay)
  ezwrite(x = npx_wide,
          output_dir = opt$outdir,
          output_filename = 'olink_data.tsv')
}, 'ERROR: failed! Check for missing/corrupt headers?')

#Converting to long format
tryTo(paste0('INFO: Converting to long format'), {
  dat_long=melt_intensity_table(npx_wide)
  dat_all_long=melt_intensity_table(npx_wide_all)
}, 'ERROR: failed! Check for missing/corrupt headers?')

#### QC ############################################################################################
## MAKE DIRS
QC_dir=paste0(opt$outdir, '/QC/')
if(! dir.exists(QC_dir)){
  dir.create(QC_dir, recursive = T)
}
## Plotting intensity distribution
tryTo('INFO: Plotting intensity distribution',{
  plot_olink_intensities(DT=my_npx,output_dir=QC_dir)
}, 'ERROR: failed!')


## Plotting the number of protein over than LOD
tryTo('INFO: Plotting the number of protein over than LOD',{
  olink_lod_fillter(DT=my_npx,output_dir=QC_dir)
}, 'ERROR: failed!')

## Plotting the number of protein over than NC
tryTo('INFO: Plotting the number of protein over than NC',{
  olink_data_fillter_nc(DT=my_npx,output_dir=QC_dir)
}, 'ERROR: failed!')


#### CLUSTERING ####################################################################################
cluster_dir=paste0(opt$outdir, '/Clustering/')
if(! dir.exists(cluster_dir)){
  dir.create(cluster_dir, recursive = T)
}

# Hierarchical Clustering
tryTo('INFO: running Hierarchical Clustering',{
  plot_hierarchical_cluster(DT = npx_wide, output_dir = cluster_dir)
}, 'ERROR: failed!')

# PCA
tryTo('INFO: running PCA and plotting first two components',{
  olink_plot_PCs(DT =npx_wide,DT_all=npx_wide_all,condition_file =condition ,output_dir = cluster_dir)
}, 'ERROR: failed!')


# UMAP
tryTo('INFO: running UMAP',{
  umap=soma_get_umap(dat,condition)
  ezwrite(umap, cluster_dir, 'UMAP.tsv')
  soma_plot_umap(umap, cluster_dir, 'UMAP.pdf')
}, 'ERROR: failed!')


#### DIFFERENTIAL INTENSITY########################################################################
if (all(is.na(condition$Sample_Group)))  {
  print('Condition file or adat SampleGroup info needed')
}else{
  DI_dir <- paste0(opt$outdir, '/Differential_Intensity/')
  if(! dir.exists(DI_dir)){
    dir.create(DI_dir, recursive = T)
  }
  
  EA_dir <- paste0(opt$outdir, '/Enrichiment_Analysis/')
  if(! dir.exists(EA_dir)){
    dir.create(EA_dir, recursive = T)
  }
  tryTo('INFO: Running differential intensity t-tests and pathway analysis',{
    conditions=na.omit(unique(condition$Sample_Group))
    for (i in 1:(length(conditions) - 1)) {
      for (j in (i + 1):length(conditions)) {
        treat <- conditions[i]
        control <- conditions[j]
        print(paste0(treat, ' vs ', control, ' DE analysis'))
        treatment_sample_names <- intersect(colnames(npx_wide), 
                                            condition$Sample_ID[condition$Sample_Group == treat])
        if (length(treatment_sample_names)==0) {next}
        control_sample_names <- intersect(colnames(npx_wide), 
                                          condition$Sample_ID[condition$Sample_Group == control])
        t_test <- olink_ttest(DT = npx_wide, 
                             treatment_samples = treatment_sample_names, 
                             control_samples = control_sample_names)
        ezwrite(t_test[order(t_test$adj.P.Val),], DI_dir, paste0(treat, '_vs_', control, '.tsv'))
        
        soma_plot_volcano(DT.original = t_test,
                          out_dir =DI_dir ,
                          output_filename = paste0(treat, '_vs_', control),
                          lfc_threshold = opt$lfc_threshold,
                          fdr_threshold = opt$fdr_threshold,
                          labelgene =opt$labelgene )
        soma_box_plot(soma_adat = my_adat,
                      treatment =treat,
                      control = control,
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
  
}