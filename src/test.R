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


source('src/functions.R')
opt <- parse_args(OptionParser(option_list=option_list))
print(opt)