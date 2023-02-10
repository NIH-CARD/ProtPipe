#!/bin/bash
#bash for diann software and downstream analysis

#data(!!!this is the only thing you need to chage by youself!!!!)
raw_data_path='/data/liz36/data/diann_pro_pip/rawdata/'
fasta_file='/data/liz36/data/diann_pro_pip/rawdata/uniprot-proteome_Human_UP000005640_20191105.fasta'
outdir='results'
project_name='test'

#DIANN
sbatch --mem=200g --cpus-per-task=20 --time=2-12:00:00 jobscript 
