#!/bin/bash

#data(!!!this is the only thing you need to chage by youself!!!!)
raw_data_path='/data/liz36/data/diann_pro_pip/rawdata/'
fasta_file='/data/liz36/data/diann_pro_pip/rawdata/uniprot-proteome_Human_UP000005640_20191105.fasta'
outdir='./result'
project_name='test'
design_matrix='design_matrix.csv' #csv file of design matrix for comparison


#DIA-NN
if [[ ! -e $outdir ]]; then
    mkdir $outdir
module load diann
diann --dir $raw_data_path   --lib  --threads 24 --verbose 1 --out $outdir/$project_name.report.tsv --qvalue 0.01 --matrices --out-lib $outdir/$project_name.report-lib.tsv --gen-spec-lib --predictor --fasta $fasta_file --fasta-search --min-fr-mz 200 --max-fr-mz 2000 --met-excision --cut K*,R* --missed-cleavages 2 --min-pep-len 7 --max-pep-len 52 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4 --unimod4 --var-mods 5 --var-mod UniMod:35,15.994915,M --var-mod UniMod:1,42.010565,*n --monitor-mod UniMod:1 --reanalyse --relaxed-prot-inf --smart-profiling --peak-center --no-ifs-removal

#R script
module load R
pro_input=$outdir/$project_name.report.pg_matrix.tsv
pep_input=$outdir/$project_name.report.pr_matrix.tsv

Rscript diann_pro_pip.R --pro_input $pro_input --pep_input $pep_input -p $project_name -o $outdir --design_matrix $design_matrix -r $raw_data_path