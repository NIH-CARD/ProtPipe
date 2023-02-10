# README

![workflow-image](src/workflow.png)

This repository facilitates quick, simple, and reproducible access to Data Independent Acquisition (DIA) proteomics workflows with minimal command line experience.

We include a convenient wrapper script for running DIA-NN inside a pre-built singularity image to first estimate protein abundance from raw mass spec output. Protein abundance estimates (accepting estimates from both `DIA-NN` and `Spectronaut`) can  be processed in a preconfigured R environment, generating QC reports, and various analyses and visualizations.


# Quick-start

1. Ensure `singularity` is installed and accessible on your system. Many HPCs (including NIH Biowulf) come with this pre-installed as a module. If your HPC has singularity installed, it will be automatically detected and loaded when necessary.
2. Clone this repository, i.e. `git clone https://github.com/cory-weller/ProtPipe.git`
3. If you are predicting protein abundances from raw mass spec output, look over and edit any custom `DIA-NN` parameters inside [`config.txt`](config.txt). You can either edit `config.txt` directly (and it will be used by default), or make a copy and save it to a different file name, then reference it with `--config newfilename.txt` when running the wrapper script.
4. Use the wrapper [`run-diann.sh`](src/run-diann.sh), either submitting to SLURM or running it directly.
5. R processing TBI

```bash
# Run directly
src/run-diann.sh \
    --fasta infile.fasta \
    --mzml infile.mzML \
    --out test_output &

# Submit to SLURM
sbatch src/run-diann.sh \
    --fasta example/uniprot-proteome_Human_UP000005640_20191105.fasta \
    --mzml example/raw_MS_mzML/HREC_ETIS_2.mzML \
    --out HREC_ETIS_2
```



# Installing Singularity

This workflow requires that [`Singularity`](https://sylabs.io/singularity) be available, which runs natively on a Linux system. `Singularity` is containerization software that allows an entire pre-configured computing environment to be accessed--reducing installation headaches and improving reproducibility. 

*We highly recommend making use a workstation or HPC with a native Linux installation.* Not only does this simplify the usage of `singularity`, it also would likely provide greater resources for DIA-NN's intensive computation.

To run on your personal/local non-Linux machine, Mac users need to first install a number of dependencies. Windows users would either need to use a virtual machine, or run things through the Windows Subsystem for Linux (WSL). Explaining the installation of `singularity` on these non-Linux systems is beyond the scope of this guide, so we defer to [the documentation here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

# Predicting Protein Abundances (running DIA-NN)

```bash
# Submit to SLURM
sbatch src/run-diann.sh \
    --config config.txt \
    --mzml ./example/raw_MS_mzML/HREC_ETIS_1.mzML \
    --fasta ./example/uniprot-proteome_Human_UP000005640_20191105.fasta \
    --out example/

# Run Locally
src/run-diann.sh \
    --config config.txt \
    --mzml ./example/raw_MS_mzML/HREC_ETIS_1.mzML \
    --fasta ./example/uniprot-proteome_Human_UP000005640_20191105.fasta \
    --out example/
```

<details><summary>Re-running on generated spectral library</summary>

For regenerating final outputs without the long computational steps. Requires the .speclib files.

```
singularity exec \
--cleanenv -H /home/wellerca/ProtPipe ./src/diann-1.8.1.sif diann \
--fasta example/uniprot-proteome_Human_UP000005640_20191105.fasta \
--reannotate \
--f example/raw_MS_mzML/HREC_ETIS_2.mzML \
--threads 4 \
--out-lib test \
--qvalue 0.01 \
--min-fr-mz 200 \
--max-fr-mz 2000  \
--cut K*,R* \
--missed-cleavages 2 \
--min-pep-len 7 \
--max-pep-len 52 \
--min-pr-mz 300 \
--max-pr-mz 1800 \
--min-pr-charge 1 \
--max-pr-charge 4 \
--var-mods 5 \
--monitor-mod UniMod:1 \
--var-mod UniMod:35,15.994915,M \
--var-mod UniMod:1,42.010565,*n \
--smart-profiling \
--peak-center \
--no-ifs-removal \
--met-excision  \
--matrices \
--lib test.speclib \
--out test
```

</details>

# Subsetting mzML file for testing purposes
from [here](https://rformassspectrometry.github.io/Spectra/articles/Spectra.html#exporting-spectra):
```R
# In R/4.2
BiocInstaller::install('mzR')
BiocInstaller::install('Spectra')

library(Spectra)

# load mzML using `Spectra`
sp <- Spectra('HREC_ETIS_1.mzML')

# write first 800 records using `export` and `MsBackendMzR()`
export(sp[1:8000], MsBackendMzR(), file='test.mzML')
```

# Processing DIA Estimates

# To Do
* update R in singularity image? (as `Spectra` package requires R/4.2)




<details><summary>extra to be changed</summary>
```bash
# Submit to SLURM
sbatch src/process-dia.sh \
    --config config.txt \
    --mzml ./example/raw_MS_mzML/HREC_ETIS_1.mzML \
    --fasta ./example/uniprot-proteome_Human_UP000005640_20191105.fasta \
    --out example/

# Run Locally
src/run-diann.sh \
    --config config.txt \
    --mzml ./example/raw_MS_mzML/HREC_ETIS_1.mzML \
    --fasta ./example/uniprot-proteome_Human_UP000005640_20191105.fasta \
    --out example/
```

Executing [`run.sh`](src/run.sh) runs the pipeline outlind below. Briefly, it
1. Retrieves the required pre-built singularity image 
2. Defines input/output parameters
3. Executes the processing script [`processing.R`](src/processing.R) within a singularity container

# Notes
Expect ~2 hours of runtime per sample with 20 cores

The output files of interest are
| Filename                  | Description |
| --------                  | ----------- |
| `<specfile>.quant`        | ?? |
| `report-lib.tsv`          | ?? |
| `report-lib.tsv.speclib`  | ?? |
| `report.pr_matrix.tsv`    | precursor ion quantities |
| `report.pg_matrix.tsv`    | protein group quantities |
| `report.gg_matrix.tsv`    | gene group quantities |
| `report.stats.tsv`        | gene group quantities |
| `report.tsv.tsv`          |  precursor ions identified, quantities, quality metrics and annotations |

---




## Running the processing R script
After ensuring the singularity image is available, [`run.sh`](src/run.sh) defines input/output
parameters and executes [`processing.R`](src/processing.R) within the container, generating
plots or tables within the defined output directories.

```bash
pro_input='input/spe_Report_Proteins.csv' #csv of the protein groups quntification generated by SN
project_name='test' #the name of this project
outdir='./output' #output dir
design_matrix='input/design_matrix.csv' #csv file of design matrix for comparison
r_script='processing.R'

# unused / NYI:
# pep_input='' #the csv file of the peptide intensity from the SN
# normalization='' # 'T' or 'F'

module load singularity

singularity run -H $PWD:/home src/${img}.sif \
    Rscript src/${r_script} \
    --pro_input ${pro_input} \
    -p ${project_name} \
    -o ${outdir} \
    --design_matrix ${design_matrix}
```




## Proteomics data anylysis 

This workflow performs the following tasks:
- [Spectronaut](https://biognosys.com/resources/spectronaut-the-deepest-proteome-coverage-available/) or [DIA-NN](https://github.com/vdemichev/DiaNN) for Protein Analysis and Quantification
- Quality control and data filltering
- Normalization and imputation(OPTIONAL)
- Data Cluster
- Differential expression analysis and visualization
- Downstream functional enrichment analysis
- Run all steps from start to finish

## System requirements
- Windows for Spectronaut
- Linux for DIA-NN
- R 



## DIA-NN for Protein Analysis and Quantification
We uses the DIA-NN software for protein analysis and quantification on biowulf or Windows systerm.

Use the following commands to run DIA-NN on prebuilt module on biowulf：
``` bash
# TODO: change diann to work as singularity container
module load diann
diann \
    --f ../20210208_KLOF_DIA_FAIMS_35V_d0_1.mzML   \
    --lib  \
    --threads 24 \
    --verbose 1 \
    --out ./report.tsv \
    --matrices \
    --out-lib ./report-lib.tsv \
    --gen-spec-lib \
    --predictor \
    --fasta ../uniprot-proteome_Human_UP000005640_20191105.fasta \
    --fasta-search \
    --min-fr-mz 200 \
    --max-fr-mz 2000 \
    --met-excision \
    --cut K*,R* \
    --missed-cleavages 2 \
    --min-pep-len 7 \
    --max-pep-len 52 \
    --min-pr-mz 300 \
    --max-pr-mz 1800 \
    --min-pr-charge 1 \
    --max-pr-charge 4 \
    --unimod4 \
    --var-mods 5 \
    --var-mod UniMod:35,15.994915,M \
    --var-mod UniMod:1,42.010565,*n \
    --monitor-mod UniMod:1 \
    --reanalyse \
    --relaxed-prot-inf \
    --smart-profiling \
    --peak-center \
    --no-ifs-removal  \
``` 

## Quality control and data filltering
We used the R to process the MS data, visualize the samples quality and fillter some sample with poor quality. This step will provide the figures including the total number of identified and quantified protein groups, the distribution of protein intensity, and correlation among protein abundance of the biological replicates.

Runing the code:
``` bash
# QC.R does not exist
Rscript src/QC.R \
    --pro_input $pro_input \
    --pep_input $pep_input \
    -p $project_name \
    -o $outdir
``` 

## Data Cluster
The data cluster include HC-cluster, PCA, and UMAP.

``` bash
# TODO: cluster_plot.R renamed cluster.R?
Rscript cluster_plot.R -i $pro_input -p $project_name -o $outdir
``` 


## Differential expression analysis(T-test) and downstream functional enrichment analysis
DE analysis is done using the t-test. Downstream functional enrichment analysis for Differential expression gene. It takes the MS quatification from Spectronaut as an input:
```
# TODO: DE_enrichment.R does not exist
Rscript DE_enrichment.R -i $pro_input -c $control -o $outdir 
```

## Bash comand line

```
# TODO: spe_pro.bash does not exist
bash spe_pro.bash
```


## Tiny example

### Generating tiny example files
Wanting to subset down to 5% of FASTA file and 1% of mzML file. Subsetting the fasta is easy.
`wc -l example/uniprot-proteome_Human_UP000005640_20191105.fasta` yields 539137 lines. First 4999
lines evenly splits between entries, and is ~1% of the file.

```bash
head -n 4999 example/uniprot-proteome_Human_UP000005640_20191105.fasta > tiny/tiny.fasta
```

Splitting the mzML file is more difficult due to XML formatting. One of the files 
`example/uniprot-proteome_Human_UP000005640_20191105.fasta` contains 78301 `<pectrum index>` fields.
We can subset by excluding the lines for 99% of the indices. If we only want ~700 of the indices,
we exclude indices 701-78301. Each field ends with a `</spectrum>` tag.

Based on pattern finding and line numbers, we want lines 1-(line BEFORE index 701, i.e. `<`) and from
(line CONTAINING `</spectrumList>`, i.e. `>=`) through EOF.

```bash
lower=$(grep -n '<spectrum index="701"' example/raw_MS_mzML/HREC_IRIS_1.mzML | cut -d ':' -f 1)
upper=$(grep -n '</spectrumList>' example/raw_MS_mzML/HREC_IRIS_1.mzML | cut -d ':' -f 1)
awk -v lower=${lower} \
    -v upper=${upper} \
    'NR < lower || NR >= upper' example/raw_MS_mzML/HREC_IRIS_1.mzML \
    > tiny/tiny.mxML
```

### Test run with container
```bash
src/tiny-diann-test.sh tiny/tiny.mxML tiny/tiny.fasta 
```

</details>
