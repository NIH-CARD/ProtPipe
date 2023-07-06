# README

![workflow-image](src/ProtPipe.png)

This repository facilitates quick, simple, and reproducible access to Data Independent Acquisition (DIA) proteomics workflows with minimal command line experience.

We include a convenient wrapper script for running DIA-NN inside a pre-built singularity image to first estimate protein abundance from raw mass spec output. Protein abundance estimates (accepting estimates from both `DIA-NN` and `Spectronaut`) can  be processed in a preconfigured R environment, generating QC reports, and various analyses and visualizations.


# Quick-start

1. Ensure `singularity` is installed and accessible on your system. Many HPCs (including NIH Biowulf) come with this pre-installed as a module. If your HPC has singularity installed, it will be automatically detected and loaded when necessary.
2. Clone this repository, i.e. execute `git clone https://github.com/cory-weller/ProtPipe.git`
3. If you are predicting protein abundances from raw mass spec output, look over and edit any custom `DIA-NN` parameters inside [`config.txt`](config.txt). You can either edit `config.txt` directly (and it will be used by default), or make a copy and save it to a different file name, then reference it with `--cfg newfilename.txt` when running the wrapper script.


# Installing Singularity

This workflow requires that [`Singularity`](https://sylabs.io/singularity) be available, which runs natively on a Linux system. `Singularity` is containerization software that allows an entire pre-configured computing environment to be accessed--reducing installation headaches and improving reproducibility. 

*We highly recommend making use a workstation or HPC with a native Linux installation.* Not only does this simplify the usage of `singularity`, it also would likely provide greater resources for DIA-NN's intensive computation.

To run on your personal/local non-Linux machine, Mac users need to first install a number of dependencies. Windows users would either need to use a virtual machine, or run things through the Windows Subsystem for Linux (WSL). Explaining the installation of `singularity` on these non-Linux systems is beyond the scope of this guide, so we defer to [the documentation here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

# Predicting Protein Abundances (running DIA-NN)
After editing the contents of [`config.txt`](config.txt), or generating a new file to specify with `--cfg newfile.txt`:
```bash
# Submit to SLURM
sbatch src/diann.sh --cfg config.txt

# Run Locally
src/diann.sh --cfg config.txt
```

# Post-analysis
## Processing total protein intensity estimates
The required documents are the csv or tsv file from the Spectronaut, design matrix csv file(example shown in folder example/design_matrix.csv).
```bash
src/analyze.sh --pgfile TEST/report.pg_matrix.tsv --design TEST/design.tsv --out TEST/
```

<details><summary>Samples for iPSCs neuron differentiation</summary>

```bash
# WITH differentiation neuron samples, Day0 as control
Rscript src/counts_processing.R --pgfile iPSC_neuron/luke.csv --design iPSC_neuron/design_matrix_iPSC_neuron.csv --out iPSC_neuron/
```

</details>

## Processing AP-MS data analysis
The required documents are the csv or tsv file from the Spectronaut, design matrix csv file(example shown in folder example/design_matrix.csv) and the gene name for the pulling down.
```bash
Rscript src/APMS.R --pgfile APMS/apms.csv --design APMS/design_matrix_apms.csv --out APMS/ --ip UNC13A
```

## Processing Immunopeptidome data analysis
The required documents are the csv or tsv file from the FragPip, HLA_typing csv file(example shown in folder example/HLA_typing.csv).

```bash
Rscript src/Peptidome.R --pepfile peptidome/combined_peptide.tsv  --out peptidome/ --hla peptidome/HLA_typing.csv
```


## Converting Mass Spec file formats

DIA-NN cannot handle some propietary file formats such as thermo fisher RAW. Thus these files must
be converted (i.e. to mzML) prior to running DIA-NN. Conversion can be done interactively or by
submitting to your HPC with `sbatch`.

Converting a single file:
```bash
src/pwiz-convert /path/to/your/MassSpecFile.raw
```

Mass spec file conversion is handled by ProteoWizard (via wine in a singularity container).
A writable sandboxed version of the container (which is required to run ProteoWizard) was built
and modified from a [docker image](docker://chambm/pwiz-skyline-i-agree-to-the-vendor-licenses) on
March 02 2023. Steps were modified from [here](https://github.com/jspaezp/elfragmentador-data#setting-up-msconvert-on-singularity-).

<details><summary>Building pwiz container</summary>


```bash
# Build writable singularity sandbox image based on docker image
singularity build --sandbox pwiz_sandbox docker://chambm/pwiz-skyline-i-agree-to-the-vendor-licenses

# Modified pwiz_sandbox/usr/bin/mywine
echo """#!/bin/sh

GLOBALWINEPREFIX=/wineprefix64
MYWINEPREFIX=/mywineprefix/

if [ ! -L "$MYWINEPREFIX"/dosdevices/z: ] ; then 
  mkdir -p "$MYWINEPREFIX"/dosdevices
  cp "$GLOBALWINEPREFIX"/*.reg "$MYWINEPREFIX"
  ln -sf "$GLOBALWINEPREFIX/drive_c" "$MYWINEPREFIX/dosdevices/c:"
  ln -sf "/" "$MYWINEPREFIX/dosdevices/z:"
  echo disable > $MYWINEPREFIX/.update-timestamp        # Line being added
  echo disable > $GLOBALWINEPREFIX/.update-timestamp    # Line being added
fi 

export WINEPREFIX=$MYWINEPREFIX
wine "$@"
""" > pwiz_sandbox/usr/bin/mywine

tar -czvf pwiz_sandbox.tar.gz pwiz_sandbox
rclone copy pwiz_sandbox.tar.gz onedrive:/singularity       # upload archive to cloud
```
</details>

# Bulk convert raw to mzml
```bash
filelist='rawfiles.txt'
cat <(find ANXA11_redux | grep raw) > ${filelist}
nfiles=$(wc -l ${filelist} | awk '{print $1}')
sbatch --array=1-${nfiles} src/pwiz-convert-array.sh ${filelist}
```
