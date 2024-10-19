# README

![workflow-image](src/ProtPipe.png)

This repository facilitates quick, simple, and reproducible access to Data Independent Acquisition (DIA) proteomics workflows with minimal command line experience.

We include a convenient wrapper script for running DIA-NN inside a pre-built singularity image to first estimate protein abundance from raw mass spec output. Protein abundance estimates (accepting estimates from both `DIA-NN` and `Spectronaut`) can  be processed in a preconfigured R environment, generating QC reports, and various analyses and visualizations.

The ProtPipe web application (http://34.42.19.73:8501/) provides a user-friendly, interactive interface for performing differential expression analysis with a single click. It is dedicated to downstream analysis following database searches.

# Quick-start

1. Ensure `singularity` is installed and accessible on your system. Many HPCs (including NIH Biowulf) come with this pre-installed as a module. If your HPC has singularity installed, it will be automatically detected and loaded when necessary.
2. Clone this repository, i.e. execute `git clone https://github.com/NIH-CARD/ProtPipe`
3. If you are predicting protein abundances from raw mass spec output, look over and edit any custom `DIA-NN` parameters inside [`config.txt`](config.txt). You can either edit `config.txt` directly (and it will be used by default), or make a copy and save it to a different file name, then reference it with `--cfg newfilename.txt` when running the wrapper script.


# Installing Singularity

This workflow requires that [`Singularity`](https://sylabs.io/singularity) be available, which runs natively on a Linux system. `Singularity` is containerization software that allows an entire pre-configured computing environment to be accessed--reducing installation headaches and improving reproducibility. 

*We highly recommend using a workstation or HPC with a native Linux installation.* Not only does this simplify the usage of `singularity`, it also would likely provide greater resources for DIA-NN's intensive computation.

To run on your personal/local non-Linux machine, Mac users need to first install a number of dependencies. Windows users would either need to use a virtual machine, or run things through the Windows Subsystem for Linux (WSL). Explaining the installation of `singularity` on these non-Linux systems is beyond the scope of this guide, so we defer to [the documentation here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

## Converting Mass Spec file formats

DIA-NN cannot handle some propietary file formats such as thermo fisher RAW. Thus these files must
be converted (i.e. to mzML) prior to running DIA-NN. Conversion can be done with the included script
[pwiz-convert.sh](src/pwiz-convert.sh). 

Conversion can be done by specifying either
* a single input file with `--file`
* an entire directory with `--dir`
* a text file that lists inputs, one per line, with `--list`
Along with a single output directory with `--out`.

For example:
```
bash src/pwiz-convert.sh --file myfile.raw --out mzml_outdir
bash src/pwiz-convert.sh --dir path/to/rawfiles/ --out mzml_outdir
bash src/pwiz-convert.sh --list rawfiles.txt --out mzml_outdir

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

# Predicting Protein Abundances (running DIA-NN)
You will need to manually build your singularity container for running DIA-NN, which requires agreement to their license terms (see their [GitHub Page](https://github.com/vdemichev/DiaNN/) for more info).

The singularity definition file is contained in [`diann.def`](src/diann-1.8.1.def). With this file, you can run `sudo singularity build diann.sif diann-1.8.1.def`. After the build completes, ensure the `sif` file is moved to 'src/diann-1.8.1.sif'.
 
After editing the contents of [`config.txt`](config.txt), or generating a new file to specify with `--cfg newfile.txt`:
```bash
# Submit to SLURM
sbatch src/diann.sh --cfg config.txt

# Run Locally
src/diann.sh --cfg config.txt
```

# Post-analysis

First you must retrieve the pre-built singularity image with the required `R` version and package dependencies. You can retrieve the image by executing:

```bash
singularity pull src/R.sif docker://quay.io/datatecnica/protpipe:latest
```

## Running containerized interactive `R` session
You can start an interactive R session within the container as follows:

```bash
./protpipe.sh R
```

The above is shorthand for executing the followiing:
```bash
singularity exec -B ${PWD} R
```

where `singularity exec R` starts the `R` session, while `-B ${PWD}` binds the
current directory within the container. Without binding, the current directory's
files would not be visible inside the container.

## Basic MS data analysis: `./protpipe.sh basic`

For performing QC and running differential abundance or enrichment analysis for typical mass spec data. The required inputs are
- protein intensity estimates from DIA-NN or Spectronaut
- experimental design matrix csv file

```bash
./protpipe.sh basic \
    --pgfile EXAMPLES/DIFF_ABUNDANCE/iPSC.csv \
    --design EXAMPLES/DIFF_ABUNDANCE/design_matrix_iPSC.csv \
    --out EXAMPLES/DIFF_ABUNDANCE/
```


## Affinity Purification Mass Spec analysis: `./protpipe.sh APMS`

Similar to `basic` but for affinity purification mass spec. Requires the user to specify which protein was used for pulldown (`--ip`)
```bash
./protpipe.sh APMS \
    --pgfile EXAMPLES/APMS/APMS.csv \
    --design EXAMPLES/APMS/design_matrix_APMS.csv \
    --ip UNC13A \
    --out EXAMPLES/APMS/
```

## Immunopeptidome analysis: `./protpipe.sh immuno`

Requires the csv or tsv output from `FragPipe` and a `csv` specifying HLA typing.
```bash
./protpipe.sh immuno \
    --pepfile EXAMPLES/IMMUNOPEPTIDOME/combined_peptide.tsv \
    --out EXAMPLES/IMMUNOPEPTIDOME/ \
    --hla EXAMPLES/IMMUNOPEPTIDOME/HLA_typing.csv
```
## Parametrization and command line options
```
--pgfile  	Input file of Protein Group Intensity (from DIA-NN or Spectronaut). *Required*
          	eg: --pgfile data/protein_groups.tsv
--design  	Comma- or tab-delimited file specifying the experimental design. 	*Required*	
          	eg:--design design/experiment_design.tsv
--ip	 	Protein name of the IP.	*Required for AMPS*		 
		eg:--ip UNC13A
--hla	  	The HLA typing information.	*Required for Immunopeptidome*			 	
		eg:--hla HLA_typing.csv
--base   	Base for log transformation of intensity data
          	default:2
          	eg:--base 2
--normalize	Method to normalize sample intensities ('shift', 'scale', 'none').
            	default:none
		eg:--normalize shift
--exclude	Semicolon-separated string of files to exclude from analysis.
		eg: --exclude sample1_name
--sds		Filter out samples with protein group counts > N standard deviations from the mean.
		default:3
		eg: --sds 3
--minintensity	Minimum linear (not log) intensity.
		default:0
		eg: --minintensity 500
--fdr		False Discovery Rate threshold for differential abundance analysis.
	 	default:0.01
		eg:--fdr 0.01
--foldchange	Minimum linear fold change for labeling protein groups in differential abundance analysis.
	 	default:2
		eg:--foldchange 2
--enrich	Cutoff p-value for gene enrichment analysis.
		default:0.01
		eg:--enrich 0.01
--gsea		Cutoff False Discovery Rate for GSEA analysis.
	 	default:0.01
		eg:--gsea 0.01

```


