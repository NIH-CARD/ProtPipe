#!/usr/bin/env bash

#### CHECK SINGULARITY #############################################################################
if ! command -v singularity &> /dev/null; then
    echo "INFO: singularity command not found, looking for singularity module"
    if ! command -v module &> /dev/null; then
        echo "ERROR: module command not found. Did you mean to run this on an HPC?"
        exit 1
    else
        if $(module avail singularity/3 2>&1 >/dev/null | grep -q 'No module'); then
            echo 'ERROR: singularity cannot be found. Recheck installation?'
            exit 1
        else
            echo 'INFO: module singularity found'
            module load singularity/3
        fi
    fi
else
    echo 'INFO: singularity command found'
fi



#### R ANALYSIS ####################################################################################
r_sif="src/R.sif"
r_sif_sha256_desired='cfab1ee7f61e2af5dff7b832ce28768ce5df2ab949c482a5bd94a91383423bb5'
if [ ! -f "${r_sif}" ]; then
    echo "INFO: pulling image from Singularity cloud"
    singularity pull --arch amd64 library://wellerca/r/4.0:sha256.cfab1ee7f61e2af5dff7b832ce28768ce5df2ab949c482a5bd94a91383423bb5
    mv 4.0_sha256.cfab1ee7f61e2af5dff7b832ce28768ce5df2ab949c482a5bd94a91383423bb5.sif "${r_sif}"
else
    echo "INFO: ${r_sif} already exists, skipping download"
fi

r_sif_sha256_actual=$(sha256sum "${r_sif}" | awk '{print $1}')

if [ ! "${r_sif_sha256_actual}" == "${r_sif_sha256_desired}" ]; then
    echo "ERROR: ${r_sif} sha256sum does not pass check. Possibly corrupted? Delete or clear singularity cache and try again."
    exit 1
else
    echo "INFO: ${r_sif} sha256sum sum passes check"
fi

r_processing_script='src/counts_processing.R'
singularity exec --cleanenv -H ${PWD} ${r_sif} Rscript ${r_processing_script} $@
