# README

## DIA-NN singularity image

The image was built from the recipe file [`diann-1.8.1.def`](diann-1.8.1.def) using the sylabs
remote builder (via web interface).

The md5 sum should be `35644c1d7217f0c65727b8fb9c8bfaae`. The image is automatically pulled and
checked when executing [`run-diann.sh`](run-diann.sh).


Note: The container must be called using the `--cleanenv` option, otherwise the container may fail 
due to collision between system and container libraries.