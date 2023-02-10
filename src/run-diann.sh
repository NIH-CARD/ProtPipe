#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --tasks-per-node 20
#SBATCH --mem 80G
#SBATCH --time 0-4:00:00
#SBATCH --partition norm

trap '[[ $? -eq 1 ]] && echo Halting execution due to errors' EXIT

SRC_DIR='./src'
SINGULARITY_IMAGE="${SRC_DIR}/diann-1.8.1.sif"
ARGS="$@"

### FUNCTIONS ######################################################################################

# Argument parsing
usage_error () { echo >&2 "$(basename $0):  $1"; exit 2; }
assert_argument () { test "$1" != "$EOL" || usage_error "$2 requires an argument"; }
if [ "$#" != 0 ]; then
    EOL=$(printf '\1\3\3\7')
    set -- "$@" "$EOL"
    while [ "$1" != "$EOL" ]; do
        opt="$1"; shift
        case "$opt" in

            # Your options go here.
            --dry-run) DRYRUN='TRUE';;
            --help) HELP='TRUE';;
            --clobber) CLOBBER='TRUE';;
            --fasta) assert_argument "$1" "$opt"; FASTA_INPUT="$1"; shift;;
            --mzml) assert_argument "$1" "$opt"; MZML="$1"; shift;;
            --raw) assert_argument "$1" "$opt"; RAW="$1"; shift;;
            --dia) assert_argument "$1" "$opt"; DIA="$1"; shift;;
            --out) assert_argument "$1" "$opt"; OUTPUT_DIR="$1"; shift;;
            --config) assert_argument "$1" "$opt"; CONFIG="$1"; shift;;
      
            # Arguments processing. You may remove any unneeded line after the 1st.
            -|''|[!-]*) set -- "$@" "$opt";;                                          # positional argument, rotate to the end
            --*=*)      set -- "${opt%%=*}" "${opt#*=}" "$@";;                        # convert '--name=arg' to '--name' 'arg'
            -[!-]?*)    set -- $(echo "${opt#-}" | sed 's/\(.\)/ -\1/g') "$@";;       # convert '-abc' to '-a' '-b' '-c'
            --)         while [ "$1" != "$EOL" ]; do set -- "$@" "$1"; shift; done;;  # process remaining arguments as positional
            -*)         usage_error "unknown option: '$opt'";;                        # catch misspelled options
            *)          usage_error "this should NEVER happen ($opt)";;               # sanity test for previous patterns
    
        esac
    done
    shift  # $EOL
fi

initialize_args() {
    # Sets default value for BADARGS to FALSE, which is flipped to TRUE in the event future args
    # fail validation checks. Also counts the number of fasta or mass spec inputs provided within
    # initial submission command
    BADARGS='FALSE'
    # N_fasta=$(echo ${ARGS} | grep -o  -- '--fasta' | wc -l)
    # N_mzml=$(echo ${ARGS} | grep -o  -- '--mzml'  | wc -l)
    # N_raw=$(echo ${ARGS} | grep -o  -- '--raw'   | wc -l)
    # N_dia=$(echo ${ARGS} | grep -o  -- '--dia'   | wc -l)
}

print_preamble() {
    echo -e '\n\nProtPipe: An all-in-one wrapper for DIA analysis in a singularity container'
    echo -e 'See GitHub for updates: https://www.github.com/cory-weller/ProtPipe\n'
    echo -e 'Using DIA-NN version: 1.8.1 ran in singularity 3.x\n'
    echo -e "command called:\n$(echo run-diann.sh $ARGS | sed -e 's/\( --\)/ \\\n   \1/g')\n"
}

check_help() {
    if [ "${HELP}" == 'TRUE' ]; then
        print_helpmsg
        echo -e '\nexiting due to --help flag\n\n'
        exit 0
    fi
}

print_helpmsg() {
    echo -e 'Option\t\t\tDescription'
    echo -e '------\t\t\t-----------'
    echo -e '--help\t\t\tPrint this message'
    echo -e '--dry-run\t\tCheck validity of inputs, but do not run analysis or write output'
    echo -e '--fasta\t\t\tPeptide library fasta file (required)'
    echo -e '--dia,--raw,--mzml\tMass spec input to analyze (exactly one required)'
    echo -e '--out\t\t\tOutput directory (default: current directory)'
    echo -e '--config\t\tImport from a file other than the default config.txt'
    echo -e '--clobber\t\tIgnore existing files, regenerate and overrwite if necessary\n'
}

# check_mass_spec_input() {
#     let N_SPECFILES=${N_mzml}+${N_raw}+${N_dia} 
#     if [ "${N_SPECFILES}" -gt 1 ]; then     # If more than one mass spec input is given
#             echo 'ERROR: Multiple mass spec inputs provided. Provide ONE  of --mzml --dia or --raw'
#             BADARGS='TRUE'
#     else
#         if [ -z "${MZML}" ] && [ -z "${RAW}" ] && [ -z "${DIA}" ]; then
#             echo "ERROR: --mzml --raw or --dai (mass spec input) is required"
#             BADARGS='TRUE'
#         else
#             ## Only one mass spec is allowed at this point; concatenation of vars is identical to
#             ## Selecting the one provided by user
#             MASS_SPEC_INPUT="${MZML}${DIA}${RAW}"
#             ## If Mass Spec file is defined but does not exist/can't be read
#             if [ ! -r "${MASS_SPEC_INPUT}" ]; then
#                 echo "ERROR: Mass spec input ${MASS_SPEC_INPUT} does not exist or is not readable"
#                 BADARGS='TRUE'
#             fi
#         fi
#     fi
# }

check_config_input() {
    if [ -z "${CONFIG}" ]; then
        echo "INFO: --config not specified, config.txt will be imported by default"
        CONFIG='config.txt'
    fi
    if [ ! -r "${CONFIG}" ]; then
        echo "ERROR: config file ${CONFIG} cannot be read or does not exist"
        BADARGS='TRUE'
    fi
}

# check_output_dir() {
#     if [ -z "${OUTPUT_DIR}" ]; then
#         echo "INFO: --out not specified, output will be saved to current working directory"
#         OUTPUT_DIR=${PWD}
#     fi
#     if mkdir -p "$OUTPUT_DIR" &> /dev/null; then
#         echo "INFO: output directory ${OUTPUT_DIR} is writable"
#     else
#         echo "ERROR: could not create or write to output directory ${OUTPUT_DIR}"
#         BADARGS='TRUE'
#     fi
# }

check_dry_run() {
    if [ "${DRYRUN}" == 'TRUE' ]; then
        echo 'INFO: doing a --dry-run'
    fi
}

# check_fasta_input() {
#     if [ -z "${FASTA_INPUT}" ]; then 
#         echo "ERROR: --fasta <input.fa> is required"
#         BADARGS='TRUE'
#     else
#         if [ "${N_fasta}" -gt 1 ]; then
#             echo "WARNING: multiple --fasta provided. Only proceeding with the last one, ${FASTA_INPUT}"
#         fi
#         ## If FASTA_INPUT file is defined but does not exist/can't be read
#         if [ ! -r "${FASTA_INPUT}" ]; then
#             echo "ERROR: FASTA input ${FASTA_INPUT} does not exist or is not readable"
#             BADARGS='TRUE'
#         fi
#     fi
# }

check_singularity_exists() {
    # Check for singularity
    if ! command -v singularity &> /dev/null; then
        echo "INFO: singularity command not found, looking for singularity module"
        if ! command -v module &> /dev/null; then
            echo "WARNING: module command not found. Did you mean to run this on an HPC?"
            BADARGS='TRUE'
        else
            if $(module avail singularity/3 2>&1 >/dev/null | grep -q 'No module'); then
                echo 'WARNING: module singularity not found'
                echo 'ERROR: singularity cannot be found. Recheck installation?'
                BADARGS='TRUE'
            else
                echo 'INFO: module singularity found'
                SINGULARITY_EXISTS='TRUE'
                module load singularity/3
            fi
        fi
    else
        echo 'INFO: singularity command found'
        SINGULARITY_EXISTS='TRUE'

    fi
}

check_singularity_version() {
    if [ "${SINGULARITY_EXISTS}" == 'TRUE' ]; then
        if ! singularity --version | grep -q "version 3."; then
            echo "WARNING: Singularity is not version 3. This may result in unexpected behavior."
        fi
    fi
}

check_singularity_image() {
    md5_desired='35644c1d7217f0c65727b8fb9c8bfaae'
    if [ ! -f "${SINGULARITY_IMAGE}" ]; then
        singularity pull \
            --arch amd64 \
            --name ${SINGULARITY_IMAGE} \
            library://wellerca/diann/1.8.1:0.9
    fi

    if [ ! -f "${SINGULARITY_IMAGE}" ]; then
        echo -e "ERROR: ${SINGULARITY_IMAGE} does not exist. This should not happen unless the pull command failed"
        exit 1
    fi

    md5_actual=$(md5sum ${SINGULARITY_IMAGE} | awk '{print $1}')

    if [ ! "${md5_actual}" == "${md5_desired}" ]; then
        echo -e "ERROR: singularity image ${SINGULARITY_MAGE} is be incomplete or corrupt. Try removing it and trying again."
        exit 1
    fi
}

# check_arguments_valid() {
#     if [ "${BADARGS}" == 'TRUE' ]; then
#         echo -e '\nCheck arguments and try again.\n'
#         print_helpmsg
#         exit 1
#     else
#         echo -e '\nSUCCESS: arguments passed validation'
#         echo -e '\n###########\nPARAMETERS:\n###########\n'
#         echo -e "Config file:\t${CONFIG}"
#         echo -e "Spec Input:\t${MASS_SPEC_INPUT}"
#         echo -e "FASTA Input:\t${FASTA_INPUT}"
#         echo -e "Output to:\t${OUTPUT_DIR}/"
#         echo -e "Singularity:\t${SINGULARITY_IMAGE}\n"
#     fi
# }

check_clobber() {
    if [ "${CLOBBER}" == 'TRUE' ]; then
        echo -e 'INFO: ignoring and re-generating pre-existing files due to --clobber flag\n'
        echo -e 'WARNING: --clobber has been specified. Pre-existing files will be overwritten.'
        echo -e 'WARNING: Starting in 10 seconds unless interrupted.'
        for i in $(seq 10 -1 1); do echo "...$i" && sleep 1; done
        unset DIANN_USE_QUANT
    fi
}

# import_config() {
#     . ${CONFIG}
#     case "${DIANN_REANALYSE}" in
#         TRUE)  DIANN_REANALYSE='--reanalyse' ;;
#         FALSE) unset DIANN_REANALYSE ;;
#         *)     echo 'ERROR: DIANN_REANALYSE must be TRUE or FALSE' && BADCONFIG='TRUE' ;;
#     esac

#     case "${DIANN_RELAXED_PROT_INF}" in
#         TRUE)  DIANN_RELAXED_PROT_INF='--relaxed-prot-inf' ;;
#         FALSE) unset DIANN_RELAXED_PROT_INF ;;
#         *)     echo 'ERROR: DIANN_RELAXED_PROT_INF must be TRUE or FALSE' && BADCONFIG='TRUE' ;;
#     esac

#     case "${DIANN_SMART_PROFILING}" in
#         TRUE)  DIANN_SMART_PROFILING='--smart-profiling' ;;
#         FALSE) unset DIANN_SMART_PROFILING ;;
#         *)     echo 'ERROR: DIANN_SMART_PROFILING must be TRUE or FALSE' && BADCONFIG='TRUE' ;;
#     esac

#     case "${DIANN_PEAK_CENTER}" in
#         TRUE)  DIANN_PEAK_CENTER='--peak-center' ;;
#         FALSE) unset DIANN_PEAK_CENTER ;;
#         *)     echo 'ERROR: DIANN_PEAK_CENTER must be TRUE or FALSE' && BADCONFIG='TRUE' ;;
#     esac

#     case "${DIANN_NO_IFS_REMOVAL}" in
#         TRUE)  DIANN_NO_IFS_REMOVAL='--no-ifs-removal' ;;
#         FALSE) unset DIANN_NO_IFS_REMOVAL ;;
#         *)     echo 'ERROR: DIANN_NO_IFS_REMOVAL must be TRUE or FALSE' && BADCONFIG='TRUE' ;;
#     esac

#     case "${DIANN_MET_EXCISION}" in
#         TRUE)  DIANN_MET_EXCISION='--met-excision' ;;
#         FALSE) unset DIANN_MET_EXCISION ;;
#         *)     echo 'ERROR: DIANN_MET_EXCISION must be TRUE or FALSE' && BADCONFIG='TRUE' ;;
#     esac

#     case "${DIANN_FASTA_SEARCH}" in
#         TRUE)  DIANN_FASTA_SEARCH='--fasta-search' ;;
#         FALSE) unset DIANN_FASTA_SEARCH ;;
#         *)     echo 'ERROR: DIANN_FASTA_SEARCH must be TRUE or FALSE' && BADCONFIG='TRUE' ;;
#     esac

#     case "${DIANN_USE_QUANT}" in
#         TRUE)  DIANN_USE_QUANT='--use-quant' ;;
#         FALSE) unset DIANN_USE_QUANT ;;
#         *)     echo 'ERROR: DIANN_USE_QUANT must be TRUE or FALSE' && BADCONFIG='TRUE' ;;
#     esac

#     case "${DIANN_GEN_SPEC_LIB}" in
#         TRUE)  DIANN_GEN_SPEC_LIB='--gen-spec-lib' ;;
#         FALSE) unset DIANN_GEN_SPEC_LIB ;;
#         *)     echo 'ERROR: DIANN_GEN_SPEC_LIB must be TRUE or FALSE' && BADCONFIG='TRUE' ;;
#     esac

#     case "${DIANN_PREDICTOR}" in
#         TRUE)  DIANN_PREDICTOR='--predictor' ;;
#         FALSE) unset DIANN_PREDICTOR ;;
#         *)     echo 'ERROR: DIANN_PREDICTOR must be TRUE or FALSE' && BADCONFIG='TRUE' ;;
#     esac

#     case "${DIANN_MATRICES}" in
#         TRUE)  DIANN_MATRICES='--matrices' ;;
#         FALSE) unset DIANN_MATRICES ;;
#         *)     echo 'ERROR: DIANN_MATRICES must be TRUE or FALSE' && BADCONFIG='TRUE' ;;
#     esac

#     case "${DIANN_REANNOTATE}" in
#         TRUE)  DIANN_REANNOTATE='--reannotate' ;;
#         FALSE) unset DIANN_REANNOTATE ;;
#         *)     echo 'ERROR: DIANN_REANNOTATE must be TRUE or FALSE' && BADCONFIG='TRUE' ;;
#     esac

#     if [ "${BADCONFIG}" == 'TRUE' ]; then
#         echo 'Check configuration and try again.'
#         exit 1
#     fi
# }


# format_args() {
# DIANN_ARGS="\
# ${DIANN_THREADS/#/--threads } \
# ${DIANN_RAW_MS_FILES[@]/#/--f } \
# ${DIANN_FASTAS[@]/#/--fasta } \
# ${OUTPUT_DIR/#/--out }/${DIANN_OUT} \
# ${DIANN_OUT_LIB/#/--out-lib } \
# ${DIANN_QVALUE/#/--qvalue } \
# ${DIANN_MIN_FR_MZ/#/--min-fr-mz } \
# ${DIANN_MAX_FR_MZ/#/--max-fr-mz } \
# ${DIANN_LIB/#/--lib } \
# ${DIANN_CUT/#/--cut } \
# ${DIANN_MISSED_CLEAVAGES/#/--missed-cleavages } \
# ${DIANN_MIN_PEP_LEN/#/--min-pep-len } \
# ${DIANN_MAX_PEP_LEN/#/--max-pep-len } \
# ${DIANN_MIN_PR_MZ/#/--min-pr-mz } \
# ${DIANN_MAX_PR_MZ/#/--max-pr-mz } \
# ${DIANN_MIN_PR_CHARGE/#/--min-pr-charge } \
# ${DIANN_MAX_PR_CHARGE/#/--max-pr-charge } \
# ${DIANN_VAR_MODS/#/--var-mods } \
# ${DIANN_MONITOR_MOD/#/--monitor-mod } \
# ${DIANN_VAR_MOD_LIST[@]/#/--var-mod } \
# ${DIANN_REANALYSE} \
# ${DIANN_RELAXED_PROT_INF} \
# ${DIANN_SMART_PROFILING} \
# ${DIANN_PEAK_CENTER} \
# ${DIANN_NO_IFS_REMOVAL} \
# ${DIANN_MET_EXCISION} \
# ${DIANN_USE_QUANT} \
# ${DIANN_FASTA_SEARCH} \
# ${DIANN_GEN_SPEC_LIB} \
# ${DIANN_PREDICTOR} \
# ${DIANN_MATRICES} \
# "
# }

# check_in_silico_lib() {
#     if [ ! -f "${OUTPUT_DIR}/report-lib.predicted.speclib" ]; then 
#         build_in_silico_lib && echo -e 'INFO: finished building in silico spectral library\n'
#     elif [ "${CLOBBER}" == 'TRUE' ]; then
#         echo -e 'WARNING: re-building in silico spectral library due to --clobber flag.'
#         echo -e 'WARNING: This will over-write the existing file.'
#         echo -e 'WARNING: Starting in 15 seconds unless interrupted.'
#         sleep 16
#         build_in_silico_lib && echo -e 'INFO: finished building in silico spectral library\n'
#     else
#         echo -e 'INFO: in silico spectral library already exists'
#         echo -e 'INFO: rerun with --clobber to delete and re-generate existing files\n'
#     fi
# }

# build_in_silico_lib() {
# echo -e 'INFO: starting generation of in silico spectral library\n'
# echo -e "calling command:\n singularity exec --cleanenv -H ${PWD} ${SINGULARITY_IMAGE} diann --fasta ${FASTA_INPUT} ${DIANN_ARGS}"
# singularity exec --cleanenv -H ${PWD} ${SINGULARITY_IMAGE} \
#     diann \
#     --fasta ${FASTA_INPUT} \
#     ${DIANN_ARGS}
# }

test_and_run() {
    echo -e 'INFO: validating arguments and starting DIA-NN\n'
    python3 src/validate-args.py ${CONFIG} ${SINGULARITY_IMAGE}
}

# run_diann() {
#     echo -e 'INFO: starting DIA-NN\n'
#     echo -e "calling command:\n singularity exec --cleanenv -H ${PWD} ${SINGULARITY_IMAGE} diann ${DIANN_ARGS}"
#     singularity exec --cleanenv -H ${PWD} ${SINGULARITY_IMAGE} \
#         diann ${DIANN_ARGS}
# }


# check_spec_sample() {
#     if [ ! -f "${OUTPUT_DIR}/report-lib.tsv" ]; then 
#         analyze_spec_sample 
#     elif [ "${CLOBBER}" == 'TRUE' ]; then
#         echo -e 'WARNING: re-analyzing mass spec sample due to --clobber flag.'
#         echo -e 'WARNING: This will over-write the existing file.'
#         echo -e 'WARNING: Starting in 15 seconds unless interrupted.'
#         sleep 16
#         analyze_spec_sample && echo -e 'INFO: finished building in silico spectral library\n'
#     else
#         echo -e 'INFO: mass spec sample already analyzed'
#         echo -e 'INFO: rerun with --clobber to delete and re-generate existing files\n'

#     fi
# }

# analyze_spec_sample() {
# echo -e 'INFO: starting analyzing mass spec sample\n'
# echo -e "calling command:\n singularity exec --cleanenv -H ${PWD} ${SINGULARITY_IMAGE} diann --f ${MASS_SPEC_INPUT} --fasta ${FASTA_INPUT} ${DIANN_ARGS}"
# singularity exec --cleanenv -H ${PWD} ${SINGULARITY_IMAGE} \
#     diann \
#     --f ${MASS_SPEC_INPUT} \
#     --fasta ${FASTA_INPUT} \
#     ${DIANN_ARGS}
# }

stop_if_dryrun() {
    if [ "${DRYRUN}" == 'TRUE' ]; then
        echo -e 'INFO: halting process before DIA-NN due to --dry-run flag\n'
        echo -e 'Shutting down...\n'
        exit 0
    fi
}

# print_config() {
#     echo -e "Imported configuration from ${CONFIG}:\n\n$DIANN_ARGS\n" | sed 's/ --/ \\\n--/g'
# }




### RUN ############################################################################################


# Initialization
print_preamble                  # Print welcome banner
check_help                      # Print help message if --help supplied
initialize_args                 # Set starting values

# Check provided inputs
check_dry_run                   # if --dry-run provided, stop before analysis
#check_fasta_input               # Confirm only 1 fasta input is given, and can be read
#check_mass_spec_input           # Confirm only 1 mass spec input is given, and can be read
check_config_input              # Use config.txt if unspecified; confirm file exists/can be read
# check_output_dir                # Use $PWD if unspecified; confirm write privilege
check_singularity_exists
check_singularity_version
check_singularity_image         # Pull .sif if necessary; confirm md5sum
# check_arguments_valid           # Print help and exit if any args invalid to this point

# Import run-specific configuration and prepare to start DIA-NN
# import_config                   # If --config not provided, defaults to config.txt
# format_args                     # Formats provided configuration as single-line args for DIA-NN
check_clobber                   # Must be before print_config to unset USE_QUANT if clobber=TRUE
# print_config                    # Print validated DIA-NN configuration to terminal
stop_if_dryrun                  # Quit before proceeding if --dry-run provided

# Perform DIA-NN steps
test_and_run                       # Build in silico library and analyze input

exit 0