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


BADARGS='FALSE'

if [ ! -n "${DRYRUN}" ]; then
    DRYRUN='FALSE'
fi


print_preamble() {
    echo -e '\n\nProtPipe: An all-in-one wrapper for DIA analysis in a singularity container'
    echo -e 'See GitHub for updates: https://www.github.com/cory-weller/ProtPipe\n'
    echo -e 'Using DIA-NN version: 1.8.1 ran in singularity\n'
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
    echo -e '--config\t\tImport from a file other than the default config.txt'
}


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

check_dry_run() {
    if [ "${DRYRUN}" == 'TRUE' ]; then
        echo 'INFO: doing a --dry-run'
    fi
}

check_singularity_exists() {
    # Check for singularity
    if ! command -v singularity &> /dev/null; then
        echo "INFO: singularity command not found, looking for singularity module"
        if ! command -v module &> /dev/null; then
            echo "WARNING: module command not found. Did you mean to run this on an HPC?"
            BADARGS='TRUE'
        else
            if $(module avail singularity 2>&1 >/dev/null | grep -q 'No module'); then
                echo 'WARNING: module singularity not found'
                echo 'ERROR: singularity cannot be found. Recheck installation?'
                BADARGS='TRUE'
            else
                echo 'INFO: module singularity found'
                SINGULARITY_EXISTS='TRUE'
                module load singularity
            fi
        fi
    else
        echo 'INFO: singularity command found'
        SINGULARITY_EXISTS='TRUE'

    fi
}

check_singularity_version() {
    if [ "${SINGULARITY_EXISTS}" == 'TRUE' ]; then
        if ! singularity --version | grep -q "version 4."; then
            echo "WARNING: Singularity is not version 4. This may result in unexpected behavior."
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

check_clobber() {
    if [ "${CLOBBER}" == 'TRUE' ]; then
        echo -e 'INFO: ignoring and re-generating pre-existing files due to --clobber flag\n'
        echo -e 'WARNING: --clobber has been specified. Pre-existing files will be overwritten.'
        echo -e 'WARNING: Starting in 10 seconds unless interrupted.'
        for i in $(seq 10 -1 1); do echo "...$i" && sleep 1; done
        unset DIANN_USE_QUANT
    fi
}

test_and_run() {
    echo -e 'INFO: validating arguments and starting DIA-NN\n'
    python3 src/validate-args.py ${CONFIG} ${SINGULARITY_IMAGE} ${DRYRUN}
}

# stop_if_dryrun() {
#     if [ "${DRYRUN}" == 'TRUE' ]; then
#         echo -e 'INFO: halting process before DIA-NN due to --dry-run flag\n'
#         echo -e 'Shutting down...\n'
#         exit 0
#     fi
# }




### RUN ############################################################################################


# Initialization
print_preamble                  # Print welcome banner
check_help                      # Print help message if --help supplied

# Check provided inputs
check_dry_run                   # if --dry-run provided, stop before analysis
check_config_input              # Use config.txt if unspecified; confirm file exists/can be read
check_singularity_exists
check_singularity_version
check_singularity_image         # Pull .sif if necessary; confirm md5sum

# Perform DIA-NN steps
test_and_run                    # validate args and run DIA-NN if checks pass; stop if dryrun

exit 0