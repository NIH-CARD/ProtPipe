#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --tasks-per-node 20
#SBATCH --mem 50G
#SBATCH --time 0-2:00:00
#SBATCH --partition quick

trap '[[ $? -eq 1 ]] && echo Halting execution due to errors' EXIT

SRC_DIR='./src'
SINGULARITY_IMAGE="${SRC_DIR}/diann-1.8.1.sif"

ARGS="$@"

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
            --fasta) assert_argument "$1" "$opt"; FASTA_INPUT="$1"; shift;;
            --mzml) assert_argument "$1" "$opt"; MZML="$1"; shift;;
            --raw) assert_argument "$1" "$opt"; RAW="$1"; shift;;
            --dia) assert_argument "$1" "$opt"; DIA="$1"; shift;;
            --out) assert_argument "$1" "$opt"; OUTPUT_DIR="$1"; shift;;
      
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


# Argument validation
BADARGS='FALSE' # initialize
N_fasta=$(echo ${ARGS} | grep -o  -- '--fasta' | wc -l)
N_mzml=$(echo ${ARGS} | grep -o  -- '--mzml'  | wc -l)
N_raw=$(echo ${ARGS} | grep -o  -- '--raw'   | wc -l)
N_dia=$(echo ${ARGS} | grep -o  -- '--dia'   | wc -l)


## If fasta file is unset or empty
if [ -z "${FASTA_INPUT}" ]; then 
    echo "ERROR: --fasta <input.fa> is required"
    BADARGS='TRUE'
elif [ "${N_fasta}" -gt 1 ]; then
    echo "WARNING: multiple --fasta provided. Only proceeding with the last one, ${FASTA_INPUT}"
fi


## If no mass spec input given
if [ -z "${MZML}" ] && [ -z "${RAW}" ] && [ -z "${DIA}" ]; then
    echo "ERROR: A mass spec input (--mzml --raw or --dai) is required"
    BADARGS='TRUE'
fi


## If more than one mass spec input is given
let N_SPECFILES=${N_mzml}+${N_raw}+${N_dia}
if [ "${N_SPECFILES}" -gt 1 ]; then
        echo 'ERROR: Multiple mass spec inputs provided. Provide ONE  of --mzml --dia or --raw'
        BADARGS='TRUE'
fi


## If output dir is not set, use current directory
if [ -z "${OUTPUT_DIR}" ]; then
    echo "INFO: --out not specified, output will be saved to current working directory"
    OUTPUT_DIR=${PWD}
fi


## If output dir is not writable
if [ ! -w "${OUTPUT_DIR}" ]; then
    echo "ERROR: output directory "${OUTPUT_DIR}" is not writable"
    exit 1
fi


## If any of the above checks failed
if [ "${BADARGS}" == 'TRUE' ]; then
    echo "Check arguments and try again."
    echo -e "command called:\n$(echo $0 $ARGS | sed -e 's/\( --\)/ \\\n   \1/g')"
    exit 1
else
    echo "INFO: arguments passed validation"
fi


# Check for singularity
if ! command -v singularity &> /dev/null; then
    echo "WARNING: singularity command not found, attempting to load singularity module"
    if ! command -v module &> /dev/null; then
        echo "WARNING: module command not found. Did you mean to run this on an HPC?"
        exit 1
    else
        if $(module avail singularity 2>&1 >/dev/null | grep -q 'No module'); then
            echo 'WARNING: module singularity not found'
            echo 'ERROR: singularity cannot be found. Recheck installation?'
            exit 1
        else
            echo 'INFO: module singularity found'
            module load singularity
        fi
    fi
else
    echo 'INFO: singularity command found'
fi


MASS_SPEC_INPUT="${MZML}${DIA}${RAW}"
echo -e "Spec Input:\t${MASS_SPEC_INPUT}"
echo -e "FASTA Input:\t${FASTA_INPUT}"
echo -e "Output to:\t${OUTPUT_DIR}/"
echo -e "Singularity:\t${SINGULARITY_IMAGE}"


singularity exec --cleanenv -H ${PWD} ${SINGULARITY_IMAGE} \
    diann \
    --f ${MASS_SPEC_INPUT}   \
    --lib  \
    --threads 24 \
    --verbose 1 \
    --out ${OUTPUT_DIR}/report.tsv \
    --qvalue 0.01 \
    --matrices \
    --out-lib ${OUTPUT_DIR}/report-lib.tsv \
    --gen-spec-lib \
    --predictor \
    --fasta ${FASTA_INPUT} \
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
    --no-ifs-removal 