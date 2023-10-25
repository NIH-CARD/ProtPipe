#!/usr/bin/env bash
#SBATCH --mem 8G
#SBATCH --nodes 1
#SBATCH --time 2:00:00
#SBATCH --ntasks 2
#SBATCH --partition quick,norm




trap '[[ $? -eq 1 ]] && echo Halting execution due to errors' EXIT


# RAWFILE=$(head -n ${N} ${filelist} | tail -n 1)
# RAWFILE_BASENAME=$(basename $RAWFILE)
PWIZ='src/pwiz_sandbox'


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
            --file) assert_argument "$1" "$opt"; FILE="$1"; shift;;
            --dir) assert_argument "$1" "$opt"; DATADIR="$1"; shift;;
            --list) assert_argument "$1" "$opt"; LIST="$1"; shift;;
            --out) assert_argument "$1" "$opt"; OUT="$1"; shift;;
            
      
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


if [ ! -z "${FILE+x}" ] && [ ! -z "${DIR+x}" ]; then
    TOOMANY=True
fi

if [ -v "${FILE+x}" ] && [ -v "${LIST+x}" ]; then
    TOOMANY=True
fi

if [ -v "${DATADIR+x}" ] && [ -v "${LIST+x}" ]; then
    TOOMANY=True
fi

if [ -z "${DATADIR}" ] && [ -z "${FILE}" ] && [ -z "${LIST}" ]; then
    echo "ERROR: No inputs given!"
    echo "Provide ONE of --file, --dir or --list" 
    exit 1
fi

if [ "${TOOMANY}" == 'True' ]; then
    echo "ERROR: Multiple input modes provided."
    echo "Only provide ONE of --file, --dir or --list"
    exit 1
fi



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

#### CHECK PWIZ_SANDBOX ############################################################################
URL='https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2124988&authkey=AFC4xp1k-1JhHEk'


if [ -d src/pwiz_sandbox ]; then
    echo "INFO: ProteoWizard singularity sandbox wizard already configured at ${PWIZ}"
else
    echo "INFO: ProteoWizard singularity sandbox wizard not found at src/pwiz_sandbox, attempting download"
    wget -O pwiz_sandbox.tar.gz ${URL}
    tar -zxv -f pwiz_sandbox.tar.gz -C src/
fi



#### CONVERT TO MZML ###############################################################################


run_pwiz() {
    local datadir=${1}
    local outdir=${2}
    local infile=${3}
    echo "INFO: Converting ${infile} from ${datadir} to .mzML"
    echo "INFO: Output saved to ${outdir}"
    echo "INFO: /etc/localtime mount error can be ignored"

    echo "run_pwiz ${datadir} ${outdir} ${infile}"

    [[ $infile =~ .*\.(raw$|RAW$) ]] || { echo "$infile is not RAW file! Skipping..."; return 0; }

    singularity exec \
        -B ${datadir}:/data \
        -B ${outdir}:/mnt \
        -B `mktemp -d /dev/shm/wineXXX`:/mywineprefix \
        -w src/pwiz_sandbox \
        mywine msconvert \
            --32 \
            --filter "peakPicking vendor msLevel=1-" \
            -o /mnt/ --verbose \
            /data/${infile}
}
export -f run_pwiz

if [ -z "${OUT}" ]; then
    echo "INFO: Saving mzML to ./mzML"
    echo "      Override with --out DIRNAME"
    OUT="${PWD}/mzML"
fi

mkdir -p ${OUT}

if [ ! -z "${FILE}" ]; then
    DATADIR=$(dirname ${FILE})
    run_pwiz ${DATADIR} ${OUT} ${FILE}
elif [ ! -z "${DATADIR}" ]; then
    FILES=$(ls ${DATADIR})
    parallel -j 1 run_pwiz ::: ${DATADIR} ::: ${OUT} ::: ${FILES[@]}
elif [ ! -z "${LIST}" ]; then
    readarray -t FILES <${LIST}
    while read FILE; do
        DATADIR=$(dirname ${FILE})
        run_pwiz ${DATADIR} ${OUT} ${FILE}
    done
else
    echo "No input given!"
fi

echo "Reached the end... Exiting"
