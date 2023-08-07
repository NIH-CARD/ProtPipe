#!/usr/bin/env bash
#SBATCH --mem 8G
#SBATCH --nodes 1
#SBATCH --time 0:59:00
#SBATCH --ntasks 2
#SBATCH --partition quick,norm
N=${SLURM_ARRAY_TASK_ID}
filelist=${1}
RAWFILE=$(head -n ${N} ${filelist} | tail -n 1)
DATADIR="$(dirname ${RAWFILE})/"
RAWFILE_BASENAME=$(basename $RAWFILE)
PWIZ='src/pwiz_sandbox'


trap '[[ $? -eq 1 ]] && echo Halting execution due to errors' EXIT

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
            --dir) assert_argument "$1" "$opt"; DIR="$1"; shift;;
      
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
    wget -O pwiz.tar.gz ${URL}
    tar -zxv -f pwiz_sandbox.tar.gz -C src/
fi



#### CONVERT TO MZML ###############################################################################

echo "INFO: Converting ${RAWFILE} to .mzML"
echo "INFO: /etc/localtime mount error can be ignored"

singularity exec \
    -B ${DATADIR}:/data \
    -B `mktemp -d /dev/shm/wineXXX`:/mywineprefix \
    -w src/pwiz_sandbox \
    mywine msconvert \
        --32 \
        --filter "peakPicking vendor msLevel=1-" \
        -o /data/ --verbose \
        /data/${RAWFILE_BASENAME}
