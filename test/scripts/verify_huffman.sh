#!/bin/bash
# Taken from verify_s5.sh


if ! command -v realpath &> /dev/null
then
    realpath() {
        [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
    }
fi

# Get the path of this script
SCRIPTPATH=$(realpath  $(dirname "$0"))

ZSTD=""
REGIONS=""
SAMPLES=""
ZSTD_LEVEL=""
unset -v NO_KEEP

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    -f|--filename)
    FILENAME="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [ -z "${FILENAME}" ]
then
    echo "Specify a filename with --filename, -f <filename>"
    exit 1
fi

echo "FILENAME        = ${FILENAME}"

function exit_fail_rm_tmp {
    echo "Removing directory : ${TMPDIR}"
    rm -r ${TMPDIR}
    exit 1
}

"${SCRIPTPATH}"/../../gp_compression -f ${FILENAME} -t || { echo "Failed!"; exit_fail_rm_tmp; }
exit 0
